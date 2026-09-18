# Read-only diagnostics from the unchanged official MuSiC implementation.
# Run from project root: Rscript --vanilla this.R selection.json output triple control|full
# This adapter never alters the official package or replaces a fitted component.
# A read-only exit trace records convergence; controls check it against untraced calls.
source("environment/music/runtime.R")
options(digits=17)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4L, as.character(packageVersion("MuSiC"))=="1.0.0")
manifest <- jsonlite::fromJSON(args[1]); out <- args[2]; triple <- args[3]; mode <- args[4]
stopifnot(mode %in% c("control", "full"))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
started <- proc.time()["elapsed"]
input <- jsonlite::fromJSON(manifest$input_manifest)
types <- input$cell_types; data_dir <- input$data_directory
cells <- read.delim(file.path(data_dir,"cells.tsv"), check.names=FALSE, colClasses="character")
refs <- read.delim(manifest$references_file, check.names=FALSE, colClasses="character")
cases <- read.delim(manifest$cases_file, check.names=FALSE, colClasses="character")
targets <- read.delim(manifest$targets_file, check.names=FALSE, colClasses="character")
all_targets <- read.delim(file.path(data_dir,"targets.tsv"), check.names=FALSE, colClasses="character")
inventory <- read.csv(manifest$inventory_file, check.names=FALSE, colClasses="character")
genes <- readLines(file.path(data_dir,"genes.tsv"))
stopifnot(length(genes)==30172L, !anyDuplicated(genes), length(types)==6L,
          nrow(targets)==112L, all(table(targets$held_out)==8L))
read_matrix <- function(path) {
  h <- gzfile(path,"rb"); on.exit(close(h))
  methods::as(Matrix::readMM(h),"CsparseMatrix")
}
bulk_all <- read_matrix(file.path(data_dir,"targets.mtx.gz"))
dimnames(bulk_all) <- list(genes,all_targets$target_name)
refs <- refs[refs$triple_key==triple,,drop=FALSE]
cases <- cases[cases$triple_key==triple,,drop=FALSE]
stopifnot(nrow(refs)==36L, nrow(cases)==144L, !anyDuplicated(refs$reference_id),
          setequal(refs$level,c("balanced","ratio10")))
needed <- strsplit(refs$reference_donors[1],"|",fixed=TRUE)[[1]]
held <- unique(cases$held_out)
stopifnot(length(needed)==3L, length(held)==4L, !any(held %in% needed),
          all(refs$reference_donors==refs$reference_donors[1]))
targets <- targets[targets$held_out %in% held,,drop=FALSE]
stopifnot(nrow(targets)==32L, !anyDuplicated(targets$target_name),
          all(table(targets$held_out)==8L), all(targets$target_name %in% colnames(bulk_all)))
# Preserve the target order in the frozen selection; do not construct new targets.
bulk_all <- as.matrix(bulk_all[,targets$target_name,drop=FALSE])
counts_by_donor <- lapply(needed,function(d) {
  path <- input$donor_matrices$path[input$donor_matrices$donor==d]
  stopifnot(length(path)==1L)
  x <- read_matrix(path); cm <- cells[cells$donor==d,,drop=FALSE]
  stopifnot(nrow(x)==length(genes), ncol(x)==nrow(cm))
  dimnames(x) <- list(genes,cm$cell_id); x
}); names(counts_by_donor) <- needed
refs <- refs[order(as.integer(refs$block),as.integer(refs$budget),
                   refs$level!="balanced",refs$dominant_donor),,drop=FALSE]
profile_cache <- new.env(parent=emptyenv())
.mechanism_trace <- new.env(parent=emptyenv());.mechanism_trace$items <- character()
.trace_active <- FALSE
set_trace <- function(enabled) {
  if(enabled && !.trace_active) {
    invisible(trace("music.iter",where=asNamespace("MuSiC"),print=FALSE,exit=quote({
      ans <- returnValue();st <- get(".mechanism_trace",envir=.GlobalEnv)
      st$items <- c(st$items,ans$converge)
    })))
    .trace_active <<- TRUE
  } else if(!enabled && .trace_active) {
    invisible(untrace("music.iter",where=asNamespace("MuSiC")))
    .trace_active <<- FALSE
  }
}
basis_maxdiff <- c(theta=0, sigma=0, cell_size=0, mean_cell_size=0, design=0)
type_inventory <- setNames(c("B","cM","T4","T8","ncM","NK"),types)
profile_key <- function(d,k,b,n) paste(d,k,b,n,sep="\t")
append_csv <- function(x,path) {
  # Explicit 17 significant digits avoid write.table's default numeric rounding.
  for(k in names(x)) if(is.double(x[[k]])) {
    keep_na <- is.na(x[[k]]); x[[k]] <- sprintf("%.17g",x[[k]]); x[[k]][keep_na] <- NA_character_
  }
  write.table(x,file=path,sep=",",row.names=FALSE,col.names=!file.exists(path),
              append=file.exists(path),na="NA",quote=TRUE)
}
expected_outputs <- if(mode=="control") c("control_predictions.csv","controls.json",
  "control_weights.csv.gz","basis.csv") else c("predictions.csv","basis.csv",
  "weights.csv","weight_pairs.csv","profile_dispersion.csv","fold.json")
stopifnot(!any(file.exists(file.path(out,expected_outputs))))
make_sce <- function(ref) {
  ids <- as.integer(strsplit(ref$reference_columns_R,"|",fixed=TRUE)[[1]])
  stopifnot(all(ids>=1L),all(ids<=nrow(cells)))
  cm <- cells[ids,,drop=FALSE]
  quota <- table(cm$donor,cm$cell_type)
  expected <- ifelse(rownames(quota)==ref$dominant_donor,
                     as.integer(ref$major_count),as.integer(ref$minor_count))
  stopifnot(length(ids)==as.integer(ref$budget)*6L,length(unique(ids))==length(ids),
            setequal(unique(cm$donor),needed),all(quota==expected),
            all(colSums(quota)==as.integer(ref$budget)),!any(held %in% cm$donor))
  pieces <- lapply(needed,function(d) counts_by_donor[[d]][,cm$cell_id[cm$donor==d],drop=FALSE])
  x <- do.call(cbind,pieces); x <- x[,cm$cell_id,drop=FALSE]; rownames(cm) <- cm$cell_id
  stopifnot(max(abs(Matrix::colSums(x)-as.numeric(cm$total_counts)))==0)
  SingleCellExperiment(assays=list(counts=x),colData=S4Vectors::DataFrame(cm))
}
cache_profiles <- function(sce,ref) {
  cm <- as.data.frame(SummarizedExperiment::colData(sce))
  for(d in needed) for(k in types) {
    ids <- cm$cell_id[cm$donor==d & cm$cell_type==k]; n <- length(ids)
    key <- profile_key(d,k,ref$block,n)
    if(exists(key,envir=profile_cache,inherits=FALSE)) {
      stopifnot(identical(get(key,envir=profile_cache)$ids,ids))
    } else {
      # All 30,172 genes remain in these profiles, including zeros.
      totals <- Matrix::rowSums(counts_by_donor[[d]][,ids,drop=FALSE])
      total <- sum(totals); stopifnot(total>0)
      theta <- totals/total; s <- total/n
      stopifnot(all(is.finite(theta)),abs(sum(theta)-1)<1e-12,is.finite(s),s>0)
      assign(key,list(ids=ids,theta=theta,S=s,donor=d,cell_type=k,
                      block=as.integer(ref$block),n=n),envir=profile_cache)
    }
  }
}
check_basis <- function(sce,ref) {
  cache_profiles(sce,ref)
  official <- music_basis(sce,non.zero=TRUE,markers=NULL,clusters="cell_type",
                          samples="donor",select.ct=types,cell_size=NULL,ct.cov=FALSE,verbose=FALSE)
  support <- rownames(official$M.theta)
  stopifnot(!anyDuplicated(support),all(support %in% genes),
            all(is.finite(official$M.theta)),all(is.finite(official$Sigma)),
            all(is.finite(official$S)),all(official$Sigma>=0))
  rows <- list()
  for(k in types) {
    profiles <- lapply(needed,function(d) {
      n <- if(d==ref$dominant_donor) as.integer(ref$major_count) else as.integer(ref$minor_count)
      get(profile_key(d,k,ref$block,n),envir=profile_cache)
    }); names(profiles) <- needed
    theta <- do.call(cbind,lapply(profiles,function(p)p$theta))
    direct_mean <- rowMeans(theta)
    direct_var <- rowSums((theta-direct_mean)^2)/(ncol(theta)-1L)
    direct_s <- vapply(profiles,function(p)p$S,numeric(1))
    excluded <- setdiff(genes,support)
    stopifnot(all(direct_mean[excluded]==0),all(direct_var[excluded]==0))
    delta <- c(theta=max(abs(direct_mean[support]-official$M.theta[,k])),
               sigma=max(abs(direct_var[support]-official$Sigma[,k])),
               cell_size=max(abs(direct_s[rownames(official$S)]-official$S[,k])),
               mean_cell_size=abs(mean(direct_s)-official$M.S[k]),
               design=max(abs(direct_mean[support]*mean(direct_s)-official$Disgn.mtx[,k])))
    # Unname scalar elements to prevent compound names such as mean_cell_size.B cells.
    names(delta) <- names(basis_maxdiff)
    stopifnot(all(is.finite(delta)),all(delta<=1e-10))
    basis_maxdiff <<- pmax(basis_maxdiff,delta)
    rows[[length(rows)+1L]] <- data.frame(reference_id=ref$reference_id,triple_key=triple,
      block=as.integer(ref$block),budget=as.integer(ref$budget),level=ref$level,
      dominant_donor=ref$dominant_donor,cell_type=k,sum_sigma=sum(official$Sigma[,k]),
      theta_squared_norm=sum(official$M.theta[,k]^2),mean_cell_size=unname(official$M.S[k]),
      basis_genes=length(support),fixed_profile_genes=length(genes),check.names=FALSE)
  }
  append_csv(do.call(rbind,rows),file.path(out,"basis.csv"))
}
fit_once <- function(sce,tm,traced=TRUE) {
  set_trace(traced);.mechanism_trace$items <- character()
  fit <- music_prop(bulk.mtx=bulk_all[,tm$target_name,drop=FALSE],sc.sce=sce,
    markers=NULL,clusters="cell_type",samples="donor",select.ct=types,cell_size=NULL,
    ct.cov=FALSE,verbose=FALSE,iter.max=1000,nu=0.0001,eps=0.01,centered=FALSE,normalize=FALSE)
  pred <- list(music_weighted=fit$Est.prop.weighted[tm$target_name,types,drop=FALSE],
               music_nnls=fit$Est.prop.allgene[tm$target_name,types,drop=FALSE])
  for(p in pred) stopifnot(all(is.finite(p)),min(p)>=-1e-10,max(abs(rowSums(p)-1))<1e-10)
  # Weight.gene contains NA for target-excluded genes. Never turn those NAs into zero weights.
  w <- fit$Weight.gene[,tm$target_name,drop=FALSE]
  stopifnot(!anyDuplicated(rownames(w)),all(rownames(w) %in% genes),
            all(is.na(w) | (is.finite(w) & w>0)),all(colSums(!is.na(w))>1L),
            all(is.finite(fit$Var.prop)),all(is.finite(fit$r.squared.full)))
  convergence <- if(traced) .mechanism_trace$items else NULL
  if(traced) stopifnot(length(convergence)==nrow(tm),
                       !any(convergence=="Reach Maxiter"),all(grepl("^Converge at [0-9]+$",convergence)))
  if(traced) names(convergence) <- tm$target_name
  list(pred=pred,weights=w,convergence=convergence)
}
positive_weights <- function(w,target) {
  x <- w[,target]; names(x) <- rownames(w)
  x[is.finite(x)&x>0]
}
format_fit <- function(ref,fit,tm=targets) {
  cfg <- cases[cases$reference_id==ref$reference_id,,drop=FALSE]; stopifnot(nrow(cfg)==4L)
  cm <- cfg[match(tm$held_out,cfg$held_out),,drop=FALSE]; stopifnot(!anyNA(cm$case_id))
  truth <- as.matrix(data.frame(lapply(tm[paste0("true_",0:5)],as.numeric))); colnames(truth) <- types
  stopifnot(max(abs(rowSums(truth)-1))<1e-12)
  do.call(rbind,lapply(names(fit$pred),function(m) {
    pred <- fit$pred[[m]][tm$target_name,types,drop=FALSE]; err <- pred-truth
    tab <- data.frame(case_id=cm$case_id,reference_id=ref$reference_id,held_out=tm$held_out,
      triple_id=as.integer(cm$triple_id),triple_key=triple,block=as.integer(ref$block),
      budget=as.integer(ref$budget),level=ref$level,dominant_donor=ref$dominant_donor,
      target_name=tm$target_name,mixture_id=as.integer(tm$mixture_id),method=m,
      mae=rowMeans(abs(err)),rmse=sqrt(rowMeans(err^2)),check.names=FALSE)
    for(k in 1:6) {tab[[paste0("true_",k-1L)]] <- truth[,k];tab[[paste0("pred_",k-1L)]] <- pred[,k]}
    tab
  }))
}
weight_summary <- function(ref,w) {
  do.call(rbind,lapply(seq_len(nrow(targets)),function(i) {
    v <- positive_weights(w,targets$target_name[i]); p <- v/sum(v); n <- length(v)
    data.frame(reference_id=ref$reference_id,triple_key=triple,block=as.integer(ref$block),
      budget=as.integer(ref$budget),level=ref$level,dominant_donor=ref$dominant_donor,
      held_out=targets$held_out[i],target_name=targets$target_name[i],
      mixture_id=as.integer(targets$mixture_id[i]),n_weight_genes=n,
      normalized_weight_effective_fraction=(1/sum(p^2))/n,
      top1pct_weight_share=sum(sort(p,decreasing=TRUE)[seq_len(ceiling(0.01*n))]),check.names=FALSE)
  }))
}
paired_weights <- function(ref,w,balanced,balanced_id) {
  do.call(rbind,lapply(seq_len(nrow(targets)),function(i) {
    a <- positive_weights(balanced,targets$target_name[i]); b <- positive_weights(w,targets$target_name[i])
    common <- intersect(names(a),names(b)); union <- union(names(a),names(b))
    stopifnot(length(common)>1L)
    aa <- a[common];bb <- b[common]
    rho <- if(sd(aa)==0 || sd(bb)==0) NA_real_ else cor(aa,bb,method="spearman")
    data.frame(reference_id=ref$reference_id,balanced_reference_id=balanced_id,triple_key=triple,
      block=as.integer(ref$block),budget=as.integer(ref$budget),level=ref$level,
      dominant_donor=ref$dominant_donor,held_out=targets$held_out[i],
      target_name=targets$target_name[i],mixture_id=as.integer(targets$mixture_id[i]),
      intersection_n=length(common),union_n=length(union),Jaccard=length(common)/length(union),
      spearman_shared=rho,normalized_weight_tv_intersection=0.5*sum(abs(aa/sum(aa)-bb/sum(bb))),
      balanced_weight_fraction_on_intersection=sum(aa)/sum(a),
      unequal_weight_fraction_on_intersection=sum(bb)/sum(b),
      spearman_defined=is.finite(rho),check.names=FALSE)
  }))
}
profile_dispersion <- function() {
  rows <- list(); ns <- c(5L,10L,20L,25L,40L,50L,100L,250L)
  for(d in needed) for(k in types) for(n in ns) {
    ps <- lapply(0:2,function(b)get(profile_key(d,k,b,n),envir=profile_cache))
    mat <- do.call(cbind,lapply(ps,function(p)p$theta)); avg <- rowMeans(mat)
    ss <- vapply(ps,function(p)p$S,numeric(1))
    pair_tv <- vapply(list(c(1L,2L),c(1L,3L),c(2L,3L)),
                       function(pair)0.5*sum(abs(mat[,pair[1]]-mat[,pair[2]])),numeric(1))
    N <- as.integer(inventory[inventory$donor_id==d,type_inventory[[k]]])
    stopifnot(length(N)==1L,N>=n,all(is.finite(mat)),all(is.finite(ss)))
    rows[[length(rows)+1L]] <- data.frame(donor=d,cell_type=k,reference_cells=n,
      source_inventory_N=N,sampling_fraction=n/N,draws=3L,fixed_profile_genes=length(genes),
      theta_variance_trace=sum((mat-avg)^2)/(ncol(mat)-1L),
      theta_pairwise_tv=mean(pair_tv),S_cv=sd(ss)/mean(ss),check.names=FALSE)
  }
  result <- do.call(rbind,rows);stopifnot(nrow(result)==144L)
  append_csv(result,file.path(out,"profile_dispersion.csv"));nrow(result)
}
if(mode=="control") {
  dominant <- sort(needed)[1]
  chosen <- refs[refs$block=="0" & (refs$level=="balanced" |
                   (refs$level=="ratio10" & refs$dominant_donor==dominant)),,drop=FALSE]
  stopifnot(nrow(chosen)==6L)
  checks <- list();raw <- list()
  for(i in seq_len(nrow(chosen))) {
    ref <- chosen[i,,drop=FALSE];sce <- make_sce(ref);check_basis(sce,ref)
    pooled <- fit_once(sce,targets)
    plain <- fit_once(sce,targets,traced=FALSE)
    plain_pred_diff <- max(vapply(names(pooled$pred),function(m)
      max(abs(pooled$pred[[m]]-plain$pred[[m]])),numeric(1)))
    plain_weight_diff <- 0;plain_support_equal <- TRUE
    for(t in targets$target_name) {
      a <- positive_weights(pooled$weights,t);b <- positive_weights(plain$weights,t)
      plain_support_equal <- plain_support_equal && setequal(names(a),names(b))
      if(!setequal(names(a),names(b))) stop("Trace control target weight support mismatch")
      plain_weight_diff <- max(plain_weight_diff,max(abs(a-b[names(a)])))
    }
    plain_passed <- plain_pred_diff<=1e-10 && plain_weight_diff<=1e-10 && plain_support_equal
    stopifnot(plain_passed);rm(plain)
    tab <- format_fit(ref,pooled);tab$batch_mode <- "pooled"
    append_csv(tab,file.path(out,"control_predictions.csv"))
    pred_diff <- 0;weight_diff <- 0;support_equal <- TRUE;convergence_equal <- TRUE
    for(d in unique(targets$held_out)) {
      tm <- targets[targets$held_out==d,,drop=FALSE];separate <- fit_once(sce,tm)
      tab <- format_fit(ref,separate,tm);tab$batch_mode <- "separate"
      append_csv(tab,file.path(out,"control_predictions.csv"))
      pred_diff <- max(pred_diff,vapply(names(pooled$pred),function(m)
        max(abs(pooled$pred[[m]][tm$target_name,,drop=FALSE]-separate$pred[[m]])),numeric(1)))
      convergence_equal <- convergence_equal && identical(pooled$convergence[tm$target_name],separate$convergence)
      for(t in tm$target_name) {
        a <- positive_weights(pooled$weights,t);b <- positive_weights(separate$weights,t)
        support_equal <- support_equal && setequal(names(a),names(b))
        if(!setequal(names(a),names(b))) stop("Control target weight support mismatch")
        weight_diff <- max(weight_diff,max(abs(a-b[names(a)])))
      }
      # The first frozen target per held-out donor, both raw paths, is retained for audit.
      t <- tm$target_name[1]
      for(style in c("pooled","separate")) {
        v <- positive_weights(if(style=="pooled")pooled$weights else separate$weights,t)
        raw[[length(raw)+1L]] <- data.frame(reference_id=ref$reference_id,held_out=d,
          target_name=t,batch_mode=style,gene=names(v),weight=as.numeric(v),check.names=FALSE)
      }
    }
    passed <- pred_diff<=1e-10 && weight_diff<=1e-10 && support_equal && convergence_equal && plain_passed
    checks[[i]] <- list(reference_id=ref$reference_id,prediction_maximum_difference=pred_diff,
      weight_maximum_difference=weight_diff,positive_weight_support_equal=support_equal,
      convergence_equal=convergence_equal,convergence_status_counts=as.list(table(pooled$convergence)),
      untraced_prediction_maximum_difference=plain_pred_diff,untraced_weight_maximum_difference=plain_weight_diff,
      untraced_positive_weight_support_equal=plain_support_equal,untraced_equivalence_passed=plain_passed,
      equivalence_passed=passed)
    stopifnot(passed)
    cat("Mechanism control",i,"/ 6; prediction maxdiff",pred_diff,"weight maxdiff",weight_diff,"\n");flush.console()
  }
  raw <- do.call(rbind,raw)
  # Gzip avoids a large readable weight dump; no matrix of all production weights is kept.
  h <- gzfile(file.path(out,"control_weights.csv.gz"),"wt")
  raw$weight <- sprintf("%.17g",raw$weight)
  write.table(raw,h,sep=",",row.names=FALSE,quote=TRUE);close(h)
  jsonlite::write_json(list(status="controls completed",triple_key=triple,checks=checks,
    batch_equivalence_passed=all(vapply(checks,function(x)x$equivalence_passed,logical(1))),
    trace_equivalence_passed=all(vapply(checks,function(x)x$untraced_equivalence_passed,logical(1))),
    prediction_rows=768L,prediction_rows_per_path=384L,representative_raw_weight_rows=nrow(raw),
    additional_untraced_weighted_control_fits=192L,traced_fits=384L,maxiter_count=0L,
    basis_checks_maximum_difference=as.list(basis_maxdiff),basis_checks_passed=all(basis_maxdiff<=1e-10),
    official_algorithm_unchanged=TRUE,tracing_used=TRUE,trace_scope="read-only convergence exit value",
    elapsed_seconds=unname(proc.time()["elapsed"]-started)),file.path(out,"controls.json"),
    auto_unbox=TRUE,pretty=TRUE,digits=17)
} else {
  total <- 0L;weight_total <- 0L;pair_total <- 0L;undefined_rho <- 0L;convergence_all <- character()
  balanced <- NULL;balanced_id <- NULL;balanced_key <- NULL;support <- c(Inf,-Inf)
  for(i in seq_len(nrow(refs))) {
    ref <- refs[i,,drop=FALSE];sce <- make_sce(ref);check_basis(sce,ref)
    fit <- fit_once(sce,targets);tab <- format_fit(ref,fit)
    convergence_all <- c(convergence_all,unname(fit$convergence))
    stopifnot(nrow(tab)==64L)
    append_csv(tab,file.path(out,"predictions.csv"));total <- total+nrow(tab)
    weights <- weight_summary(ref,fit$weights);stopifnot(nrow(weights)==32L)
    append_csv(weights,file.path(out,"weights.csv"));weight_total <- weight_total+nrow(weights)
    support <- c(min(support[1],weights$n_weight_genes),max(support[2],weights$n_weight_genes))
    key <- paste(ref$block,ref$budget,sep="|")
    if(ref$level=="balanced") {
      balanced <- fit$weights;balanced_key <- key;balanced_id <- ref$reference_id
    } else {
      stopifnot(identical(key,balanced_key),!is.null(balanced))
      pairs <- paired_weights(ref,fit$weights,balanced,balanced_id);stopifnot(nrow(pairs)==32L)
      append_csv(pairs,file.path(out,"weight_pairs.csv"));pair_total <- pair_total+nrow(pairs)
      undefined_rho <- undefined_rho+sum(!pairs$spearman_defined)
    }
    if(i%%3L==0L || i==nrow(refs)) {
      progress <- list(triple_key=triple,references_completed=i,total_references=nrow(refs),
        prediction_rows=total,elapsed_seconds=unname(proc.time()["elapsed"]-started))
      jsonlite::write_json(progress,file.path(out,"progress.json"),auto_unbox=TRUE,pretty=TRUE,digits=17)
      cat(triple,i,"/",nrow(refs),"mechanism references; elapsed",round(progress$elapsed_seconds,1),"s\n");flush.console()
    }
    rm(fit,sce)
  }
  dispersion_rows <- profile_dispersion()
  stopifnot(total==2304L,weight_total==1152L,pair_total==864L)
  jsonlite::write_json(list(status="completed",triple_key=triple,reference_matrices=nrow(refs),
    logical_cases=nrow(cases),prediction_rows=total,weighted_fits=total/2L,basis_rows=6L*nrow(refs),
    traced_fits=length(convergence_all),
    weight_rows=weight_total,weight_pair_rows=pair_total,profile_dispersion_rows=dispersion_rows,
    convergence_status_counts=as.list(table(convergence_all)),maxiter_count=sum(convergence_all=="Reach Maxiter"),
    cached_donor_type_draw_profiles=length(ls(profile_cache)),effective_gene_range=support,
    undefined_spearman_rows=undefined_rho,finite_predictions_and_observed_weights=TRUE,
    basis_checks_maximum_difference=as.list(basis_maxdiff),basis_checks_passed=all(basis_maxdiff<=1e-10),
    official_algorithm_unchanged=TRUE,tracing_used=TRUE,trace_scope="read-only convergence exit value",
    elapsed_seconds=unname(proc.time()["elapsed"]-started)),file.path(out,"fold.json"),
    auto_unbox=TRUE,pretty=TRUE,digits=17)
}
set_trace(FALSE)
