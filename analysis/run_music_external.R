# Full external adapter and batching controls. The official package is unchanged.
source("environment/music/runtime.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=4L,as.character(packageVersion("MuSiC"))=="1.0.0")
input<-jsonlite::fromJSON(args[1]);out<-args[2];triple<-args[3];mode<-args[4]
batch_mode<-if(length(args)>=5L)args[5] else "pooled"
stopifnot(mode %in% c("controls","full"),batch_mode %in% c("pooled","separate"))
dir.create(out,recursive=TRUE,showWarnings=FALSE);started<-proc.time()["elapsed"]
selection<-jsonlite::fromJSON(input$selection_manifest);types<-input$cell_types;data_dir<-input$data_directory
cells<-read.delim(file.path(data_dir,"cells.tsv"),check.names=FALSE,colClasses="character")
refs<-read.delim(input$reference_file,check.names=FALSE,colClasses="character")
cases<-read.delim(input$case_file,check.names=FALSE,colClasses="character")
targets<-read.delim(file.path(data_dir,"targets.tsv"),check.names=FALSE,colClasses="character")
genes<-readLines(file.path(data_dir,"genes.tsv"))
read_matrix<-function(path){h<-gzfile(path,"rb");on.exit(close(h));methods::as(Matrix::readMM(h),"CsparseMatrix")}
bulk_all<-read_matrix(file.path(data_dir,"targets.mtx.gz"));dimnames(bulk_all)<-list(genes,targets$target_name)
refs<-refs[refs$triple_key==triple,,drop=FALSE];stopifnot(nrow(refs)==90L)
cases<-cases[cases$triple_key==triple,,drop=FALSE];stopifnot(nrow(cases)==360L)
needed<-strsplit(refs$reference_donors[1],"|",fixed=TRUE)[[1]]
held<-unique(cases$held_out);stopifnot(length(held)==4L,!any(held %in% needed))
targets<-targets[targets$held_out %in% held,,drop=FALSE];stopifnot(nrow(targets)==240L,all(table(targets$held_out)==60L))
bulk_all<-as.matrix(bulk_all[,targets$target_name,drop=FALSE])
counts_by_donor<-lapply(needed,function(d){
  path<-input$donor_matrices$path[input$donor_matrices$donor==d];stopifnot(length(path)==1L)
  x<-read_matrix(path);cm<-cells[cells$donor==d,,drop=FALSE]
  stopifnot(nrow(x)==length(genes),ncol(x)==nrow(cm));dimnames(x)<-list(genes,cm$cell_id);x
});names(counts_by_donor)<-needed
.music_trace<-new.env(parent=emptyenv());.music_trace$items<-list()
invisible(trace("music.iter",where=asNamespace("MuSiC"),print=FALSE,exit=quote({
  res<-returnValue();st<-get(".music_trace",envir=.GlobalEnv)
  st$items[[length(st$items)+1L]]<-list(convergence=res$converge,n_features=nrow(D),
    variance_finite=all(is.finite(res$var.p)),nnls_coefficient_sum=sum(res$q.nnls),weighted_coefficient_sum=sum(res$q.weight))
})))
make_sce<-function(ref){
  ids<-as.integer(strsplit(ref$reference_columns_R,"|",fixed=TRUE)[[1]]);cm<-cells[ids,,drop=FALSE]
  quota<-table(cm$donor,cm$cell_type);expected<-ifelse(rownames(quota)==ref$dominant_donor,as.integer(ref$major_count),as.integer(ref$minor_count))
  stopifnot(length(ids)==as.integer(ref$budget)*6L,length(unique(ids))==length(ids),setequal(unique(cm$donor),needed),
    all(quota==expected),all(colSums(quota)==as.integer(ref$budget)),!any(held %in% cm$donor))
  pieces<-lapply(needed,function(d)counts_by_donor[[d]][,cm$cell_id[cm$donor==d],drop=FALSE])
  x<-do.call(cbind,pieces);x<-x[,cm$cell_id,drop=FALSE];rownames(cm)<-cm$cell_id
  stopifnot(max(abs(Matrix::colSums(x)-as.numeric(cm$total_counts)))==0)
  SingleCellExperiment(assays=list(counts=x),colData=S4Vectors::DataFrame(cm))
}
fit_once<-function(sce,tm){
  .music_trace$items<-list()
  fit<-music_prop(bulk.mtx=bulk_all[,tm$target_name,drop=FALSE],sc.sce=sce,markers=NULL,clusters="cell_type",samples="donor",
    select.ct=types,cell_size=NULL,ct.cov=FALSE,verbose=FALSE,iter.max=1000,nu=0.0001,eps=0.01,centered=FALSE,normalize=FALSE)
  diagnostics<-do.call(rbind,lapply(.music_trace$items,as.data.frame));stopifnot(nrow(diagnostics)==nrow(tm))
  diagnostics$target_name<-tm$target_name
  pred<-list(music_weighted=fit$Est.prop.weighted[tm$target_name,types,drop=FALSE],music_nnls=fit$Est.prop.allgene[tm$target_name,types,drop=FALSE])
  for(p in pred)stopifnot(all(is.finite(p)),min(p)>=-1e-10,max(abs(rowSums(p)-1))<1e-10)
  list(pred=pred,diag=diagnostics)
}
fit_reference<-function(sce,style){
  if(style=="pooled")return(fit_once(sce,targets))
  pieces<-lapply(unique(targets$held_out),function(d)fit_once(sce,targets[targets$held_out==d,,drop=FALSE]))
  pred<-lapply(names(pieces[[1]]$pred),function(m)do.call(rbind,lapply(pieces,function(x)x$pred[[m]]))[targets$target_name,,drop=FALSE])
  names(pred)<-names(pieces[[1]]$pred)
  diag<-do.call(rbind,lapply(pieces,function(x)x$diag));diag<-diag[match(targets$target_name,diag$target_name),,drop=FALSE]
  list(pred=pred,diag=diag)
}
format_fit<-function(ref,fit){
  cfg<-cases[cases$reference_id==ref$reference_id,,drop=FALSE];stopifnot(nrow(cfg)==4L)
  cm<-cfg[match(targets$held_out,cfg$held_out),,drop=FALSE];stopifnot(!anyNA(cm$case_id))
  truth<-as.matrix(data.frame(lapply(targets[paste0("true_",0:5)],as.numeric)));colnames(truth)<-types
  stopifnot(max(abs(rowSums(truth)-1))<1e-12)
  rows<-lapply(names(fit$pred),function(m){
    pred<-fit$pred[[m]][targets$target_name,types,drop=FALSE];err<-pred-truth
    tab<-data.frame(case_id=cm$case_id,reference_id=ref$reference_id,held_out=targets$held_out,triple_id=as.integer(cm$triple_id),
      triple_key=triple,block=as.integer(ref$block),budget=as.integer(ref$budget),level=ref$level,dominant_donor=ref$dominant_donor,
      target_name=targets$target_name,mixture_id=as.integer(targets$mixture_id),method=m,mae=rowMeans(abs(err)),rmse=sqrt(rowMeans(err^2)),check.names=FALSE)
    for(k in 1:6){tab[[paste0("true_",k-1L)]]<-truth[,k];tab[[paste0("pred_",k-1L)]]<-pred[,k]};tab
  })
  diag<-fit$diag[match(targets$target_name,fit$diag$target_name),,drop=FALSE]
  diag$reference_id<-ref$reference_id;diag$case_id<-cm$case_id;diag$held_out<-targets$held_out;diag$mixture_id<-as.integer(targets$mixture_id)
  list(rows=do.call(rbind,rows),diag=diag)
}
append_csv<-function(x,path)write.table(x,file=path,sep=",",row.names=FALSE,col.names=!file.exists(path),append=file.exists(path),na="NA")
if(mode=="controls"){
  first<-selection$donors[1];dominant<-selection$donors[2]
  testcases<-cases[cases$held_out==first & cases$triple_id=="0" & cases$block=="0",,drop=FALSE]
  chosen<-testcases[(testcases$budget=="60" & testcases$level=="balanced") |
    (testcases$budget=="120" & testcases$level=="ratio4" & testcases$dominant_donor==dominant) |
    (testcases$budget=="300" & testcases$level=="ratio10" & testcases$dominant_donor==dominant),,drop=FALSE]
  stopifnot(nrow(chosen)==3L);checks<-list()
  for(i in seq_len(nrow(chosen))){
    ref<-refs[refs$reference_id==chosen$reference_id[i],,drop=FALSE];sce<-make_sce(ref)
    pooled<-fit_reference(sce,"pooled");separate<-fit_reference(sce,"separate")
    difference<-max(sapply(names(pooled$pred),function(m)max(abs(pooled$pred[[m]]-separate$pred[[m]]))))
    dc<-c("convergence","n_features","variance_finite","target_name")
    diag_equal<-all(sapply(dc,function(k)identical(pooled$diag[[k]],separate$diag[[k]])))
    coef_diff<-max(abs(as.matrix(pooled$diag[c("nnls_coefficient_sum","weighted_coefficient_sum")])-as.matrix(separate$diag[c("nnls_coefficient_sum","weighted_coefficient_sum")])))
    for(style in c("pooled","separate")){
      tab<-format_fit(ref,if(style=="pooled")pooled else separate)
      append_csv(tab$rows,file.path(out,paste0(style,"_predictions.csv")));append_csv(tab$diag,file.path(out,paste0(style,"_diagnostics.csv")))
    }
    checks[[i]]<-list(reference_id=ref$reference_id,prediction_maximum_difference=difference,
      diagnostic_discrete_equal=diag_equal,diagnostic_coefficient_maximum_difference=coef_diff,
      equivalence_passed=difference<=1e-10 && diag_equal && coef_diff<=1e-10)
    cat("Batch control",i,"/ 3; max difference",difference,"\n")
  }
  jsonlite::write_json(list(status="controls completed",checks=checks,batch_equivalence_passed=all(sapply(checks,function(x)x$equivalence_passed)),
    prediction_rows_per_path=1440L,elapsed_seconds=unname(proc.time()["elapsed"]-started)),file.path(out,"controls.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
}else{
  total<-0L;limits<-0L;badvar<-0L;support<-c(Inf,-Inf)
  for(i in seq_len(nrow(refs))){
    ref<-refs[i,,drop=FALSE];fit<-fit_reference(make_sce(ref),batch_mode);tab<-format_fit(ref,fit)
    stopifnot(nrow(tab$rows)==480L,nrow(tab$diag)==240L)
    append_csv(tab$rows,file.path(out,"predictions.csv"));append_csv(tab$diag,file.path(out,"diagnostics.csv"))
    total<-total+nrow(tab$rows);limits<-limits+sum(tab$diag$convergence=="Reach Maxiter");badvar<-badvar+sum(!tab$diag$variance_finite)
    support<-c(min(support[1],tab$diag$n_features),max(support[2],tab$diag$n_features))
    if(i%%5L==0L || i==nrow(refs)){
      progress<-list(triple_key=triple,references_completed=i,total_references=nrow(refs),prediction_rows=total,elapsed_seconds=unname(proc.time()["elapsed"]-started))
      jsonlite::write_json(progress,file.path(out,"progress.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
      cat(triple,i,"/",nrow(refs),"references; elapsed",round(progress$elapsed_seconds,1),"s\n");flush.console()
    }
  }
  stopifnot(total==43200L)
  jsonlite::write_json(list(status="completed",triple_key=triple,batch_mode=batch_mode,reference_matrices=nrow(refs),logical_cases=nrow(cases),
    prediction_rows=total,weighted_fits=total/2L,maxiter_count=limits,nonfinite_variance_count=badvar,effective_gene_range=support,
    elapsed_seconds=unname(proc.time()["elapsed"]-started)),file.path(out,"fold.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
}
