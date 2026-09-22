# Technical gate only. Official MuSiC and scientific protocol are unchanged.
source("environment/music/runtime.R")
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L, as.character(packageVersion("MuSiC"))=="1.0.0")
input <- jsonlite::fromJSON(args[1]); out <- args[2]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
started <- proc.time()["elapsed"]
selection <- jsonlite::fromJSON(input$selection_manifest)
types <- input$cell_types
read_matrix <- function(path) {
  handle <- gzfile(path,"rb"); on.exit(close(handle))
  methods::as(Matrix::readMM(handle),"CsparseMatrix")
}
official <- function(sce,bulk) {
  music_prop(bulk.mtx=bulk,sc.sce=sce,markers=NULL,clusters="cell_type",samples="donor",
    select.ct=types,cell_size=NULL,ct.cov=FALSE,verbose=FALSE,
    iter.max=1000,nu=0.0001,eps=0.01,centered=FALSE,normalize=FALSE)
}
.music_trace <- new.env(parent=emptyenv())
run_control <- function(case_id,ref_counts,ref_meta,bulk,targets) {
  stopifnot(identical(colnames(ref_counts),ref_meta$cell_id),length(unique(ref_meta$donor))==3L,
    !any(targets$held_out %in% ref_meta$donor), ncol(bulk)==3L)
  rownames(ref_meta)<-ref_meta$cell_id
  sce<-SingleCellExperiment(assays=list(counts=ref_counts),colData=S4Vectors::DataFrame(ref_meta))
  .music_trace$items<-list()
  invisible(trace("music.iter",where=asNamespace("MuSiC"),print=FALSE,exit=quote({
    result<-returnValue(); st<-get(".music_trace",envir=.GlobalEnv)
    st$items[[length(st$items)+1L]]<-list(convergence=result$converge,n_features=nrow(D),
      variance_finite=all(is.finite(result$var.p)),weighted_coefficient_sum=sum(result$q.weight),
      nnls_coefficient_sum=sum(result$q.nnls))
  })))
  fit<-official(sce,bulk); diagnostics<-.music_trace$items
  invisible(untrace("music.iter",where=asNamespace("MuSiC")))
  plain<-official(sce,bulk)
  difference<-max(abs(fit$Est.prop.weighted-plain$Est.prop.weighted),abs(fit$Est.prop.allgene-plain$Est.prop.allgene))
  stopifnot(difference==0,length(diagnostics)==3L)
  truth<-as.matrix(data.frame(lapply(targets[paste0("true_",0:5)],as.numeric)));colnames(truth)<-types
  stopifnot(max(abs(rowSums(truth)-1))<1e-12)
  outputs<-list(music_weighted=fit$Est.prop.weighted,music_nnls=fit$Est.prop.allgene)
  rows<-lapply(names(outputs),function(method) {
    prediction<-outputs[[method]][targets$target_name,types,drop=FALSE]
    stopifnot(all(is.finite(prediction)),min(prediction)>=-1e-10,max(abs(rowSums(prediction)-1))<1e-10)
    errors<-prediction-truth
    tab<-data.frame(case_id=case_id,target_name=targets$target_name,held_out=targets$held_out,
      mixture_id=as.integer(targets$mixture_id),method=method,mae=rowMeans(abs(errors)),rmse=sqrt(rowMeans(errors^2)),check.names=FALSE)
    for (k in 1:6) {tab[[paste0("true_",k-1L)]]<-truth[,k];tab[[paste0("pred_",k-1L)]]<-prediction[,k]}
    tab
  })
  predictions<-do.call(rbind,rows);write.csv(predictions,file.path(out,paste0(case_id,".csv")),row.names=FALSE)
  diag<-do.call(rbind,lapply(seq_along(diagnostics),function(i)data.frame(case_id=case_id,target_name=targets$target_name[i],as.data.frame(diagnostics[[i]]))))
  write.csv(diag,file.path(out,paste0(case_id,"_diagnostics.csv")),row.names=FALSE)
  list(case_id=case_id,prediction_rows=6L,reference_cells=ncol(ref_counts),targets=3L,
    traced_vs_plain_maximum_difference=difference,maximum_row_sum_error=max(sapply(outputs,function(p)max(abs(rowSums(p)-1)))),
    maxiter_count=sum(diag$convergence=="Reach Maxiter"),nonfinite_variance_count=sum(!diag$variance_finite),support_range=range(diag$n_features))
}
old<-jsonlite::fromJSON("results/music_input/20260917T090438982999Z/input_manifest.json")
olddata<-old$data_directory; oldgenes<-readLines(file.path(olddata,"genes.tsv"))
oldcells<-read.delim(file.path(olddata,"cells.tsv"),check.names=FALSE,colClasses="character")
oldcases<-read.delim(file.path(olddata,"cases.tsv"),check.names=FALSE,colClasses="character")
oldtargets<-read.delim(file.path(olddata,"targets.tsv"),check.names=FALSE,colClasses="character")
oldcounts<-read_matrix(file.path(olddata,"counts.mtx.gz"));dimnames(oldcounts)<-list(oldgenes,oldcells$cell_id)
oldbulk<-read_matrix(file.path(olddata,"targets.mtx.gz"));dimnames(oldbulk)<-list(oldgenes,oldtargets$target_name)
oc<-oldcases[oldcases$case_id=="101_t00_b0_balanced",,drop=FALSE]
oid<-as.integer(strsplit(oc$reference_columns_R,"|",fixed=TRUE)[[1]])
ot<-oldtargets[oldtargets$held_out=="101" & as.integer(oldtargets$mixture_id)<3L,,drop=FALSE]
controls<-list(run_control("old_replay",oldcounts[,oid,drop=FALSE],oldcells[oid,,drop=FALSE],as.matrix(oldbulk[,ot$target_name,drop=FALSE]),ot))
rm(oldcounts,oldbulk);gc()
data_dir<-input$data_directory
cells<-read.delim(file.path(data_dir,"cells.tsv"),check.names=FALSE,colClasses="character")
refs<-read.delim(input$reference_file,check.names=FALSE,colClasses="character")
cases<-read.delim(input$case_file,check.names=FALSE,colClasses="character")
genes<-readLines(file.path(data_dir,"genes.tsv"))
targets<-read.delim(file.path(data_dir,"targets.tsv"),check.names=FALSE,colClasses="character")
bulk<-read_matrix(file.path(data_dir,"targets.mtx.gz"));dimnames(bulk)<-list(genes,targets$target_name)
held<-selection$donors[1];dominant<-selection$donors[2] # first forward-cyclic reference donor
targets<-targets[targets$held_out==held & as.integer(targets$mixture_id)<3L,,drop=FALSE]
bulk<-as.matrix(bulk[,targets$target_name,drop=FALSE])
base<-cases[cases$held_out==held & cases$triple_id=="0" & cases$block=="0",,drop=FALSE]
chosen<-base[(base$budget=="60" & base$level=="balanced") |
    (base$budget=="120" & base$level=="ratio4" & base$dominant_donor==dominant) |
    (base$budget=="300" & base$level=="ratio10" & base$dominant_donor==dominant),,drop=FALSE]
chosen<-chosen[order(as.integer(chosen$budget)),,drop=FALSE];stopifnot(nrow(chosen)==3L)
needed<-strsplit(refs$reference_donors[match(chosen$reference_id[1],refs$reference_id)],"|",fixed=TRUE)[[1]]
counts_by_donor<-lapply(needed,function(d){
  path<-input$donor_matrices$path[input$donor_matrices$donor==d];stopifnot(length(path)==1L)
  x<-read_matrix(path);cm<-cells[cells$donor==d,,drop=FALSE]
  stopifnot(nrow(x)==length(genes),ncol(x)==nrow(cm));dimnames(x)<-list(genes,cm$cell_id);x
});names(counts_by_donor)<-needed
for (i in seq_len(nrow(chosen))) {
  c<-chosen[i,,drop=FALSE];ref<-refs[refs$reference_id==c$reference_id,,drop=FALSE]
  ids<-as.integer(strsplit(ref$reference_columns_R,"|",fixed=TRUE)[[1]]);rmdata<-cells[ids,,drop=FALSE]
  pieces<-lapply(needed,function(d)counts_by_donor[[d]][,rmdata$cell_id[rmdata$donor==d],drop=FALSE])
  counts<-do.call(cbind,pieces);counts<-counts[,rmdata$cell_id,drop=FALSE]
  quota<-table(rmdata$donor,rmdata$cell_type);expected<-ifelse(rownames(quota)==ref$dominant_donor,as.integer(ref$major_count),as.integer(ref$minor_count))
  stopifnot(length(ids)==6L*as.integer(c$budget),length(unique(ids))==length(ids),all(quota==expected),all(colSums(quota)==as.integer(c$budget)))
  controls[[length(controls)+1L]]<-run_control(c$case_id,counts,rmdata,bulk,targets)
  cat("External technical control",i,"/ 3 completed\n")
}
jsonlite::write_json(list(status="technical controls completed",controls=controls,external_prediction_rows=18L,
    old_replay_prediction_rows=6L,scientific_grid_prediction_rows=0L,held_out=held,dominant_donor=dominant,
    elapsed_seconds=unname(proc.time()["elapsed"]-started)),file.path(out,"technical.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
