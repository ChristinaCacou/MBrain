#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(edgeR))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4L) stop("usage: Rscript edgeR_EGR1_all_donors.R counts.csv.gz metadata.csv outdir batch_number_1_to_8")
raw<-as.matrix(read.csv(gzfile(args[1]),row.names=1,check.names=FALSE))
md<-read.csv(args[2],stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(raw)==60676L,ncol(raw)==80L, all(is.finite(raw)),all(raw>=0),all(raw==round(raw)),
  !anyDuplicated(colnames(raw)),!anyDuplicated(rownames(raw)),!anyDuplicated(md$sample),setequal(colnames(raw),md$sample))
md<-md[match(colnames(raw),md$sample),,drop=FALSE]
m<-md[md$source_name_ch1=="DLPFC",,drop=FALSE]
stopifnot(nrow(m)==40L, all(sort(m$sample)==sort(paste0("HBA",1:40))),
          all(table(m$sex,m$diagnosis)==10L))
batch<-as.integer(args[4]);stopifnot(!is.na(batch),batch>=1L,batch<=8L)
dir.create(args[3],recursive=TRUE,showWarnings=FALSE)
target_ids<-((batch-1L)*5L+1L):(batch*5L)
allout<-list()
for(id in target_ids) {
  omit<-paste0("HBA",id)
  mm<-m[m$sample!=omit,,drop=FALSE]
  mm$sex<-factor(mm$sex,levels=c("male","female"))
  mm$diagnosis<-factor(mm$diagnosis,levels=c("CONT","OUD"))
  for(v in c("age","pmi","rpf_rin")) {
    mm[[v]]<-as.numeric(mm[[v]])
    stopifnot(all(is.finite(mm[[v]])),sd(mm[[v]])>0)
    mm[[v]]<-as.numeric(scale(mm[[v]]))
  }
  x<-model.matrix(~sex+diagnosis+age+pmi+rpf_rin,data=mm)
  stopifnot(nrow(x)==39L,qr(x)$rank==ncol(x), "diagnosisOUD" %in% colnames(x))
  y<-DGEList(counts=raw[,mm$sample,drop=FALSE])
  keep<-filterByExpr(y,design=x)
  y<-y[keep,,keep.lib.sizes=FALSE]
  y<-calcNormFactors(y,method="TMM")
  y<-estimateDisp(y,x,robust=TRUE)
  fit<-glmQLFit(y,x,robust=TRUE)
  test<-glmQLFTest(fit,coef=match("diagnosisOUD",colnames(x)))
  t<-topTags(test,n=Inf,sort.by="none")$table
  t$gene_id<-rownames(t)
  stopifnot(isTRUE(all.equal(as.numeric(t$FDR),as.numeric(p.adjust(t$PValue,method="BH")),tolerance=1e-10)))
  e<-t[t$gene_id=="ENSG00000120738",,drop=FALSE]
  stopifnot(nrow(e)==1L)
  allout[[length(allout)+1L]]<-data.frame(
    omitted_sample=omit,donor_id=id,deleted_sex=m$sex[match(omit,m$sample)],
    deleted_diagnosis=m$diagnosis[match(omit,m$sample)],n=39L,
    tested_genes=nrow(t),fdr_lt_05=sum(t$FDR<.05),
    EGR1_logFC=e$logFC,EGR1_logCPM=e$logCPM,
    EGR1_F=e$F,EGR1_PValue=e$PValue,EGR1_FDR=e$FDR,
    edgeR_version=as.character(packageVersion("edgeR")),stringsAsFactors=FALSE)
  t<-t[,c("gene_id","logFC","PValue","FDR")]
  write.csv(t,gzfile(file.path(args[3],paste0("DLPFC_omit_",id,"_full_gene_table.csv.gz"))),row.names=FALSE)
  cat("Done omit",id,"EGR1 FC",e$logFC,"p",e$PValue,"q",e$FDR,"tested",nrow(t),"\n")
}
z<-do.call(rbind,allout)
write.csv(z,file.path(args[3],paste0("LOO_EGR1_batch_",batch,"_summary.csv")),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[3],paste0("sessionInfo_batch_",batch,".txt")))
stopifnot(nrow(z)==5L,all(z$EGR1_logFC<0))
print(z)
