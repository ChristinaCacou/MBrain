#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(edgeR))
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=3L) stop("Usage: Rscript run_discovery_edgeR.R counts.csv.gz matched_metadata.csv outdir")
out <- a[3]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
ct <- as.matrix(read.csv(gzfile(a[1]),row.names=1,check.names=FALSE))
md <- read.csv(a[2],stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(is.numeric(ct), all(is.finite(ct)), all(ct>=0), all(ct==round(ct)),
  !anyDuplicated(rownames(ct)), !anyDuplicated(colnames(ct)))
stopifnot(!anyDuplicated(md$sample),setequal(md$sample,colnames(ct)))
md <- md[match(colnames(ct),md$sample),,drop=FALSE]
stopifnot(all(md$diagnosis%in%c("CONT","OUD")),all(md$sex%in%c("male","female")))
results <- list()
for(region in c("DLPFC","NAC")) {
  m <- md[md$source_name_ch1==region,,drop=FALSE]
  stopifnot(nrow(m)==40L,all(table(factor(m$sex,levels=c("male","female")),
     factor(m$diagnosis,levels=c("CONT","OUD")))==10L))
  m$sex <- factor(m$sex,levels=c("male","female"))
  m$diagnosis <- factor(m$diagnosis,levels=c("CONT","OUD"))
  m$age <- as.numeric(m$age); m$pmi <- as.numeric(m$pmi)
  m$rpf_rin <- as.numeric(m$rpf_rin)
  stopifnot(all(is.finite(as.matrix(m[,c("age","pmi","rpf_rin")]))))
  m$age <- as.numeric(scale(m$age)); m$pmi <- as.numeric(scale(m$pmi))
  m$rpf_rin <- as.numeric(scale(m$rpf_rin))
  for(contrast in c("additive_OUD","sex_by_OUD")) {
    f <- if(contrast=="additive_OUD") ~ sex+diagnosis+age+pmi+rpf_rin
      else ~ sex*diagnosis+age+pmi+rpf_rin
    design <- model.matrix(f,data=m)
    stopifnot(qr(design)$rank==ncol(design),nrow(design)-ncol(design)>=10)
    y <- DGEList(counts=ct[,m$sample,drop=FALSE])
    keep <- filterByExpr(y,design=design)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y,method="TMM")
    y <- estimateDisp(y,design,robust=TRUE)
    fit <- glmQLFit(y,design,robust=TRUE)
    coefficient <- if(contrast=="additive_OUD") "diagnosisOUD"
      else "sexfemale:diagnosisOUD"
    stopifnot(coefficient%in%colnames(design))
    tst <- glmQLFTest(fit,coef=match(coefficient,colnames(design)))
    tab <- topTags(tst,n=Inf,sort.by="none")$table
    tab$gene_id <- rownames(tab)
    tab <- tab[,c("gene_id",setdiff(colnames(tab),"gene_id"))]
    stopifnot(all(is.finite(tab$PValue)),all(is.finite(tab$FDR)))
    write.csv(tab,file.path(out,paste(region,contrast,"edgeR_QL.csv",sep="_")),row.names=FALSE)
    egr <- tab[tab$gene_id=="ENSG00000120738",,drop=FALSE]
    results[[length(results)+1]] <- data.frame(
      region=region,contrast=contrast,genes_tested=nrow(tab),
      fdr05=sum(tab$FDR<0.05),egr1_tested=nrow(egr)==1L,
      egr1_logFC=if(nrow(egr)) egr$logFC else NA_real_,
      egr1_p=if(nrow(egr)) egr$PValue else NA_real_,
      egr1_q=if(nrow(egr)) egr$FDR else NA_real_,
      edgeR_version=as.character(packageVersion("edgeR")))
    png(file.path(out,paste(region,contrast,"BCV.png",sep="_")),width=900,height=650)
    plotBCV(y,main=paste(region,contrast,"edgeR BCV"));dev.off()
  }
}
summary <- do.call(rbind,results)
write.csv(summary,file.path(out,"edgeR_VERIFIED_SUMMARY.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"R_SESSION_INFO.txt"))
print(summary)
