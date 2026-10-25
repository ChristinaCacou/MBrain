#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(edgeR))
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=3L) stop('Usage: Rscript donor22_edgeR_QL.R raw_counts.csv.gz matched_metadata.csv outdir')
dir.create(a[3],recursive=TRUE,showWarnings=FALSE)
ct <- as.matrix(read.csv(gzfile(a[1]),row.names=1,check.names=FALSE))
md <- read.csv(a[2],stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(ct)==60676L,ncol(ct)==80L,all(is.finite(ct)),all(ct>=0),all(ct==round(ct)),!anyDuplicated(rownames(ct)),!anyDuplicated(colnames(ct)),!anyDuplicated(md$sample),setequal(md$sample,colnames(ct)))
md <- md[match(colnames(ct),md$sample),,drop=FALSE]
stopifnot(md$diagnosis[md$sample=='HN22']=='CONT',md$sex[md$sample=='HN22']=='male',md$diagnosis[md$sample=='HBA22']=='CONT',md$sex[md$sample=='HBA22']=='male')
output <- list()
for(region in c('DLPFC','NAC')) {
  m <- md[md$source_name_ch1==region,,drop=FALSE]
  stopifnot(nrow(m)==40L)
  exclude <- if(region=='DLPFC') 'HBA22' else 'HN22'
  m <- m[m$sample!=exclude,,drop=FALSE]
  stopifnot(nrow(m)==39L)
  m$sex <- factor(m$sex,levels=c('male','female'))
  m$diagnosis <- factor(m$diagnosis,levels=c('CONT','OUD'))
  for(v in c('age','pmi','rpf_rin')) {
    m[[v]] <- as.numeric(m[[v]])
    stopifnot(all(is.finite(m[[v]])),sd(m[[v]])>0)
    m[[v]] <- as.numeric(scale(m[[v]]))
  }
  for(contrast in c('additive_OUD','sex_by_OUD')) {
    f <- if(contrast=='additive_OUD') ~ sex+diagnosis+age+pmi+rpf_rin else ~ sex*diagnosis+age+pmi+rpf_rin
    design <- model.matrix(f,data=m)
    stopifnot(qr(design)$rank==ncol(design),nrow(design)-ncol(design)>=10)
    y <- DGEList(counts=ct[,m$sample,drop=FALSE])
    keep <- filterByExpr(y,design=design)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y,method='TMM')
    y <- estimateDisp(y,design,robust=TRUE)
    fit <- glmQLFit(y,design,robust=TRUE)
    coefficient <- if(contrast=='additive_OUD') 'diagnosisOUD' else 'sexfemale:diagnosisOUD'
    stopifnot(coefficient %in% colnames(design))
    tst <- glmQLFTest(fit,coef=match(coefficient,colnames(design)))
    tab <- topTags(tst,n=Inf,sort.by='none')$table
    tab$gene_id <- rownames(tab)
    tab <- tab[,c('gene_id',setdiff(colnames(tab),'gene_id'))]
    stopifnot(all(is.finite(tab$PValue)),all(is.finite(tab$FDR)))
    write.csv(tab,file.path(a[3],paste(region,contrast,'donor22_deleted_edgeR_QL.csv',sep='_')),row.names=FALSE)
    for(g in c('ENSG00000120738','ENSG00000269028','ENSG00000196126')) {
      e <- tab[tab$gene_id==g,,drop=FALSE]
      output[[length(output)+1L]] <- data.frame(region=region,contrast=contrast,excluded=exclude,n=39L,genes_tested=nrow(tab),number_q05=sum(tab$FDR<.05),gene_id=g,logFC=if(nrow(e)==1L)e$logFC else NA_real_,PValue=if(nrow(e)==1L)e$PValue else NA_real_,FDR=if(nrow(e)==1L)e$FDR else NA_real_,edgeR_version=as.character(packageVersion('edgeR')))
    }
  }
}
s <- do.call(rbind,output)
write.csv(s,file.path(a[3],'donor22_deleted_selected_gene_summary.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(a[3],'R_SESSION_INFO.txt'))
print(s)
