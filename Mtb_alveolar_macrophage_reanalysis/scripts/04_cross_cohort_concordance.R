#!/usr/bin/env Rscript
# 04_cross_cohort_concordance.R
# Final cross-study comparison used for manuscript Figure 20 and summary tables.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly=TRUE)
arg <- function(flag,default=NULL){i<-match(flag,args);if(is.na(i))default else args[i+1]}
g189 <- arg("--gse189996","GSE189996/results/primary_DE_121")
g223 <- arg("--gse223863","GSE223863/results/gene_level")
p189 <- arg("--camera189","GSE189996/results/pathway_CAMERA")
p223 <- arg("--camera223","GSE223863/results/gene_level/09_pathway_CAMERA")
out  <- arg("--outdir","cross_dataset_results/reproducible")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

contrasts <- c("Mtb_vs_Control_2h","Mtb_vs_Control_24h","Mtb_vs_Control_72h",
               "Response_24h_vs_2h","Response_72h_vs_2h","Response_72h_vs_24h")

priority <- c("IL1B","SERPINB2","IFIT1","IFIT2","IFIT3","IFIT5","STAT1","STAT2","IRF1","IRF7",
              "IFI44L","ISG15","MX1","MX2","OAS2","OAS3","OASL","USP18","RSAD2","CMPK2",
              "HERC5","GBP1","GBP4","GBP5")

read_gene <- function(root,cn,study) {
  candidates <- c(file.path(root,paste0(cn,"_all_genes.tsv")),
                  file.path(root,"02_primary_DE",paste0(cn,"_all_genes.tsv")),
                  file.path(root,"03_time_interactions",paste0(cn,"_all_genes.tsv")))
  fp <- candidates[file.exists(candidates)][1]
  if(is.na(fp)) stop(study,": missing ",cn)
  z <- fread(fp)
  gcol <- intersect(c("gene_name","symbol","SYMBOL"),names(z))[1]
  if(is.na(gcol)) stop(study,": gene symbol column missing")
  setnames(z,gcol,"gene_name")
  # remove ambiguous duplicated symbols rather than arbitrarily choosing one
  dup <- z[!is.na(gene_name) & gene_name!="", .N, by=gene_name][N>1,gene_name]
  z <- z[!gene_name%in%dup & !is.na(gene_name) & gene_name!=""]
  z[,study:=study]
  z
}

gsum <- list(); pall <- list()
for(cn in contrasts) {
  a <- read_gene(g189,cn,"GSE189996")
  b <- read_gene(g223,cn,"GSE223863")
  m <- merge(a[,.(gene_name,logFC_189=logFC,FDR_189=adj.P.Val)],
             b[,.(gene_name,logFC_223=logFC,FDR_223=adj.P.Val)],by="gene_name")
  m[,same_direction:=sign(logFC_189)==sign(logFC_223)]
  m[,contrast:=cn]
  fwrite(m,file.path(out,paste0(cn,"_shared_genes.tsv")),sep="\t")

  sig_either <- m$FDR_189<0.05 | m$FDR_223<0.05
  sig_both   <- m$FDR_189<0.05 & m$FDR_223<0.05
  gsum[[cn]] <- data.table(
    contrast=cn,n_shared_genes=nrow(m),
    pearson_logFC=cor(m$logFC_189,m$logFC_223,method="pearson"),
    spearman_logFC=cor(m$logFC_189,m$logFC_223,method="spearman"),
    direction_concordance_all=mean(m$same_direction),
    n_FDR05_GSE189996=sum(m$FDR_189<0.05),
    n_FDR05_GSE223863=sum(m$FDR_223<0.05),
    n_FDR05_both=sum(sig_both),
    direction_concordance_FDR05_both=if(sum(sig_both)) mean(m$same_direction[sig_both]) else NA_real_
  )
  pp <- m[gene_name%in%priority]
  pall[[cn]] <- pp
}
fwrite(rbindlist(gsum),file.path(out,"gene_concordance_summary.tsv"),sep="\t")
fwrite(rbindlist(pall),file.path(out,"priority_gene_concordance.tsv"),sep="\t")

read_path <- function(root,cn) {
  cand <- c(file.path(root,paste0(cn,"_Hallmark_CAMERA.tsv")),
            file.path(root,paste0(cn,"_Hallmark_CAMERA_FDR05.tsv")),
            file.path(root,paste0(cn,"_Hallmark_cameraPR.tsv")))
  fp <- cand[file.exists(cand)][1]
  if(is.na(fp)) stop("Missing Hallmark: ",root," ",cn)
  z <- fread(fp)
  setnames(z,intersect(c("Pathway","pathway"),names(z))[1],"Pathway")
  if(!"FDR"%in%names(z)) z[,FDR:=p.adjust(PValue,"BH")]
  z
}

psum <- list(); allp <- list()
for(cn in contrasts) {
  a<-read_path(p189,cn); b<-read_path(p223,cn)
  m<-merge(a[,.(Pathway,Direction_189=Direction,FDR_189=FDR)],
           b[,.(Pathway,Direction_223=Direction,FDR_223=FDR)],by="Pathway")
  m[,same_direction:=Direction_189==Direction_223]
  m[,contrast:=cn]
  both<-m$FDR_189<0.05 & m$FDR_223<0.05
  either<-m$FDR_189<0.05 | m$FDR_223<0.05
  psum[[cn]]<-data.table(
    contrast=cn,n_shared_pathways=nrow(m),
    n_same_direction=sum(m$same_direction),
    direction_concordance=mean(m$same_direction),
    n_FDR05_189=sum(m$FDR_189<0.05),
    n_FDR05_223=sum(m$FDR_223<0.05),
    n_FDR05_both=sum(both),
    n_FDR05_either=sum(either),
    concordance_FDR05_either=if(sum(either))mean(m$same_direction[either]) else NA_real_,
    concordance_FDR05_both=if(sum(both))mean(m$same_direction[both]) else NA_real_)
  allp[[cn]]<-m
}
fwrite(rbindlist(psum),file.path(out,"Hallmark_concordance_summary.tsv"),sep="\t")
fwrite(rbindlist(allp),file.path(out,"Hallmark_all.tsv"),sep="\t")

cat("Cross-cohort concordance complete.\n")
