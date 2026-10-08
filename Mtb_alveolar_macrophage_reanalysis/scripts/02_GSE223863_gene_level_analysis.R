#!/usr/bin/env Rscript
# 02_GSE223863_gene_level_analysis.R
# Final donor-blocked gene-level replication workflow

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly=TRUE)
arg <- function(flag, default=NULL) {
  i <- match(flag,args); if (is.na(i)) return(default); args[i+1]
}
counts_file <- arg("--counts")
meta_file   <- arg("--metadata")
ann_file    <- arg("--annotation", NULL)
outdir      <- arg("--outdir","reproducible_results/GSE223863_gene_level")
if (is.null(counts_file) || is.null(meta_file))
  stop("Usage: Rscript 02_GSE223863_gene_level_analysis.R --counts counts.rds --metadata metadata.tsv [--annotation annotation.tsv]")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

counts <- if (grepl("\\.rds$",counts_file,ignore.case=TRUE)) readRDS(counts_file) else as.matrix(fread(counts_file),rownames=1)
if (inherits(counts,"DGEList")) counts <- counts$counts
counts <- as.matrix(counts)
meta <- fread(meta_file)

pick <- function(cand) { z<-cand[cand%in%names(meta)]; if(!length(z)) stop(paste(cand,collapse="/")); z[1] }
setnames(meta,
 c(pick(c("SRR","run_accession","sample_id","sample")),
   pick(c("donor","Donor")),
   pick(c("time","time_point","Time")),
   pick(c("condition","group","sample_group","Condition"))),
 c("sample","donor","time","condition"))

meta[,time:=gsub("hr","h",as.character(time),ignore.case=TRUE)]
meta[,time:=factor(time,levels=c("2h","24h","72h"))]
meta[,condition:=factor(fifelse(grepl("mtb|infect",condition,ignore.case=TRUE),"Mtb","Control"),
                        levels=c("Control","Mtb"))]
meta[,donor:=factor(donor)]

stopifnot(all(meta$sample %in% colnames(counts)))
counts <- counts[,meta$sample,drop=FALSE]

group <- interaction(meta$time,meta$condition,drop=TRUE)
y <- DGEList(counts=counts)
keep <- filterByExpr(y,group=group)
y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
cat("Genes after filter:",nrow(y),"\n")

design <- model.matrix(~ time*condition,data=meta)
v0 <- voom(y,design,plot=FALSE)
dc <- duplicateCorrelation(v0,design,block=meta$donor)
cat("Consensus donor correlation:",dc$consensus.correlation,"\n")
v <- voom(y,design,block=meta$donor,correlation=dc$consensus.correlation,plot=FALSE)
fit0 <- lmFit(v,design,block=meta$donor,correlation=dc$consensus.correlation)

findcoef <- function(p) { z<-grep(p,colnames(design),value=TRUE); if(length(z)!=1) stop(p); z }
m <- findcoef("^conditionMtb$")
i24 <- findcoef("^time24h:conditionMtb$|^conditionMtb:time24h$")
i72 <- findcoef("^time72h:conditionMtb$|^conditionMtb:time72h$")

C <- matrix(0,ncol(design),6,dimnames=list(colnames(design),
 c("Mtb_vs_Control_2h","Mtb_vs_Control_24h","Mtb_vs_Control_72h",
   "Response_24h_vs_2h","Response_72h_vs_2h","Response_72h_vs_24h")))
C[m,1]<-1
C[m,2]<-1; C[i24,2]<-1
C[m,3]<-1; C[i72,3]<-1
C[i24,4]<-1
C[i72,5]<-1
C[i72,6]<-1; C[i24,6]<--1

fit <- eBayes(contrasts.fit(fit0,C),robust=TRUE)
saveRDS(list(y=y,v=v,design=design,correlation=dc$consensus.correlation,contrasts=C,fit=fit),
        file.path(outdir,"GSE223863_gene_level_model.rds"))

ann <- NULL
if (!is.null(ann_file) && file.exists(ann_file)) ann <- fread(ann_file)

summ <- list()
for (cn in colnames(C)) {
  tt <- topTable(fit,coef=cn,n=Inf,sort.by="P")
  tt$gene_id <- rownames(tt)
  if (!is.null(ann)) {
    gid <- intersect(c("gene_id","Geneid","ENSEMBL"),names(ann))[1]
    gnm <- intersect(c("gene_name","symbol","SYMBOL"),names(ann))[1]
    if(!is.na(gid)&&!is.na(gnm)) tt$gene_name <- ann[[gnm]][match(tt$gene_id,ann[[gid]])]
  }
  fwrite(tt,file.path(outdir,paste0(cn,"_all_genes.tsv")),sep="\t")
  fwrite(tt[tt$adj.P.Val<0.05,],file.path(outdir,paste0(cn,"_FDR05.tsv")),sep="\t")
  summ[[cn]] <- data.table(contrast=cn,tested=nrow(tt),FDR05=sum(tt$adj.P.Val<0.05))
}
fwrite(rbindlist(summ),file.path(outdir,"DE_summary_all_contrasts.tsv"),sep="\t")

# Focus genes used across the project
focus <- c("IFIT1","IFIT2","IFIT3","IFIT5","STAT1","STAT2","IRF1","IRF7","IRF9","IFI44L",
           "ISG15","MX1","MX2","OAS2","OAS3","OASL","USP18","RSAD2","CMPK2","HERC5",
           "GBP1","GBP4","GBP5","SERPINB2","IL1B","IDO1")
for (cn in colnames(C)) {
  f <- fread(file.path(outdir,paste0(cn,"_all_genes.tsv")))
  if("gene_name"%in%names(f)) fwrite(f[gene_name%in%focus],file.path(outdir,paste0(cn,"_focus_genes.tsv")),sep="\t")
}

# Optional local MSigDB CAMERA/cameraPR pathway analysis
read_gmt <- function(path) {
  x<-readLines(path); z<-lapply(x,function(s){f<-strsplit(s,"\t",fixed=TRUE)[[1]]; list(n=f[1],g=f[-c(1,2)])})
  setNames(lapply(z,`[[`,"g"),vapply(z,`[[`,"","n"))
}
if (!is.null(ann)) {
  gid <- intersect(c("gene_id","Geneid","ENSEMBL"),names(ann))[1]
  gnm <- intersect(c("gene_name","symbol","SYMBOL"),names(ann))[1]
  syms <- ann[[gnm]][match(rownames(y),ann[[gid]])]
  gmt_root <- Sys.getenv("MSIGDB_DIR",unset="reference/MSigDB")
  for (spec in list(Hallmark=file.path(gmt_root,"h.all.v2025.1.Hs.symbols.gmt"),
                    Reactome=file.path(gmt_root,"c2.cp.reactome.v2025.1.Hs.symbols.gmt"))) {
    nm<-names(spec); fp<-spec[[1]]
    if(file.exists(fp)) {
      gs<-read_gmt(fp)
      ix<-lapply(gs,function(g) which(syms%in%g)); ix<-ix[lengths(ix)>=5]
      for(j in seq_len(ncol(C))) {
        stat <- fit$t[,j]; names(stat)<-rownames(y)
        cp <- cameraPR(stat,index=ix,inter.gene.cor=0.01,use.ranks=FALSE)
        cp$Pathway<-rownames(cp); cp$FDR<-p.adjust(cp$PValue,"BH")
        fwrite(cp,file.path(outdir,paste0(colnames(C)[j],"_",nm,"_cameraPR.tsv")),sep="\t")
      }
    }
  }
}

cat("GSE223863 gene-level analysis complete.\n")
