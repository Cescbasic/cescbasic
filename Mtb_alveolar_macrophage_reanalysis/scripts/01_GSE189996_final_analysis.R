#!/usr/bin/env Rscript
# 01_GSE189996_final_analysis.R
# Final locked GSE189996 analysis used in the manuscript
# R 4.4.1; edgeR 4.4.2; limma 3.62.2

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

BASE <- Sys.getenv("BASE", unset = getwd())
OUT  <- file.path(BASE, "reproducible_results", "GSE189996")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(BASE, "GSE189996/results/primary_121/GSE189996_STAR_raw_counts_primary_121.rds")
meta_file   <- file.path(BASE, "GSE189996/results/primary_121/GSE189996_metadata_primary_121.tsv")
ann_file    <- file.path(BASE, "GSE189996/results/primary_121/GSE189996_GENCODEv19_annotation_primary.rds")

stopifnot(file.exists(counts_file), file.exists(meta_file), file.exists(ann_file))

counts <- readRDS(counts_file)
ann    <- readRDS(ann_file)
meta   <- fread(meta_file)

if (inherits(counts, "DGEList")) counts <- counts$counts
counts <- as.matrix(counts)

pick_col <- function(x, candidates) {
  z <- candidates[candidates %in% names(x)]
  if (!length(z)) stop("Missing required metadata column: ", paste(candidates, collapse=", "))
  z[1]
}
sample_col <- pick_col(meta, c("SRR","run_accession","sample_id","sample"))
donor_col  <- pick_col(meta, c("donor","Donor"))
time_col   <- pick_col(meta, c("time","time_point","Time"))
cond_col   <- pick_col(meta, c("condition","sample_group","group","Condition"))
sex_col    <- pick_col(meta, c("sex","Sex"))

setnames(meta, c(sample_col,donor_col,time_col,cond_col,sex_col),
               c("sample","donor","time","condition","sex"))

meta[, time := gsub("hr","h",as.character(time), ignore.case=TRUE)]
meta[, time := sub("^2$","2h",time)]
meta[, time := sub("^24$","24h",time)]
meta[, time := sub("^72$","72h",time)]
meta[, condition := fifelse(grepl("mtb|infect", condition, ignore.case=TRUE), "Mtb", "Control")]
meta[, sex := fifelse(grepl("^m", sex, ignore.case=TRUE), "Male", "Female")]

meta[, donor := factor(donor)]
meta[, time := factor(time, levels=c("2h","24h","72h"))]
meta[, condition := factor(condition, levels=c("Control","Mtb"))]
meta[, sex := factor(sex, levels=c("Female","Male"))]

stopifnot(all(meta$sample %in% colnames(counts)))
counts <- counts[, meta$sample, drop=FALSE]
stopifnot(identical(colnames(counts), meta$sample))

# -------- Filtering and TMM --------
group <- interaction(meta$time, meta$condition, drop=TRUE)
y0 <- DGEList(counts=counts)
keep <- filterByExpr(y0, group=group)
y <- y0[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")

cat("Input genes:", nrow(y0), "\n")
cat("Retained genes:", nrow(y), "\n")
stopifnot(nrow(y) == 10879)

saveRDS(y, file.path(OUT, "GSE189996_DGEList_filtered_TMM.rds"))
logCPM <- cpm(y, log=TRUE, prior.count=0.5)
saveRDS(logCPM, file.path(OUT, "GSE189996_logCPM.rds"))

# -------- Primary donor-aware model --------
design <- model.matrix(~ donor + time * condition, data=meta)
stopifnot(qr(design)$rank == ncol(design))

coef_name <- function(pattern) {
  hit <- grep(pattern, colnames(design), value=TRUE)
  if (length(hit) != 1L) stop("Could not uniquely match coefficient: ", pattern,
                              "\nAvailable: ", paste(colnames(design),collapse=", "))
  hit
}
c_mtb <- coef_name("^conditionMtb$")
c_24  <- coef_name("^time24h:conditionMtb$|^conditionMtb:time24h$")
c_72  <- coef_name("^time72h:conditionMtb$|^conditionMtb:time72h$")

C <- matrix(0, nrow=ncol(design), ncol=6,
            dimnames=list(colnames(design),
              c("Mtb_vs_Control_2h","Mtb_vs_Control_24h","Mtb_vs_Control_72h",
                "Response_24h_vs_2h","Response_72h_vs_2h","Response_72h_vs_24h")))
C[c_mtb,1] <- 1
C[c_mtb,2] <- 1; C[c_24,2] <- 1
C[c_mtb,3] <- 1; C[c_72,3] <- 1
C[c_24,4]  <- 1
C[c_72,5]  <- 1
C[c_72,6]  <- 1; C[c_24,6] <- -1

pdf(file.path(OUT,"voomWithQualityWeights_mean_variance.pdf"))
vqw <- voomWithQualityWeights(y, design, plot=TRUE)
dev.off()

fit0 <- lmFit(vqw, design)
fit  <- contrasts.fit(fit0, C)
fit  <- eBayes(fit, robust=TRUE)
saveRDS(list(design=design, contrasts=C, fit=fit, voom=vqw),
        file.path(OUT,"GSE189996_primary_model.rds"))

summ <- list()
for (cn in colnames(C)) {
  tt <- topTable(fit, coef=cn, number=Inf, sort.by="P")
  tt$gene_id <- rownames(tt)
  if (is.data.frame(ann) || is.data.table(ann)) {
    a <- as.data.table(ann)
    gid <- intersect(c("gene_id","Geneid","ENSEMBL"), names(a))[1]
    gnm <- intersect(c("gene_name","symbol","SYMBOL"), names(a))[1]
    if (!is.na(gid) && !is.na(gnm)) {
      m <- match(tt$gene_id, a[[gid]])
      tt$gene_name <- a[[gnm]][m]
    }
  }
  fwrite(tt, file.path(OUT,paste0(cn,"_all_genes.tsv")), sep="\t")
  sig <- tt[tt$adj.P.Val < 0.05,]
  fwrite(sig, file.path(OUT,paste0(cn,"_FDR05.tsv")), sep="\t")
  sig1 <- sig[abs(sig$logFC) >= 1,]
  fwrite(sig1, file.path(OUT,paste0(cn,"_FDR05_LFC1.tsv")), sep="\t")
  summ[[cn]] <- data.table(
    contrast=cn, tested=nrow(tt), FDR05=nrow(sig),
    Up=sum(sig$logFC>0), Down=sum(sig$logFC<0),
    FDR05_LFC1=nrow(sig1))
}
fwrite(rbindlist(summ), file.path(OUT,"GSE189996_primary_DE_summary.tsv"), sep="\t")

# Expected locked result counts
expected <- c(Mtb_vs_Control_2h=25, Mtb_vs_Control_24h=2604, Mtb_vs_Control_72h=4043,
              Response_24h_vs_2h=1046, Response_72h_vs_2h=2430, Response_72h_vs_24h=1559)
obs <- setNames(sapply(summ, function(z) z$FDR05), names(summ))
if (!all(obs[names(expected)] == expected)) {
  warning("DE counts differ from the locked CHPC run. Check package versions/input objects.\n",
          paste(names(expected), "expected", expected, "observed", obs[names(expected)], collapse="; "))
}

# -------- Priority genes --------
priority <- c("IFIT1","IFIT2","IFIT3","IFIT5","STAT1","STAT2","IRF1","IRF7","IRF9",
              "IFI44L","ISG15","MX1","MX2","OAS2","OAS3","OASL","USP18","RSAD2","CMPK2",
              "HERC5","GBP1","GBP4","GBP5","SERPINB2","IL1B","IDO1")
priority_out <- list()
for (cn in colnames(C)) {
  tt <- fread(file.path(OUT,paste0(cn,"_all_genes.tsv")))
  if ("gene_name" %in% names(tt))
    priority_out[[cn]] <- tt[gene_name %in% priority][, contrast:=cn]
}
if (length(priority_out))
  fwrite(rbindlist(priority_out,fill=TRUE), file.path(OUT,"GSE189996_priority_gene_primary_DE.tsv"), sep="\t")

# -------- Sensitivity A: same 121, standard voom --------
v_std <- voom(y, design, plot=FALSE)
fit_std <- eBayes(contrasts.fit(lmFit(v_std,design), C), robust=TRUE)

# -------- Sensitivity B: complete 17 donors, quality weights --------
full <- meta[, .(n=.N, nt=uniqueN(time), nc=uniqueN(condition)), by=donor][n==6 & nt==3 & nc==2, as.character(donor)]
stopifnot(length(full) == 17)
idx17 <- meta$donor %in% full
m17 <- droplevels(meta[idx17])
y17 <- y[,idx17,keep.lib.sizes=FALSE]
y17 <- calcNormFactors(y17)
d17 <- model.matrix(~ donor + time*condition, data=m17)
stopifnot(qr(d17)$rank == ncol(d17))
cn17 <- colnames(d17)
find17 <- function(p) { z<-grep(p,cn17,value=TRUE); if(length(z)!=1) stop(p); z }
m0<-find17("^conditionMtb$"); i24<-find17("^time24h:conditionMtb$|^conditionMtb:time24h$"); i72<-find17("^time72h:conditionMtb$|^conditionMtb:time72h$")
C17 <- matrix(0,ncol(d17),6,dimnames=list(cn17,colnames(C)))
C17[m0,1]<-1; C17[m0,2]<-1; C17[i24,2]<-1; C17[m0,3]<-1; C17[i72,3]<-1
C17[i24,4]<-1; C17[i72,5]<-1; C17[i72,6]<-1; C17[i24,6]<--1
v17 <- voomWithQualityWeights(y17,d17,plot=FALSE)
fit17 <- eBayes(contrasts.fit(lmFit(v17,d17),C17),robust=TRUE)

sens <- list()
for (cn in colnames(C)) {
  p <- topTable(fit,coef=cn,n=Inf,sort.by="none")
  a <- topTable(fit_std,coef=cn,n=Inf,sort.by="none")
  b <- topTable(fit17,coef=cn,n=Inf,sort.by="none")
  commonA <- intersect(rownames(p),rownames(a)); commonB <- intersect(rownames(p),rownames(b))
  psig <- rownames(p)[p$adj.P.Val<0.05]
  sens[[cn]] <- data.table(
    contrast=cn,
    primary_FDR05=sum(p$adj.P.Val<0.05),
    noQW_FDR05=sum(a$adj.P.Val<0.05),
    complete17_FDR05=sum(b$adj.P.Val<0.05),
    pearson_noQW=cor(p[commonA,"logFC"],a[commonA,"logFC"],method="pearson"),
    spearman_noQW=cor(p[commonA,"logFC"],a[commonA,"logFC"],method="spearman"),
    pearson_complete17=cor(p[commonB,"logFC"],b[commonB,"logFC"],method="pearson"),
    spearman_complete17=cor(p[commonB,"logFC"],b[commonB,"logFC"],method="spearman"),
    direction_primary_sig_noQW=mean(sign(p[psig,"logFC"])==sign(a[psig,"logFC"])),
    direction_primary_sig_complete17=mean(sign(p[psig,"logFC"])==sign(b[psig,"logFC"]))
  )
}
fwrite(rbindlist(sens), file.path(OUT,"GSE189996_DE_sensitivity_summary.tsv"), sep="\t")

# -------- GMT / CAMERA --------
read_gmt <- function(path) {
  x <- readLines(path)
  out <- lapply(x, function(line) {
    f <- strsplit(line,"\t",fixed=TRUE)[[1]]
    list(name=f[1], genes=f[-c(1,2)])
  })
  setNames(lapply(out,`[[`,"genes"), vapply(out,`[[`,"", "name"))
}
gene_names <- NULL
if (is.data.frame(ann) || is.data.table(ann)) {
  a <- as.data.table(ann)
  gid <- intersect(c("gene_id","Geneid","ENSEMBL"), names(a))[1]
  gnm <- intersect(c("gene_name","symbol","SYMBOL"), names(a))[1]
  if (!is.na(gid) && !is.na(gnm)) gene_names <- a[[gnm]][match(rownames(y),a[[gid]])]
}
if (!is.null(gene_names)) {
  hfile <- file.path(BASE,"reference/MSigDB/h.all.v2025.1.Hs.symbols.gmt")
  rfile <- file.path(BASE,"reference/MSigDB/c2.cp.reactome.v2025.1.Hs.symbols.gmt")
  for (obj in list(Hallmark=hfile,Reactome=rfile)) {
    coll <- names(obj); gfile <- obj[[1]]
    if (file.exists(gfile)) {
      gs <- read_gmt(gfile)
      idx <- lapply(gs, function(s) which(gene_names %in% s))
      idx <- idx[lengths(idx)>=5]
      for (j in seq_len(ncol(C))) {
        z <- camera(vqw, index=idx, design=design, contrast=C[,j],
                    inter.gene.cor=0.01, use.ranks=FALSE, allow.neg.cor=FALSE)
        z$Pathway <- rownames(z); z$FDR <- p.adjust(z$PValue,"BH")
        fwrite(z, file.path(OUT,paste0(colnames(C)[j],"_",coll,"_CAMERA.tsv")), sep="\t")
      }
    }
  }
}

# -------- Sex x Mtb interaction model --------
male <- as.integer(meta$sex=="Male")
t24  <- as.integer(meta$time=="24h")
t72  <- as.integer(meta$time=="72h")
mtb  <- as.integer(meta$condition=="Mtb")

baseD <- model.matrix(~ donor + time*condition, data=meta)
sexD <- cbind(baseD,
              sexM_time24=male*t24,
              sexM_time72=male*t72,
              sexM_Mtb=male*mtb,
              sexM_Mtb_time24=male*mtb*t24,
              sexM_Mtb_time72=male*mtb*t72)
stopifnot(qr(sexD)$rank == ncol(sexD))

SX <- matrix(0,ncol(sexD),6,dimnames=list(colnames(sexD),
 c("SexDifference_MtbResponse_2h","SexDifference_MtbResponse_24h","SexDifference_MtbResponse_72h",
   "SexDifference_ResponseChange_24h_vs_2h","SexDifference_ResponseChange_72h_vs_2h",
   "SexDifference_ResponseChange_72h_vs_24h")))
SX["sexM_Mtb",1] <- 1
SX["sexM_Mtb",2] <- 1; SX["sexM_Mtb_time24",2] <- 1
SX["sexM_Mtb",3] <- 1; SX["sexM_Mtb_time72",3] <- 1
SX["sexM_Mtb_time24",4] <- 1
SX["sexM_Mtb_time72",5] <- 1
SX["sexM_Mtb_time72",6] <- 1; SX["sexM_Mtb_time24",6] <- -1

vsx <- voomWithQualityWeights(y,sexD,plot=FALSE)
fsx <- eBayes(contrasts.fit(lmFit(vsx,sexD),SX),robust=TRUE)
sexsum <- list()
for (cn in colnames(SX)) {
  tt <- topTable(fsx,coef=cn,n=Inf,sort.by="P")
  tt$gene_id <- rownames(tt)
  fwrite(tt,file.path(OUT,paste0(cn,"_all_genes.tsv")),sep="\t")
  sexsum[[cn]] <- data.table(contrast=cn,FDR05=sum(tt$adj.P.Val<0.05))
}
fwrite(rbindlist(sexsum),file.path(OUT,"GSE189996_sex_interaction_summary.tsv"),sep="\t")

cat("GSE189996 final analysis complete.\n")
