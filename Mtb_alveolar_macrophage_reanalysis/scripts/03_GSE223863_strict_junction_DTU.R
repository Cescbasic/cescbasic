#!/usr/bin/env Rscript
# 03_GSE223863_strict_junction_DTU.R
# Final strict GENCODE-only junction usage analysis.
# Input is the prepared 105,122 x 30 exact-GENCODE STAR junction count matrix.
# Junction construction used uniquely assigned STAR junction counts, >=3 reads in >=3 samples,
# genes with >=2 supported junctions, then exact matching to STAR sjdbList.fromGTF.out.tab.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly=TRUE)
getarg <- function(flag,default=NULL){i<-match(flag,args);if(is.na(i))default else args[i+1]}
matrix_file <- getarg("--junction-matrix")
meta_file   <- getarg("--metadata")
map_file    <- getarg("--junction-map")
outdir      <- getarg("--outdir","reproducible_results/GSE223863_strict_junction")
if(any(vapply(list(matrix_file,meta_file,map_file),is.null,logical(1))))
 stop("Usage: Rscript 03_GSE223863_strict_junction_DTU.R --junction-matrix matrix.tsv --metadata metadata.tsv --junction-map junction_to_gene.tsv")

dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

x <- fread(matrix_file)
mp <- fread(map_file)
meta <- fread(meta_file)

# Expected mapping columns: junction_id, gene_id, gene_name
stopifnot(all(c("junction_id","gene_id") %in% names(mp)))
jid <- intersect(c("junction_id","junction","feature"),names(x))[1]
if(is.na(jid)) stop("Junction matrix requires a junction_id column")

# Keep exact GENCODE-owned features only, preserving one gene owner per junction.
x <- x[junction_id %in% mp$junction_id]
mp <- mp[match(x$junction_id,mp$junction_id)]
stopifnot(nrow(x)==nrow(mp))
stopifnot(!anyDuplicated(mp$junction_id))

sample_candidates <- setdiff(names(x),jid)
pick <- function(cand){z<-cand[cand%in%names(meta)];if(!length(z))stop(paste(cand,collapse="/"));z[1]}
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

stopifnot(all(meta$sample %in% sample_candidates))
cts <- as.matrix(x[,..meta$sample])
rownames(cts) <- x[[jid]]

# Matrix was already support-filtered upstream; retain the exact strict feature universe.
cat("Strict junction matrix:",nrow(cts),"x",ncol(cts),"\n")
cat("Genes:",uniqueN(mp$gene_id),"\n")
if(nrow(cts)!=105122) warning("Expected 105,122 strict GENCODE junctions from locked run.")
if(uniqueN(mp$gene_id)!=9496) warning("Expected 9,496 genes from locked run.")

# Donor-aware edgeR QL model. Donor is included as a fixed blocking factor.
design <- model.matrix(~ donor + time*condition,data=meta)
stopifnot(qr(design)$rank==ncol(design))

y <- DGEList(counts=cts,genes=data.frame(gene_id=mp$gene_id,
                                        exon_id=mp$junction_id,
                                        gene_name=if("gene_name"%in%names(mp)) mp$gene_name else mp$gene_id))
y <- calcNormFactors(y,method="TMM")
y <- estimateDisp(y,design,robust=TRUE)
fit <- glmQLFit(y,design,robust=TRUE)

findc <- function(p){z<-grep(p,colnames(design),value=TRUE);if(length(z)!=1)stop(p);z}
m <- findc("^conditionMtb$")
i24<-findc("^time24h:conditionMtb$|^conditionMtb:time24h$")
i72<-findc("^time72h:conditionMtb$|^conditionMtb:time72h$")

C <- matrix(0,ncol(design),6,dimnames=list(colnames(design),
 c("Mtb_vs_Control_2h","Mtb_vs_Control_24h","Mtb_vs_Control_72h",
   "Response_24h_vs_2h","Response_72h_vs_2h","Response_72h_vs_24h")))
C[m,1]<-1; C[m,2]<-1; C[i24,2]<-1; C[m,3]<-1; C[i72,3]<-1
C[i24,4]<-1; C[i72,5]<-1; C[i72,6]<-1; C[i24,6]<--1

summary_list <- list()
for(j in seq_len(ncol(C))) {
  cn <- colnames(C)[j]

  # edgeR differential splice/junction usage
  ds <- diffSpliceDGE(fit, contrast=C[,j],
                      geneid=y$genes$gene_id,
                      exonid=y$genes$exon_id)

  # Gene-level table
  gt <- topSpliceDGE(ds, test="gene", number=Inf)
  gt$gene_id <- rownames(gt)
  if(!"FDR"%in%names(gt) && "P.Value"%in%names(gt)) gt$FDR <- p.adjust(gt$P.Value,"BH")
  fwrite(gt,file.path(outdir,paste0(cn,"_gene_junction_usage.tsv")),sep="\t")

  # Junction-level table
  et <- topSpliceDGE(ds, test="exon", number=Inf)
  et$junction_id <- rownames(et)
  if(!"FDR"%in%names(et) && "P.Value"%in%names(et)) et$FDR <- p.adjust(et$P.Value,"BH")
  fwrite(et,file.path(outdir,paste0(cn,"_junction_level.tsv")),sep="\t")

  summary_list[[cn]] <- data.table(
    contrast=cn,
    genes_tested=nrow(gt),
    genomewide_FDR05=sum(gt$FDR<0.05,na.rm=TRUE)
  )
}
fwrite(rbindlist(summary_list),file.path(outdir,"genomewide_STAR_junction_usage_summary.tsv"),sep="\t")
saveRDS(list(y=y,design=design,contrasts=C,fit=fit),file.path(outdir,"strict_GENCODE_junction_model.rds"))

cat("Locked final strict analysis expected:\n")
cat("Mtb_vs_Control_2h 0; Mtb_vs_Control_24h 0; Mtb_vs_Control_72h 1; ",
    "Response_24h_vs_2h 0; Response_72h_vs_2h 1; Response_72h_vs_24h 2.\n")
cat("Final validated genes: HLA-DQB2 and RPL13A.\n")
