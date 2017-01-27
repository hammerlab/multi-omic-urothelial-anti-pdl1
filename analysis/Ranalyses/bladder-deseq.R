library(tximport)
library(readr)
library(DESeq2)
library(gage)
library(gageData)

# Load ensembl information to map transcript ids to gene ids
mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
t2g <- dplyr::rename(t2g, TXNAME = ensembl_transcript_id,
                     GENEID = ensembl_gene_id)

# Load kallisto data
data <- read.csv('~/src/hammerlab/bladder-analyses/analyses/Ranalyses/sleuth_input.csv', stringsAsFactors = FALSE, row.names = 1)
files <- file.path(data$path, "abundance.tsv")
names(files) <- rownames(data)

# Aggregate to gene level using tximport
txi <- tximport(files, type = "kallisto", tx2gene = t2g, reader = read_tsv)
write.csv(as.data.frame(txi$counts), 'bladder-kallisto-tximport.csv', row.names=TRUE)
rownames(data) <- colnames(txi$counts)


###
# Perform gene level differential analysis using DESeq2
###

# Convert the transcript level data to gene-level for use with DESeq
dds <- DESeqDataSetFromTximport(txi, data, ~is_benefit)
dds = DESeq(dds)

# Compare expression based on 'is_benefit' and order by p-value
deseq_results = results(dds, contrast=c("is_benefit", "True", "False"))
deseq_results = deseq_results[order(deseq_results$pvalue),]

# Convert from Ensembl IDs to gene symbols
deseq_results$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(deseq_results), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
deseq_results$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(deseq_results), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
deseq_results$name = mapIds(org.Hs.eg.db,
                    keys=row.names(deseq_results), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
summary(deseq_results)

library(pathview)
library(GSVAdata)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

foldchanges = deseq_results$log2FoldChange
names(foldchanges) = deseq_results$entrez

keggres <- gage(foldchanges, gsets = kegg.sets.hs)

# Look at differentially expressed pathways
lapply(keggres, head)

library(Biobase)
library(genefilter)
library(limma)
library(GSVA)

###
# Perform GSVA analysis using the Hallmark gene sets
###

# Get Ensembl to Entrez gene ID map
egIDs <- id2eg(row.names(txi$counts), category = "ENSEMBL", org = "Hs", pkg.name = NULL)

# Rename the expression rows to their Entrez ID
rownames(txi$counts) <- egIDs[, 2]

# Load the Hallmark gene sets
hallmark <- GSEABase::getGmt("~/h.all.v5.1.entrez.gmt")
r <- gsva(txi$counts, hallmark, min.sz=10, max.sz=500, verbose=TRUE, rnaseq=TRUE)$es.obs
design <- model.matrix(~ data$is_benefit)
fit <- lmFit(r, design)
fit <- eBayes(fit)

adjPvalueCutoff <- 0.001
# Find the top differentially expressed pathways with multiple hypothesis correction
allGeneSets <- topTable(fit, coef="data$is_benefitTrue", number=Inf)
DEgeneSets <- topTable(fit, coef="data$is_benefitTrue", number=Inf, p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)

