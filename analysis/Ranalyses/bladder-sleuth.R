library('sleuth')
library("biomaRt")

# Get transcript_id to gene name mappings
mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Read in the data 
data <- read.csv('/Users/arahuja/src/hammerlab/bladder-analyses/analyses/Ranalyses/sleuth_input.csv', stringsAsFactors = FALSE)

# Kallisto needs a sample columns
data$sample <- data$patient_id

# Fit sleuth model
# TODO: Add more covariates/confounders to this
so <- sleuth_prep(data, ~ is_benefit, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'is_benefitTrue')
sleuth_live(so)


