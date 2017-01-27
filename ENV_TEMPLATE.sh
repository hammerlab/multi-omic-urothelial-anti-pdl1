# Parent directory for the bladder data set
# Common ones:
# - Demeter:
#export DATA_DIR=/demeter/scratch/datasets/mskcc/bladder/
# - GCloud:
#export DATA_DIR=/mnt/checkpoint-trials/bladder
# - Your own:
export DATA_DIR=

# By default, `analyses/util/data.py` creates/uses the folder `data`
# within the bladder-analyses folder, so there is no absolute need
# to set the `CACHE_DATA_DIR` if you don't need to customize this.
# 
# If you would like to re-use already existent cache files/folders,
# either uncomment the following line (this is not recommended!)
# export CACHE_DATA_DIR=$DATA_DIR/cohorts-cache
# or copy the folder above over to a custom location
# and set that as your own cache dir:
#export CACHE_DATA_DIR=

# Change these only if you are using a custom directory structure
export VCF_DATA_DIR=$DATA_DIR/vcfs-indelrealigned-bqsr
export BAM_DNA_DATA_DIR=$DATA_DIR/bams
export BAM_RNA_DATA_DIR=$DATA_DIR/rna-aligned/star-aligned-bams
export KALLISTO_DATA_DIR=$DATA_DIR/rna-aligned/kallisto
export CUFFLINKS_DATA_DIR=$DATA_DIR/rna-aligned/star-cufflinks
export TCR_DATA_DIR=$DATA_DIR/tcr
export POLYPHEN_DATA_DUMP=$DATA_DIR/annotation/polyphen2/polyphen-2.2.2-whess-2011_12.sqlite
