###############################################################################
## Establishing a pipeline for calculating signatures scores from bulk RNA seq data
###############################################################################
##
## attempting to follow along with: 
##    https://www.bioconductor.org/packages/devel/bioc/vignettes/signifinder/inst/doc/signifinder.html
##
###############################################################################

setwd("~/Dropbox (UFL)/GitHub/bulkRNAseq")

envirodir <- "~/Dropbox (UFL)/GitHub/bulkRNAseq/bin"
datadir <- "/Users/ferrallm/Dropbox (UFL)/collaborations/Rinaldi-Ramos/DC RNAseq Data - CRR/2023 Raw Data/03.Result_X202SC23100299-Z01-F001_Mus_musculus/Result_X202SC23100299-Z01-F001_Mus_musculus/3.Quant/1.Count"
resultdir <- "/Users/ferrallm/Dropbox (UFL)/collaborations/Rinaldi-Ramos/DC RNAseq Data - CRR/additional-analysis"

#### starting with count matrices from original analysis
## gene_count is the readcount list with the annotation of each gene
## gene_fpkm is the FPKM list with the annotation of each gene
## gene_fpkm_group is the list of FPKM values for each group; FPKM value: the adjusted read count (for biol rep)

### data input for package needs to be a SummarizedExperiment with the
###   normalized RNAseq counts

install.packages('BiocManager', lib=envirodir)
library('BiocManager', lib=envirodir)

BiocManager::install('signifinder', lib=envirodir)

install.packages('dplyr',lib=envirodir)
install.packages('tzdb',lib=envirodir)
install.packages('backports',lib=envirodir)

## instructions for constructing SummarizedExperiment: https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html#constructing-a-summarizedexperiment

# loading packages
library('tzdb', lib=envirodir)
library('backports', lib=envirodir)
library('signifinder', lib=envirodir)
library('SummarizedExperiment', lib=envirodir)
library('dplyr', lib=envirodir)

initDat <- read.csv(paste(datadir,"/gene_fpkm.csv",sep=""),header=TRUE)
nrows <- dim(initDat)[1]      ## number of genes
ncols <- 12                   ## 12 samples
counts <- initDat[,2:13]
meta <- initDat[,14:22]

exp <- SummarizedExperiment(assays=list(norm_expr=counts), metadata = meta)
rownames(exp) <- initDat[,1]

availSigns <- availableSignatures()
## ones that seem useful
## ExpandedImmune_Ayers, Chemokines_Messina, ImmunoScore_Hao, ImmuneCyt_Rooney,
## Tinflam_Ayers, ImmunoScore_Roh, ImmuneCyt_Davoli, IFN_Ayers

sigList <- c("ImmunoScore_Hao", "ImmunoScore_Roh", 
             "ImmuneCyt_Davoli", "ImmuneCyt_Rooney",
             "ExpandedImmune_Ayers", 
             "Chemokines_Messina", 
             "Tinflam_Ayers", 
             "IFN_Ayers")

exp <- multipleSign(dataset = exp, 
                     inputType = "rnaseq",
                     signature = sigList)





