# HUGO_map_prep.R
#
# Purpose: Sample processing of gene expression studies with RNA seq and
#            microarray platforms
# Version: 1.2
# Date:    2018 02 16
# Author:  Gregory Huang <gregory.huang@mail.utoronto.ca>
#
#
# ToDo:
# Notes:  R Code adapted from RPR-GEO2R  ABC learning unit.
#
#         Quantile normalization needs to cite:
#           Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)
#           A Comparison of Normalization Methods for High Density
#           Oligonucleotide Array Data Based on Bias and Variance.
#           Bioinformatics 19(2) ,pp 185-193.
#
#         Should we calculate DE instead and normalize that? There is otherwise
#         a danger of ML fitting the noise more than the otherwise sparse
#         signal. We could use (test - control) * pDE and set all insignificant
#         changes to 0.?
#
# ==============================================================================

# ====  PACKAGES  ==============================================================

if (! require(Biobase, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biobase")
  library(Biobase)
}

if (! require(GEOquery, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("GEOquery")
  library(GEOquery)
}

# for quantile normalization ...
if (! require(preprocessCore, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("preprocessCore")
  library(preprocessCore)
}

# ==============================================================================
#
# ==== Microarray data =========
#
# ==============================================================================

# Load in GSE gene expression data from GEO
GSE54017<- getGEO("GSE54017", GSEMatrix =TRUE, AnnotGPL=TRUE)
GSE54017 <- GSE54017[[1]]

# Save so there's no need to go to GEO every time
save(GSE54017, file="GSE54017.RData")
load("GSE54017.RData")

# ==== What data do we have?
nrow(exprs(GSE54017))      # 22215 rows
ncol(exprs(GSE54017))      # 8 (4 control replicates, 4 treatment replicates)
colnames(exprs(GSE54017))  # sample names, in the order of c1-t1-c2-t2-c3-t3-c4-t4

# Assess distributions
# define a color scheme: controls: greens, test: purples
c1 <- colorRampPalette(c("#C27E7E", "#816EBA"))(1)   # test
c2 <- colorRampPalette(c("#758AC9", "#82C9B6"))(1)   # ctrl

# a color vector for the 24 samples ...
myArrayCols <- c(rep(c2[1], 4),   # ctrl  4 reps
                 rep(c1[1], 4))   # test  4 reps

#reorder - 1,3,5,7 are control, and 2,4,6,8 are treatment
iReorder <- c(1,3,5,7,2,4,6,8)

boxplot(log(exprs(GSE54017)[ , iReorder]),
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE54017",
        outline = FALSE,
        col = myArrayCols)

# ==== extract columnd in new order

myEx <- exprs(GSE54017)[ , iReorder]
colnames(myEx) <- c("ca", "cb", "cc", "cd",
                    "ta", "tb", "tc", "td")

boxplot(log(myEx),
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE54017",
        outline = FALSE,
        col = myArrayCols)

# How to normalize? One way is to quantile normalize all replicate sets
# individually, but keep the trend differences between them.
myEx[ ,   1:4] <- normalize.quantiles(myEx[ ,   1:4], copy = FALSE)
myEx[ ,   5:8] <- normalize.quantiles(myEx[ ,   5:8], copy = FALSE)

boxplot(log(myEx),
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE54017-normQuant",
        outline = FALSE,
        col = myArrayCols)

# ==== Prepare annotation data
str(GSE54017@featureData@data)   # Have annotations been properly loaded ?
myAnnot <- GSE54017@featureData@data[ , c("ID", "Gene symbol")]
str(myAnnot)                                        # confirm
colnames(myAnnot) <- c("probeIDs", "symbols")       # rename
myAnnot$probeIDs <- as.character(myAnnot$probeIDs)  # convert to character
myAnnot$symbols  <- as.character(myAnnot$symbols)   # convert to character

any(is.na(myAnnot$symbols))       # FALSE
sum(myAnnot$symbols == "")        # 1069 probes are not annotated with a
# HUGO symbol. We will throw them out.

sum(grepl("/", myAnnot$symbols))  # 1223 probes are annotated with more than
# one symbol. We will throw them out too.

#just to keep track of number of rows
original_myEx <- nrow(myEx)
original_myAnnot <- nrow(myAnnot)

remove_dup <- grepl("/", myAnnot$symbols) #get index where there are dups (abc///abd)
myEx <- myEx[!remove_dup == TRUE,] # get rid of rows with dups in myEx
myAnnot <- myAnnot[!(remove_dup) == TRUE,] # same for myAnnot

remove_empty <- myAnnot$symbols == "" # get index where there are no symbols
myEx <- myEx[!(remove_empty)==TRUE,] # get rid of rows with no symbols in myEx
myAnnot <- myAnnot[!(remove_empty)==TRUE,] # same for myAnnot

#double check again, should be both zero
sum(myAnnot$symbols == "")   # 0
sum(grepl("/", myAnnot$symbols)) # 0


sum(duplicated(myAnnot$symbols[myAnnot$symbols != ""])) # 7421 duplicated
# symbols (not counting the un-annotated rows).

# How many unique symbols do these contain?
x <- unique(myAnnot$symbols[duplicated(myAnnot$symbols)]) # 4499 ...
# 7421 duplicated symbols, 4499 unique ones in the set of 7421

# ==================================================

# === More considerations of the symbol annotations. How do the existing
# symbols compare to our traget vector of gene symbols?

load("./inst/extdata/HUGOsymbols.RData")
load("./inst/extdata/synMap.RData")

# How many target symbols do we have covered?
sum(HUGOsymbols %in% unique(myAnnot$symbol)) # 11920; 11920/20347 = ~58%

# make a vector of missing symbols
missingSymbols <- HUGOsymbols[!(HUGOsymbols %in% unique(myAnnot$symbol))]

# What annotation symbols do we have that are NOT in our symbol list?
x <- unique(myAnnot$symbol) #12502
extraSymbols <- x[!(x %in% HUGOsymbols)] #582
head(extraSymbols, 50)

#let's check for symbols we have that don't match HUGO
matched_symbols <- match(myAnnot$symbols,HUGOsymbols)
#and check symbols we have that match with synonyms (precautionary)
matched_synonyms <- match(myAnnot$symbols, synMap$synoyms)

#get the locations of where the myAnnot$symbols match with synMap$synonyms
need_substitute <- !is.na(matched_synonyms)

#keep the synonym matches
myAnnot$symbols[need_substitute] <- synMap$symbols[matched_synonyms[need_substitute]]

#now we deal with the HUGO Symbol mismatches. There are 758 of them.
sum(is.na(match(myAnnot$symbols, HUGOsymbols)))

#remove the non-matches
ignore <- (is.na(match(myAnnot$symbols, HUGOsymbols)))
myEx <- myEx[!ignore,]
myAnnot <- myAnnot[!ignore,]
nrow(myEx)
nrow(myAnnot) #19165 lines left for both myEx and myAnnot

#check if there are any NAs left; 0
sum(is.na(match(myAnnot$symbols,HUGOsymbols)))

#Time to deal with the HUGO duplicates
#Use two duplicated() functions combined with an OR statement
#to make sure that every duplicated item is selected
from_beginning <- duplicated(myAnnot$symbols)
from_behind <- duplicated(myAnnot$symbols, fromLast = TRUE)
combined <- (from_beginning | from_behind)

all_duplicates <- myAnnot[combined, ] #all values with duplicates
all_uniques <- myAnnot[!combined, ]   #quick sanity check; all unique values
# add up to 19165

#now that we have a key for duplicates, remove them by the control group medians.
dfCombined <- myEx[combined,]
nrow(dfCombined) #11614 rows of symbols that have duplicates

#use this data frame to deal with the duplicate controls so medians are calculated
dfCombined_control <- data.frame(dfCombined[,1:4])

#medians for the duplicates are calculated
med <- apply(dfCombined_control, 1, median)

#bind the medians and symbol names for the controls to the df
dfCombined_control <- cbind(dfCombined_control, median = med)
dfCombined_control <- cbind(dfCombined_control, symbol = all_duplicates$symbols)

#Order the df by the medians from largest to smallest
dfCombined_control <- dfCombined_control[order(dfCombined_control$median, decreasing = TRUE),]

#Top-down duplicate search, this way, the duplicates with smaller medians get removed
dfCombined_control <- dfCombined_control[!duplicated(dfCombined_control$symbol),]

nrow(dfCombined_control) #4369 rows remain; that means we deleted 7245 rows of duplicates

# list of probeIDs that will be in our resulting data frame
final_rownames <- append(rownames(dfCombined_control), rownames(all_uniques))

#prepping for final data frame
result <- myEx[final_rownames,]
result_myAnnot <- myAnnot[final_rownames,]

#Calculate the control and treatment group means
ctrl_mean <- apply(result[,1:4], 1, mean)
treatment_mean <- apply(result[,5:8], 1, mean)

#===Final product===

#setup data frames
HUGO_result <-data.frame(symbol = HUGOsymbols, stringsAsFactors = FALSE)
setup <- data.frame(symbol = result_myAnnot$symbols,ctrl_mean,treatment_mean)

#add ctrl and treatment averages to the full list of HUGO symbols
HUGO_result$GSE54017.ctrl <- setup$ctrl_mean[match(HUGO_result$symbol,
                                                   setup$symbol)]
HUGO_result$GSE54017.treatment <- setup$treatment_mean[match(HUGO_result$symbol,
                                                             setup$symbol)]

coverage = (nrow(setup)/nrow(HUGO_result))*100
coverage #~58% coverage, same with earlier result; indeed, we kept the unique entries

# [END]


