# HUGO_map_GOA.R
# Map HUGO to Gene Ontology
# Team Chihuahua
#
#==== Reading GOA data and initial cleaning ====
# Reading in .gaf file
# skip 12 lines of comments
# GAF : Gene Association Format/File
# GPI : Gene Product Identifier + GPAD: Gene Product Association Data
GOA <- read.delim("./inst/extdata/goa_human.gaf",skip = 12, stringsAsFactors = FALSE,header = FALSE)

#extracting needed columns: Object Symbol, Qualifier, ID, and Aspect
GOA <- GOA[,c(3,4,5,9)]
colnames(GOA) <-  c("symbol","qualifier","GO_ID","Aspect")

#removing any rows with "NOT" qualifier - NOT GO_ID : making especially clear this GOID does
#not belong to the DB Oject symbol
#TODO: discuss removal (or not) of rows with "colocalizes_with" and "contributes_to" qualifier
GOA <-  GOA[!grepl("NOT",GOA$qualifier),]
#GOA <- GOA[!grepl("colocalizes_with",GOA$qualifier),]
#GOA <-  GOA[!grepl("contributes_to",GOA$qualifier),]

# ==== Mapping synonyms ====
#loading synMap
load("./inst/extdata/synMap.RData")

mySynonyms <- match(GOA$symbol,synMap$synonyms)
nrow(GOA)
nrow(synMap)

#GOA$symbol with a synonyms match <- HUGO symbols at GOA$symbol-synonyms index,
# where there was a match (not NA)
GOA$symbol[!is.na(mySynonyms)] <- synMap$symbols[mySynonyms[!is.na(mySynonyms)]]

# ==== Removing non-HUGO symbol rows ====
load("./inst/extdata/HUGOsymbols.RData")

GOA <- GOA[GOA$symbol %in% HUGOsymbols,]

# ==== List initialization and GOA data spliting ====
#spliting GOA data into 3 ontologies using "Aspect" column
#Aspect (column 9)
#refers to the namespace or ontology to which the GO ID (column 5) belongs;
#one of P (biological process, BP), F (molecular function, MF) or C (cellular component, CC)
#this field is mandatory; cardinality 1 (can only be one of the three)
GOA_MF <- GOA[GOA$Aspect=="F",]
GOA_BP <- GOA[GOA$Aspect=="P",]
GOA_CC <- GOA[GOA$Aspect=="C",]

#for list-based storage of symbols
MFsymList <- list()
BPsymList <- list()
CCsymList <- list()

#for list-based storage of pipe(|)-concatenated GO_IDs
GO_MFlist <- list()
GO_BPlist <- list()
GO_CClist <- list()

#Splits each dataframe into a list of dataframes, each with a common $symbol (argument f)
#extremely fast!
GOA_MFsplit<-split(GOA_MF,f = GOA_MF$symbol)
GOA_BPsplit<-split(GOA_BP,f = GOA_BP$symbol)
GOA_CCsplit<-split(GOA_CC,f = GOA_CC$symbol)

# ==== MF dataframe ====
#From each list of dataframes:
for(i in 1:length(GOA_MFsplit)){
  #store list of symbols (element 1 of $symbol)
  MFsymList[[i]]<-GOA_MFsplit[[i]]$symbol[1]
  
  #store list of pipe-concatenated GO_IDs
  #concatenate by paste0 (no spaces), collapse = "|"
  GO_MFlist[[i]]<-paste0(GOA_MFsplit[[i]]$GO_ID,collapse = "|")
}

#rbind() list of symbols into vector
# do.call extracts element from list, returns [matrix]
MFsym<-do.call(rbind,MFsymList)

#rbind() list of pipe-concatenated GO_IDs into vector
GO_MF<-do.call(rbind,GO_MFlist)

#generate dataframe with symbols and GOs, for matching with HUGO symbols
dfMF<-data.frame(symbol = MFsym, MF = GO_MF, stringsAsFactors = FALSE)

# ==== BP dataframe ====
#see MF dataframe section for comments
for(i in 1:length(GOA_BPsplit)){
  BPsymList[[i]]<-GOA_BPsplit[[i]]$symbol[1]
  GO_BPlist[[i]]<-paste0(GOA_BPsplit[[i]]$GO_ID,collapse = "|")
}
BPsym<-do.call(rbind,BPsymList)
GO_BP<-do.call(rbind,GO_BPlist)
dfBP<-data.frame(symbol = BPsym, BP = GO_BP, stringsAsFactors = FALSE)

# ==== CC dataframe ====
#see MF dataframe section for comments
for(i in 1:length(GOA_CCsplit)){
  CCsymList[[i]]<-GOA_CCsplit[[i]]$symbol[1]
  GO_CClist[[i]]<-paste0(GOA_CCsplit[[i]]$GO_ID,collapse = "|")
}
CCsym<-do.call(rbind,CCsymList)
GO_CC<-do.call(rbind,GO_CClist)
dfCC<-data.frame(symbol = CCsym, CC = GO_CC, stringsAsFactors = FALSE)

# ==== HUGO (output) dataframe ====
#initialize output dataframe with HUGO symbols
HUGO_GO<-data.frame(symbol = HUGOsymbols, stringsAsFactors = FALSE)

#match symbols with HUGO_GO and store GOs from each ontology dataframe
HUGO_GO$MF <- dfMF$MF[match(HUGO_GO$symbol,dfMF$symbol)]
HUGO_GO$BP <- dfBP$BP[match(HUGO_GO$symbol,dfBP$symbol)]
HUGO_GO$CC <- dfCC$CC[match(HUGO_GO$symbol,dfCC$symbol)]
