library(biomaRt)
library(tidyr)
library(dplyr)
library(FDb.InfiniumMethylation.hg19)
library(AnnotationDbi)
library(hgu95av2.db)    ##for Human

### after download the test datasets and clinical information

# annotate genes in gene level
# read clinical data
au_pdac.clinical <- read.csv("PDAC_AU/donor1.tsv",sep = "\t")
au_pdac.specimen <- read.csv("PDAC_AU/specimen.tsv",sep = "\t")

# read mrna data
au_pdac.mrna <- read.csv("PDAC_AU/exp_seq.tsv",sep = "\t")


# get database
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# get annotation from database
genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),mart = ensembl)


##### australia pdac #####
#

### see if all the samples are from tumor tissue: no.14 is not
#try = unique(au_pdac.mrna[c("icgc_donor_id","icgc_specimen_id")])
# DO33168 has 2 specimen: SP110769,SP70820
au_pdac.mrna$icgc_donor_id = as.character(au_pdac.mrna$icgc_donor_id)
au_pdac.mrna[which(au_pdac.mrna$icgc_specimen_id=="SP110769"),"icgc_donor_id"] = "DO33168.1"
au_pdac.mrna$icgc_donor_id = as.factor(au_pdac.mrna$icgc_donor_id)

# primary tissue
primary.tissue <- unique(au_pdac.mrna$icgc_specimen_id) %in% au_pdac.specimen[au_pdac.specimen$specimen_type=="Primary tumour - solid tissue",]$icgc_specimen_id
# which specimen is not tumor tissue
non.primary <- unique(au_pdac.mrna$icgc_specimen_id)[!primary.tissue]

# which specimen is not PDAC
PDAC.tissue <- unique(au_pdac.mrna$icgc_specimen_id) %in% au_pdac.specimen[au_pdac.specimen$tumour_histological_type=="Pancreatic Ductal Adenocarcinoma",]$icgc_specimen_id
non.PDAC <- unique(au_pdac.mrna$icgc_specimen_id)[!PDAC.tissue]

# drop non-primary and non-pdac and without clinical info
au_pdac.mrna1 = au_pdac.mrna[-which(au_pdac.mrna$icgc_specimen_id %in% non.primary),]
au_pdac.mrna1 = au_pdac.mrna1[-which(au_pdac.mrna1$icgc_specimen_id %in% non.PDAC),]
au_pdac.mrna1 = au_pdac.mrna1[-which(au_pdac.mrna1$icgc_specimen_id == "SP110769"),]

# add gene annotation
au_pdac.mrna2 = left_join(au_pdac.mrna1,genes,by = c("gene_id"="ensembl_gene_id"))
# donor id and specimen id are one-by-one matched now

# keep useful info
au_pdac.mrna2 = au_pdac.mrna2[c("gene_id","raw_read_count","hgnc_symbol","icgc_donor_id")]
# 56280 gene_ids
au_pdac.mrna2 = spread(au_pdac.mrna2,icgc_donor_id,raw_read_count)

# remove gene ids without gene symbols: 56280 -> 34044 features
au_pdac.mrna2 = au_pdac.mrna2[au_pdac.mrna2$hgnc_symbol!="",]

# aggregate
au_pdac.mrna2 = aggregate(au_pdac.mrna2[c(3:69)],by = list(au_pdac.mrna2$hgnc_symbol), FUN =function(x) x=max(x) )
au_pdac.mrna3 <- au_pdac.mrna2[rowSums(au_pdac.mrna2 == 0) <= 14,]

# au pdac clinical data
au_pdac.clinical2 = au_pdac.clinical[c("icgc_donor_id","donor_sex","donor_vital_status","donor_survival_time")]
# match clinical information
id_mrna = intersect(au_pdac.clinical2$icgc_donor_id,colnames(au_pdac.mrna3))
au_pdac.mrna3 = au_pdac.mrna3[c("Group.1",id_mrna)]
#au_pdac.mrna3 = au_pdac.mrna3[-1,]
# match clinical info
au_pdac.clinical2 = au_pdac.clinical2[which(au_pdac.clinical2$icgc_donor_id %in% id_mrna),]

# save data
write.csv(au_pdac.mrna3,"pancreatic_cancer/external dataset/PDAC_AU/rnaseq.tsv")

######## au mrna array #########
# read mrna data
au_mrna.array <- read.csv("PDAC_AU/exp_array.tsv",sep = "\t")
au_array.annot = read.csv("PDAC_AU/array_gene_annot2.txt",sep = '\t')

# primary tissue
primary.tissue <- unique(au_mrna.array$icgc_specimen_id) %in% au_pdac.specimen[au_pdac.specimen$specimen_type=="Primary tumour - solid tissue",]$icgc_specimen_id
# which specimen is not tumor tissue
non.primary <- unique(au_mrna.array$icgc_specimen_id)[!primary.tissue]
# which specimen is not PDAC
PDAC.tissue <- unique(au_mrna.array$icgc_specimen_id) %in% 
  au_pdac.specimen[au_pdac.specimen$tumour_histological_type=="Pancreatic Ductal Adenocarcinoma",]$icgc_specimen_id
non.PDAC <- unique(au_mrna.array$icgc_specimen_id)[!PDAC.tissue]

# drop non-pdac & non-primary 
au_mrna.array1 = au_mrna.array[-which(au_mrna.array$icgc_specimen_id %in% non.primary),]
au_mrna.array1 = au_mrna.array1[-which(au_mrna.array1$icgc_specimen_id %in% non.PDAC),]

# merge with annotation
au_mrna.array1 <- left_join(au_mrna.array1,au_array.annot[c("ID","ILMN_Gene")],c("gene_id"="ID"))

#select useful info
au_mrna.array2 = au_mrna.array1[c("icgc_donor_id","normalized_expression_value","ILMN_Gene","gene_id")]
au_mrna.array3 = spread(au_mrna.array2,icgc_donor_id,normalized_expression_value)

# aggregate
au_mrna.array3 = aggregate(au_mrna.array3[c(3:66)],by = list(au_mrna.array3$ILMN_Gene), FUN =function(x) x=max(x) )
# save data
write.csv(au_mrna.array3,"pancreatic_cancer/external dataset/PDAC_AU/mrna_array.tsv")
au_mrna.array3 = read.csv("pancreatic_cancer/external dataset/PDAC_AU/mrna_array.tsv")

#match clinical info
id_mrna = intersect(au_pdac.clinical$icgc_donor_id,colnames(au_mrna.array3))
au_pdac.clinical1 = au_pdac.clinical[which(au_pdac.clinical$icgc_donor_id %in% colnames(au_mrna.array3)),]

##### australia pdac methylation ######

au_pdac.methy <- read.csv("PDAC_AU/meth_array.tsv",sep = "\t")
# get database
hm450 <- get450k()
# get probe id list
probenames <- unique(au_pdac.methy$probe_id)
# get nearest annotation
probes <- hm450[probenames]
symbol <- getNearestTSS(probes)
# 1500 bp within TSS
symbol = symbol[symbol$distance<1500,]
symbol$cpg_id = rownames(symbol)
au_pdac.methy = left_join(au_pdac.methy,symbol[c("cpg_id","nearestGeneSymbol")],by = c("probe_id"="cpg_id"))

# primary tissue
primary.tissue <- unique(au_pdac.methy$icgc_specimen_id) %in% au_pdac.specimen[au_pdac.specimen$specimen_type=="Primary tumour - solid tissue",]$icgc_specimen_id
# which specimen is not tumor tissue
non.primary <- unique(au_pdac.methy$icgc_specimen_id)[!primary.tissue]
# which specimen is not PDAC
PDAC.tissue <- unique(au_pdac.methy$icgc_specimen_id) %in% au_pdac.specimen[au_pdac.specimen$tumour_histological_type=="Pancreatic Ductal Adenocarcinoma",]$icgc_specimen_id
non.PDAC <- unique(au_pdac.methy$icgc_specimen_id)[!PDAC.tissue]

au_pdac.methy1 = au_pdac.methy[-which(au_pdac.methy$icgc_specimen_id %in% non.primary),]
au_pdac.methy1 = au_pdac.methy1[-which(au_pdac.methy1$icgc_specimen_id %in% non.PDAC),]

# select info
au_pdac.methy2 = au_pdac.methy1[c("icgc_specimen_id","icgc_donor_id","probe_id","methylation_value","nearestGeneSymbol")]
au_pdac.methy2 = au_pdac.methy2[!is.na(au_pdac.methy2$nearestGeneSymbol),]

# with treatment before and after surgery, filter
# DO33168:SP110770, DO32860:SP110736, DO32936:SP110751, DO32875:SP72262
no_treat = c("SP110770","SP110736","SP110751","SP72262")
au_pdac.methy3= au_pdac.methy2[!(au_pdac.methy2$icgc_specimen_id %in% no_treat),]
# duplicate measurments
au_pdac.methy3= au_pdac.methy3[!(au_pdac.methy3$icgc_donor_id %in% c("DO34240","DO33128")),]
au_pdac.methy3$icgc_specimen_id <- NULL
au_pdac.methy3 = spread(au_pdac.methy3,icgc_donor_id,methylation_value)

# filter out missing value
#  remove genes with more than 20% missing value
au_pdac.methy4 <- au_pdac.methy3[rowSums(is.na(au_pdac.methy3)) <= 13, ]
au_pdac.methy4$probe_id <-NULL
#aggregate
# impute
imputed <- impute.knn(as.matrix(au_pdac.methy4[,c(2:65)]))
au_pdac.methy4[,c(2:65)] <- as.data.frame(imputed$data)

au_pdac.methy4 <- aggregate(au_pdac.methy4[,c(2:65)],by=list(au_pdac.methy4$nearestGeneSymbol), FUN=function(x) x=max(x))

# march clinical id
common_id = intersect(au_pdac.clinical$icgc_donor_id,colnames(au_pdac.methy4))
au_pdac.methy4 = au_pdac.methy4[c("Group.1",common_id)]
au_pdac.clinical4 = au_pdac.clinical[au_pdac.clinical$icgc_donor_id %in% colnames(au_pdac.methy4),]

write.csv(au_pdac.methy4,"pancreatic_cancer/external dataset/PDAC_AU/methylation.tsv")

########## GEO mrna ##########

geo.mrna = read.csv("pancreatic_cancer/external dataset/PDAC_GEO_mrna/GSE62452_series_matrix.txt",sep = '\t')
geo.mrna_clinic = read.csv("PDAC_GEO_mrna/GSE62452_clinical.txt",sep = '\t')
geo.annot = read.csv("PDAC_GEO_mrna/GPL6244_annot.txt",sep = '\t')

# match id 
geo.mrna = left_join(geo.mrna,geo.annot[c("ID","Gene.symbol")],c("ID_REF"="ID"))
# filter gene without gene symbol
geo.mrna = geo.mrna[geo.mrna$Gene.symbol!=""]

geo.mrna1 = read.csv("PDAC_GEO_mrna/GSE62452_series_matrix1.txt",sep = '\t')
try2 = geo.mrna1[132:168]
try2 = as.data.frame(t(try2))
names(try2) = geo.mrna1$ID_REF
try3 = gather(try2,key = "ID_REF",value = "genesymbol")
try3 = try3[try3$genesymbol!="",]
geo.mrna1$ID_REF <- as.character(geo.mrna1$ID_REF)
geo.mrna1 = left_join(try3,geo.mrna)
geo.mrna1$Gene.symbol = NULL
geo.mrna2 = aggregate(geo.mrna1[c(3:132)],by = list(geo.mrna1$genesymbol), FUN =function(x) x=max(x) )
 
# clinical data
# clean data, filter sample without clinical info 
geo.mrna_clinic$time = sub("survival months: ","",geo.mrna_clinic$time)
geo.mrna_clinic$vital_status = sub("survival status: ","",geo.mrna_clinic$vital_status)
geo.mrna_clinic = geo.mrna_clinic[geo.mrna_clinic$time != "" & geo.mrna_clinic$vital_status != "",]
geo.mrna_clinic = geo.mrna_clinic[geo.mrna_clinic$time != "?" & geo.mrna_clinic$vital_status != "?",]
geo.mrna_clinic = geo.mrna_clinic[geo.mrna_clinic$Sample!="E138-Tn",]
# 
# geo.mrna3 = geo.mrna2[as.character(geo.mrna_clinic$geo_accession)]
geo.mrna3$Group.1 = geo.mrna2$Group.1
write.csv(geo.mrna3,"pancreatic_cancer/external dataset/PDAC_GEO_mrna/GSE62452_final.tsv")

######### GEO mirna ##########

# load clinical and mirna data
geo.clinic = read.csv("pancreatic_cancer/external dataset/PDAC_GEO_mirna/GSE62498_clinical.txt",sep = '\t')
geo.mirna = read.csv("pancreatic_cancer/external dataset/PDAC_GEO_mirna/GSE62498_processed_data.txt",sep = '\t',header = FALSE)

#clean clinical info, filter sample without survival data
geo.clinic$time = sub("survival months: ","",geo.clinic$time)
geo.clinic$vital_status = sub("survival status: ","",geo.clinic$vital_status)
geo.clinic = geo.clinic[geo.clinic$time != "NA" & geo.clinic$vital_status != "NA",]
geo.clinic$time = as.numeric(geo.clinic$time) / 12 * 365

# match id with omics,remove 0 before number 005->5, 042-> 42
# geo.clinic$Sample_id = sub("E0","E",geo.clinic$Sample_id) #*2

## fill in missing and remove 0 
geo.mirna = as.data.frame(t(geo.mirna))
rownames(geo.mirna) = geo.mirna$ID_REF
geo.mirna$ID_REF = NULL

# filter data with 20% missing or more
geo.mirna2 = geo.mirna[colSums(geo.mirna == "?") <= 14]
geo.mirna2[geo.mirna2=="?"] = NA
geo.mirna3 = sapply(geo.mirna2,function(x) as.numeric(levels(x)[x]) ) 

# impute missing data
imputed <- impute.knn(geo.mirna3)
geo.mirna3 <- imputed$data
geo.mirna3 <- as.data.frame(geo.mirna3)
names(geo.mirna3) <- tolower(names(geo.mirna3))

#match id of clinical and omics
# remove rebundant after "_"
geo.mirna3[,1] = sub("_[0-9]*","",geo.mirna3[,1])
# remove space
geo.mirna3[,1] = sub(" ","",geo.mirna3[,1])
geo.mirna3 = geo.mirna3[geo.mirna3[,1] %in% geo.clinic$Sample_id,]
# save data
write.csv(geo.mirna3,"pancreatic_cancer/external dataset/PDAC_GEO_mirna/mirna.tsv")
