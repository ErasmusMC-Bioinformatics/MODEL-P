
############ download data ############


# clinical data
filename_biosClin <- DownloadBiospecimenClinicalData(cancerType = "PAAD", 
                                                     saveFolderName = "./BiospecimenClinicalData", 
                                                     outputFileName = "clinical")
clinical = read.table("BiospecimenClinicalData/clinical__nationwidechildrens.org_clinical_patient_paad.txt", header = TRUE,sep = "\t")

# download PAADdata
pan_dna_methy <- DownloadMethylationData(cancerType = "PAAD",assayPlatform = "methylation_450",saveFolderName = "./pan_dna_methy")
pan_mrna <- DownloadRNASeqData(cancerType = "PAAD",assayPlatform = "gene.normalized_RNAseq",saveFolderName = "./pan_mrna")
pan_mirna <- DownloadmiRNASeqData(cancerType = "PAAD",assayPlatform = "mir_HiSeq.hg19.mirbase20",saveFolderName = "./pan_mirna")
pan_mutation <- DownloadSomaticMutationData(cancerType = "PAAD", assayPlatform = "somaticMutation_DNAseq", saveFolderName = "./pan_mutation")


library(dplyr)
# read and process data
Methylation450Data <- ProcessMethylation450Data(inputFilePath = "TCGA-Assembler/pan_dna_methy/PAAD_methylation_450.txt",
                                                outputFileName = "READ__humanmethylation450", outputFileFolder =
                                                  "ProcessedData.methy")
methy <- read.table("ProcessedData/READ__humanmethylation450.txt",sep = '\t',header = TRUE)

# miRNA
miRNASeqRawData <- ProcessmiRNASeqData(inputFilePath =  "TCGA-Assembler/pan_mirna/PAAD_mirbase20.txt",
                                       outputFileName = "READ__illuminahiseq_mirnaseq", 
                                       outputFileFolder = "ProcessedData.mirna")
mirna <- read.table("ProcessedData/READ__illuminahiseq_mirnaseq__RPM.txt",sep = '\t',header = TRUE)
rownames(mirna) <- mirna$V1
mirna$V1 <- NULL

# mutation
mutationData <- ProcessSomaticMutationData(inputFilePath = "pan_mutation/PAAD__somaticMutation_maf.txt",
                                           outputFileName = "MutationData", 
                                           outputFileFolder ="ProcessedData.mutation")
mutation.orig <- read.table("D:/1_Jie_OIO/pancreatic cancer/ProcessedData.mutation/MutationData_mutationLevel.txt",sep = "\t",header = TRUE)#




# # cut samples names
# colnames(mutation)[3:180] <- substr(colnames(mutation)[3:180],1,12)
# mutation <- mutation[mutation$GeneSymbol %in% c("KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA1", "BRCA2","TP63"),]
# # aggregate mutation to gene level
# 
# mutation$geneMutation <- paste(mutation$GeneSymbol,"-",mutation$Variant_Classification)
# try = aggregate(mutation[c(3:180)], by=list(mutation$geneMutation), sum)
# rownames(try) <- try$Group.1
# try$Group.1 <- NULL
# 
# try <- as.data.frame(t(try))
# try1 = merge(clinical_new3, try, by=0, all=TRUE)
# try1 = try1[-which(is.na(try1$bcr_patient_barcode)),]
# try1$Row.names = NULL
# try1 = try1[c(24:37)]
# try1 = na.omit(try1)
# 
# fisher.p = lapply(try1[c(2:14)], function(x) fisher.test(try1$class,x)$p.value)



# RNA
mrna <- ProcessRNASeqData(inputFilePath = "TCGA-Assembler/pan_mrna/PAAD_RNAseq.txt", 
                          outputFileName = "READ__illuminahiseq_rnaseqv2__GeneExp", 
                          outputFileFolder = "ProcessedData.mrna", 
                          dataType = "geneExp", verType = "RNASeqV2")

mrna <- read.table("ProcessedData/READ__illuminahiseq_rnaseqv2__GeneExp.txt",sep = '\t',header = TRUE)


######## dataset aggregation #########


# modify smaple id to match, drop useless
methy_id = colnames(methy)[5:199]
methy_id = substr(methy_id,1,16)

mrna_id = colnames(mrna)[3:185]
mrna_id = substr(mrna_id,1,16)

mirna_id = colnames(mirna)[2:183]
mirna_id = substr(mirna_id,1,16)

mutation <- mutation.orig[c(1,9,20:197)]
mutation_id = colnames(mutation)[3:180]
mutation_id = substr(mutation_id,1,16)
colnames(mutation)[3:180] = mutation_id

# select common id of methylation, mirna and rna
common_id = Reduce(intersect,list(methy_id,mrna_id,mirna_id,mutation_id))
# find patient samples and exclude normal
types <- gsub("TCGA.[A-Z0-9]*.[A-Z0-9]*.", "",common_id)
unique(types)
# "01A" "11A" "06A"
# 01A and 06A are patients, but only 01A is primary tumor
tumor.indices <- which(types %in% "01A")
common_id <- common_id[tumor.indices]
# remove non-primary tumor
methy_new = methy[common_id]
mirna_new = mirna[common_id]
mrna_new = mrna[common_id]

# get patient id
colnames(methy)[5:199] = substr(methy_id,1,12)
colnames(mrna)[3:185] = substr(mrna_id,1,12)
colnames(mirna)[2:183] = substr(mirna_id,1,12)
common_id = substr(common_id,1,12)

# change name of cinical id to match genomic information
clinical_id = clinical$bcr_patient_barcode[3:187]
clinical_id = gsub("-",".",clinical_id)
common_id = intersect(common_id,clinical_id)
clinical$bcr_patient_barcode = gsub("-",".",clinical$bcr_patient_barcode)

# select patient having clinical, mrna, mirna and methylation information
common_id = intersect(common_id,clinical_id)
common_id = c("GeneSymbol",common_id)
methy_new = methy_new[common_id]
mirna_new = mirna_new[common_id]
mrna_new = mrna_new[common_id]

# mutation
mutation$TCGA.HZ.A9TJ.06A <- NULL
colnames(mutation)[3:179] = substr(colnames(mutation)[3:179],1,12)
mutation_new <- mutation[c("GeneSymbol","Variant_Classification",common_id)]
mutation_new <- mutation_new[c("GeneSymbol","Variant_Classification",intersect(colnames(mutation_new),clinical_new3$bcr_patient_barcode))]
# related genes
mutation_new <- mutation_new[mutation_new$GeneSymbol %in% c("KRAS", "TP53", "CDKN2A", "SMAD4", "BRCA1", "BRCA2","CD44","TP63","PRSS1","CHI3L2"),]#
# subtype of mutation data
mutation_subtype = clinical_new3[clinical_new3$bcr_patient_barcode %in% colnames(mutation_new)[3:142],c("bcr_patient_barcode","class")]
#clinical_mutation = clinical_new3[clinical_new3$bcr_patient_barcode %in% colnames(mutation_new)[3:142],c("class","survival","vital_status")]

mutation_new1 = data.frame(t(mutation_new[3:142]))
mutation_new1$bcr_patient_barcode = rownames(mutation_new1)
mutation_subtype = left_join(mutation_subtype,mutation_new1)
#geo_surv.mut = survfit(Surv(survival, vital_status) ~ class, clinical_mutation)
#
#ggsurvplot(geo_surv.mut,pval = TRUE)


odd.ratio = rnorm(104)
pvalue = rnorm(104)

for (i in 3:106){
  if (all(mutation_subtype[[i]]==0)){
  #if (i == 18 | i == 46 | i == 64 |i == 68| i == 76  | i == 104){
    print(i)
    odd.ratio[i] = 10
    pvalue[i] = 10
    next
  }
  #x = unlist(mutation_new[i,3:142])
  #logit.mut = glm(mutation_subtype~x,family=binomial(link='logit'))
  fishertest  = fisher.test(mutation_subtype[[i]],mutation_subtype$class)
  odd.ratio[i-2] = fishertest$estimate
  #odd.ratio[i] = exp(coef(logit.mut))[2]
  pvalue[i-2] = fishertest$p.value 
  #pvalue[i] = summary(logit.mut)$coefficients["x", "Pr(>|z|)"]
}
odd.ratio
pvalue

mutation_subtype[c("survival","vital_status")] = clinical_new3[clinical_new3$bcr_patient_barcode %in% mutation_subtype$bcr_patient_barcode,c("survival","vital_status")]
summary(coxph(Surv(survival, vital_status) ~ class + X17, mutation_subtype))

clinical_new = clinical[which(clinical$bcr_patient_barcode %in% common_id),]
clinical_new = clinical_new[c("bcr_patient_barcode","vital_status","last_contact_days_to","death_days_to","cause_of_death","histologic_diagnosis","histologic_diagnosis_other")]
# modify survival time and status column to prepare for cox-ph
clinical_new$last_contact_days_to <- as.numeric(as.character(clinical_new$last_contact_days_to))
clinical_new$death_days_to <- as.numeric(as.character(clinical_new$death_days_to))

clinical_new["last_contact_days_to"][is.na(clinical_new["last_contact_days_to"])] <- 0
clinical_new["death_days_to"][is.na(clinical_new["death_days_to"])] <- 0
clinical_new$survival <- clinical_new$last_contact_days_to + clinical_new$death_days_to


######### features preprocessing #########


######## mrna #######
# filter out data with more than 20% 0 or missing value

# mrna 20530 features -> 17225 features
mrna_filter <- mrna_new[rowSums(mrna_new == 0) <= 40,]

# delete rows without genesymbols 17225 features -> 17211 features
mrna_filter <- mrna_filter[mrna_filter$GeneSymbol != "?",]

# aggragate mrna features, 17211 features -> 17188 features
mrna_aggre <- aggregate(mrna_filter[,c(2:178)],by=list(mrna_filter$GeneSymbol), FUN=function(x) x=max(x))


write.csv(mrna_aggre,"pancreatic cancer/input data/mrna.csv")
######### mirna ##########
# filter out data with more than 20% 0 or missing value
# mirna 1870 features -> 429 features
mirna_filter <- mirna_new[rowSums(mirna_new == 0) <= 40, ]


###### methylation ########
# filter out data with more than 20% 0 or missing value
# methylation 526732 variables -> 432516 variables
methy_filter <- methy_new[rowSums(is.na(methy_new)) <= 40, ]

# delete rows without genesymbols 432516 features -> 346440 features
methy_filter <- methy_filter[methy_filter$GeneSymbol != "",]
# fill in missing value in DNA methtlation
# missings <- sum(is.na(methy_filter))

# 111783 missing values
library(impute)
imputed <- impute.knn(as.matrix(methy_filter[,c(2:178)]))
methy_filled <- imputed$data
methy_filled <- as.data.frame(methy_filled)
methy_filled <- cbind("GeneSymbol"=factor(methy_filter$GeneSymbol),methy_filled)

# aggregate methylation features, 346440 features -> 20980 features 
methy_aggre <- aggregate(methy_filled[,c(2:178)],by=list(methy_filled$GeneSymbol), FUN=function(x) x=max(x))

methy_aggre2 <- aggregate(methy_filled[,c(2:178)],by=list(methy_filled$GeneSymbol), FUN=function(x) x=mean(x))


###### scaling for svm #######
write.table(mrna_aggre,"pancreatic cancer/input data/mrna.txt",row.names = FALSE, sep = "\t",col.names = FALSE)

library(SCAN.UPC)
mrna_scan = UPC_RNASeq("pancreatic cancer/input data/mrna.txt")
rownames(mrna_scan) = colnames(mrna_aggre)[2:178]
write.csv(mrna_scan,"pancreatic cancer/input data/mrna_upc.csv")

# mrna_scan_array = sapply(mrna_aggre[2:178], UPC_Generic)
# mrna_scan_array = as.data.frame(mrna_scan_array)
# mrna_scan_array$GeneSymbol = mrna_aggre$Group.1
# write.csv(mrna_scan_array,"pancreatic cancer/input data/mrna_array_upc.csv")

mirna_scan = sapply(mirna_filter[2:178], UPC_Generic)
mirna_scan = as.data.frame(mirna_scan)
mirna_scan$GeneSymbol = mirna_filter$GeneSymbol
write.csv(mirna_scan,"pancreatic cancer/input data/mirna_upc.csv")


# output centered data
write.csv(methy.train1,"svm data/methylation_train.csv")
write.csv(mrna.train1,"svm data/mrna_train.csv")
write.csv(mirna.train1,"svm data/mirna_train.csv")
write.csv(methy.test1,"svm data/methylation_test.csv")
write.csv(mrna.test1,"svm data/mrna_test.csv")
write.csv(mirna.test1,"svm data/mirna_test.csv")



