library(data.table)
library(dplyr)
library(readxl)
library(devtools)

bp_stopifnot = getFromNamespace("stopifnot", "backports")

# Load data:
setDTthreads(0L)
idats_dir<-"data-raw/20221102_SG0002/"

Sample_sheet1<- fread("metadata/Sample_sheet_206467010100.csv")
Sample_sheet2<- fread("metadata/Sample_sheet_290622.csv")

Sample_sheet <- rbind(Sample_sheet1,Sample_sheet2)

# # Make a standard sample sheet for the samples. 
# ##It must contain the following columns:
# - Sample_Name
# - Basename
# 
# ##Other recommended columns:
# - Project
# - Pool_ID
# - Sample_Plate
# - Sample_Well
# - Sample_Group
# - Sentrix_ID
# - Sentrix_Position
# 
# ##Pheno columns:
# - Gender
# - Type
# - Condition
Basename<-paste0(idats_dir,Sample_sheet$Sentrix_ID,"/",Sample_sheet$Sentrix_ID,"_",Sample_sheet$Sentrix_Position)
bp_stopifnot("Missing idats " = file.exists(paste0(Basename,"_Grn.idat")) & file.exists(paste0(Basename,"_Red.idat")))

# stri_sub_all(Sample_sheet$Sample_Name, stri_locate_all_regex(Sample_sheet$Sample_Name, ' ', omit_no_match=TRUE)) <- '_'
ss<-data.table(Sample_Name=gsub(" ", "_",Sample_sheet$Sample_Name),
               Basename=Basename,
               barcode=basename(Basename),
               Sample_Plate=Sample_sheet$Sample_Plate,
               Pool_ID=Sample_sheet$Pool_ID,
               Sample_Well =Sample_sheet$Sample_Well,
               Sample_Group=gsub(" ", "_",Sample_sheet$Sample_Group),
               Sentrix_ID = Sample_sheet$SentrixID,
               Sentrix_Position = Sample_sheet$Sentrix_Position,
               Type = as.factor(sapply(strsplit(Sample_sheet$Sample_Group," "),"[",2)),
               Condition = as.factor(sapply(strsplit(Sample_sheet$Sample_Group," "),"[",1))
)



# ids:
ids <- c("Sample_Name",  # Unique identifier for study
         "barcode",      # Unique idats identifiersentrix_ID+Sentrix_position 
         "Basename")     # Path to idats/raw files.


batch <- c("CL") # Not applicable in this setup

covs <- c(
  "Type",           # Gene: Control= no gene present, case= has gene ATF7IP
  "Condition",      # Induction: Untreated: not inoculated Treated: 96h after inoculation  
  NULL)            

mgroups <-c(
  "Sample_Group"    # condition_type
  
)
# # Pheno:
ss$CL<-as.factor(sapply(strsplit(Sample_sheet$Sample_Name," "),"[",1))
# ss$Gene<-as.factor(sapply(strsplit(Sample_sheet$Sample_Name,"_"),"[",2))
# ss$Induction<-as.factor(sapply(strsplit(Sample_sheet$Sample_Name,"_"),"[",3))
# ss$Rep<-as.factor(sapply(strsplit(Sample_sheet$Sample_Name,"_"),"[",4))
# 

ss_clean <- ss[,.SD,.SDcols=c(ids,batch,covs,mgroups)]     
category <- as.list(c(rep("ids",length(ids)),rep("batch",length(batch)),rep("covs",length(covs)),rep("mgroups",length(mgroups))))
ss_clean <- rbindlist(list(category,ss_clean))
ss_clean[Type=="Control",Sample_Group:="Control"]
ss_clean[Type=="Case"&Condition=="Treated",Sample_Group:="Treat"]
ss_clean[Type=="Case"&Condition=="Untreated",Sample_Group:="Untreat"]
saveRDS(ss_clean[1:9,],paste0("./data/ss_","H2009",".rds"))
# In this case we don't want obs10 "H358_Ø_-Dox_A"
saveRDS(ss_clean[c(1,11:17),],paste0("./data/ss_","H358",".rds"))
# setkey(ss_clean, "CL")
# sapply(levels(ss_clean$CL),function(x)saveRDS(ss_clean[x,],paste0("./data/ss_",x,".rds")))

#usethis::use_data(ss, overwrite = TRUE)