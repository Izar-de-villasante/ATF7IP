library(data.table)
library(dplyr)
library(readxl)
library(devtools)

bp_stopifnot = getFromNamespace("stopifnot", "backports")

# Load data:
setDTthreads(0L)
idats_dir<-"data-raw/20221102_SG0002/batch2/JuanMOrillas_28_08_2023/"

JuanMorillas_posicionchips_28_08_2023 <- read_excel("data-raw/20221102_SG0002/batch2/JuanMOrillas_28_08_2023/JuanMorillas_posicionchips_28_08_2023.xlsx", 
                                                    col_names = FALSE)
Sample_sheet <- data.table::as.data.table(JuanMorillas_posicionchips_28_08_2023)
colnames(Sample_sheet)<- c("Sample_Name","Sentrix_ID","Sentrix_Position")
Sample_sheet[,c("CL","Gene","Induction","Rep"):=lapply(tstrsplit(Sample_Name," "),as.factor)]
Sample_sheet[,Condition := as.factor(ifelse(Induction=="96h","Treated","Untreated"))]
Sample_sheet[,Type := as.factor(ifelse(Gene=="ATF7IP","Case","Control"))]
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
Basename<-paste0(idats_dir,Sample_sheet$Sentrix_ID,"_",Sample_sheet$Sentrix_Position)
bp_stopifnot("Missing idats " = file.exists(paste0(Basename,"_Grn.idat")) & file.exists(paste0(Basename,"_Red.idat")))

# stri_sub_all(Sample_sheet$Sample_Name, stri_locate_all_regex(Sample_sheet$Sample_Name, ' ', omit_no_match=TRUE)) <- '_'

ss_curated <- Sample_sheet[,.(Sample_Name=gsub(" ", "_",Sample_sheet$Sample_Name),
               Basename=Basename,
               barcode=basename(Basename),
               Sample_Plate=Sample_sheet$Sample_Plate,
               Pool_ID=Sample_sheet$Pool_ID,
               Sample_Well =Sample_sheet$Sample_Well,
               Sentrix_ID = Sample_sheet$Sentrix_ID,
               Sentrix_Position = Sample_sheet$Sentrix_Position

)]
# ss <- merge(ss_curated,Sample_sheet,by=c("Sentrix_ID","Sentrix_Position"),all.x=T)
common_names <- intersect(names(ss_curated),names(Sample_sheet))
# different names and index names 
diff_names <- c(setdiff(names(Sample_sheet),nm1),c("Sentrix_ID","Sentrix_Position"))
ss <- ss_curated[ss[, ..diff_names], on = .(Sentrix_ID, Sentrix_Position),  nomatch = 0]
ss[,Sample_Group:= ifelse(Type == "Control","Control",ifelse(Condition=="Treated","Treated",as.character(Type))) ]

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

ss_clean <- ss[,.SD,.SDcols=c(ids,batch,covs,mgroups)]     
category <- as.list(c(rep("ids",length(ids)),rep("batch",length(batch)),rep("covs",length(covs)),rep("mgroups",length(mgroups))))
ss_clean <- rbindlist(list(category,ss_clean))
ss_clean$arraytype <- "EPICv2"
saveRDS(ss_clean,paste0("./data/ss_","PDC11_batch2.rds"))

#usethis::use_data(ss, overwrite = TRUE)