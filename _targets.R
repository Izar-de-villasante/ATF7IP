# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(renv)
library(tarchetypes) # Load other packages as needed. # nolint
# renv::activate()
ncores=max(RcppParallel::defaultNumThreads(),120)



# Set target options:
tar_option_set(
  packages = c("tibble","foreach","S4Vectors","data.table"), # packages that your targets need to run
  #imports = "cnv.methyl",
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess"# multiprocess or multicore, LSF, SGE, Slurm etc.
)

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)

# Load the R scripts with your custom functions:
#lapply(list.files("R", full.names = TRUE, recursive = TRUE), source)
# source("other_functions.R") # Source other scripts as needed. # nolint
source("R/functions.R")
.libPaths()

# produce_data <- function() {
#   expand.grid(samplesheet = c("a", "b"), model = c("c", "d"), normalization = c(1, 2, 3))
# }
# list(
#   tar_group_by(data, produce_data(), samplesheet, model),
#   tar_target(group, data, pattern = map(data))
# )


results_folder = "./results/batch2_rep"
analysis_folder = "./analysis/"

samplesheets <- c(
  sample_Sheet1="data/ss_H2009.rds",
  sample_Sheet2="data/ss_H358.rds",
  sample_Sheet3="data/ss_PDC11.rds",#EPIC
  #batch2:
  sample_sheet4="data/ss_PDC11_batch2.rds" # EPICv2
)



values <- tibble::tibble( # Use all possible combinations of input settings.
  method_function = rlang::syms(c("noob")),#, "pq","funn","noob_pq","Em")),
  data_paths = samplesheets,
  data_names = c("H2009","H358","PDC11","PDC11_batch2"),
  arraytype = c(rep("EPIC",3),"EPICv2")
  
  )

##########################
# Workaround Sample group with interaction:
for(i in samplesheets){
  dt<-readRDS(i)
  dt$interact<-paste0(dt$Type,"_",dt$Condition)
  saveRDS(dt,i)
}
vals<-rbind(values,values)
vals$groupvar<-c(rep("Sample_Group",nrow(values)),
              rep("interact",nrow(values))
              )
values<-vals[4,]
###########################
targets <- tar_map(
  values = values,
  names = data_names, #"data_source", # Select columns from `values` for target names.
  tar_target(params,make_results_dirs(subf=data_names, results_folder = results_folder, analysis_folder = analysis_folder)),
  tar_target(samplesheet_path, data_paths, format = "file"),
  tar_target(samplesheet, readRDS(samplesheet_path)),
  tar_target(ss,samplesheet[-1,]),
  tar_target(category,samplesheet[1,]),
  # tar_target(rgSet, cnv.methyl::read.metharray.exp.par(ss,folder="ana",extended=T,copy = F,arraytype = "EPIC",ncores = 8)),
  tar_target(nrgSet, {
    rgSet <- minfi::read.metharray.exp(base = NULL, targets = ss, extended = T, recursive = FALSE, verbose = FALSE, force = T)
    if(arraytype == "EPIC")rgSet@annotation<-c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b2.hg19")
    if(arraytype == "EPICv2")rgSet@annotation<-c(array="IlluminaHumanMethylationEPICv2",annotation="20a1.hg38")
    return(rgSet)}),
  tar_target(rgSet,name_rgset(nrgSet,ss),deployment = "main"),
  # Qc report:
  tar_target(QC_plots, qc(rgSet,sampGroups = "Sample_Name",sampNames="barcode", qc_folder = params[["qc_folder"]])),

  # Calculate purity:
  # tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
  tar_target(purity, {
    if(arraytype=="EPICv2"){
        
      library("randomForest")
      library(impute)
      RFpurify_ABSOLUTE<-cnv.methyl:::RFpurify_ABSOLUTE
      betas<-minfi::getBeta(rgSet)
      rownames(betas)<-unlist(sapply(strsplit(rownames(betas),"_"),"[")[1,])
      betas <- betas[match(rownames(RFpurify_ABSOLUTE$importance), 
                           rownames(betas)), , drop = FALSE]
      rownames(betas)<-rownames(RFpurify_ABSOLUTE$importance)
      betas <- tryCatch({
        betas <- impute::impute.knn(data = betas, k = 5)$data
      }, error = {
        betas[is.na(betas)] <- mean(betas, na.rm = T)
      })
      absolute <- stats::predict(RFpurify_ABSOLUTE, t(betas))
      return(absolute)
    }else if(arraytype=="EPIC" |arraytype=="450K") cnv.methyl::purify(myLoad=rgSet)
  }),
  tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="Sample_Group",qc_folder = params[["qc_folder"]])),

  tar_target(normalize, method_function(filtered)),
  tar_target(clean, prep(normalize)),
  tar_target(ss_clean, droplevels.data.frame( cbind(clean@colData,purity))),
  tar_target(save_ss_clean,write.table(ss_clean,paste0(params[["ss_clean_path"]],"/","ss_clean.csv"),quote = F,sep = ",") ),
  
  tar_target(plotvars, c(data.table::last(colnames(ss_clean)),"predictedSex",colnames(as.matrix(category[1,]))[as.matrix(category[1,])%in%c("batch","covs")])),
  
  tar_target(ann, {
    if(arraytype=="EPICv2"){
      annot<-data.table::as.data.table(minfi::getAnnotation(clean),keep.rownames = "ProbeID")
      annot<-annot[Rep_Num==1,]
      annot[,c("ProbeID","ext_Name"):=tstrsplit(ProbeID,"_")]
      return(annot)
    }else if(arraytype=="EPIC" |arraytype=="450K") minfi::getAnnotation(clean)
    }),
  tar_target(betas, {
    
      if(arraytype=="EPICv2"){
    library(data.table)
    bval <- data.table::as.data.table(minfi::getBeta(clean),keep.rownames = "ProbeID")
    bval[,c("ProbeID","ext_Name"):=tstrsplit(ProbeID,"_")]
    bval$ext_Name<-NULL
    bval<-bval[,lapply(.SD,mean),by=ProbeID]
    b<-as.matrix(bval[,-1],rownames=bval$ProbeID)
    return(b)
      }else if(arraytype == "EPIC") minfi::getBeta(clean)
      
    }),
  tar_target(top,top_beta(betas,n=1000)),
  tar_target(pca, pca_res(top)),
 
  tar_target(pca_corrplot,corpca(beta_top100 = top,
                                 metadata=ss_clean[,plotvars],
                                 title=paste0("PC1-6 correlations with ",data_names," clinical vars"))
             ),
  tar_target(save_pca_corrplot,save_plot(object=pca_corrplot,path=params[["corrplot_folder"]],filename=paste0(data_names,"_pca_corrplot.png"))),
  
  tar_target(bplots, bplot(pca,ss=ss_clean,col="Type",s="Condition",folder = params$bplots_folder)
             ),

  tar_target(model, mod(object = betas, group_var = groupvar,
                          metadata = ss_clean)
             ),

  tar_target(dmps_mod1, DMPextr(fit = model,
                                ContrastsDM = colnames(model$contrasts),
                                beta_normalized = betas,
                                p.value = 0.05,
                                mDiff = 0.15,
                                ann = ann,
                                writeOut = F,
                                ncores=3
                                            )),
  tar_target(save_dmps, writexl::write_xlsx(dmps_mod1,paste0(params$dmp_folder,data_names,"_dmps.xlsx"))),
  # tar_target(betas_DIF,betasdmps(betas,dmps_mod1,rgSet,anno_cols=c("Name","chr","pos","UCSC_RefGene_Name"))),
  # tar_target(betas_DIF_full,betasdmps(betas,ann,rgSet,anno_cols=c("Name","chr","pos","UCSC_RefGene_Name"))),
  # tar_target(write_betas_DIF, writexl::write_xlsx(betas_DIF,paste0(params$dmp_folder,"betas_DIF_",data_names,".xlsx"))),
  # tar_target(write_betas_DIF_full, writexl::write_xlsx(betas_DIF_full,paste0(params$dmp_folder,"betas_DIF_full_",data_names,".xlsx"))),

  
  tar_target(dmps_summary,{
    data.table::setDT(dmps_mod1)
    dmps_mod1[,list(Hyper=sum(Type=="Hyper"),Hypo=sum(Type=="Hypo")),by=c("Contrast")]
    }),

  tar_target(save_dmps_summary, writexl::write_xlsx(dmps_summary,paste0(params$dmp_folder,data_names,"_summary.xlsx"))),
  tar_target(dmpplot_mod1, plotDMP(dmps_mod1,path=params[["dmpplots_folder"]])),

  tar_target(dmrs, find_dmrs(object=dmps_mod1,model=model, bcutoff = 0.10, min.cpg=3,ncores=3)),
  tar_target(save_dmrs, 
             writexl::write_xlsx(as.data.frame(dmrs),
                                 paste0(params$dmrs_folder,"_",data_names,".xlsx"))),
  tar_target(sumaries,summarize(dmrs = dmrs, dmps = dmps_mod1,path = paste0(params$results_folder,"/",data_names,"/"))),
  tar_target(dmr_pathways, get_pathways(dmrs,res.folder =paste0(params[["pathway_folder"]],"/DMRs/"),savefile=TRUE)),
  tar_target(dmrs_summary,                                                       # Summary stats for DMRs
             summary_dmrs(
               dmrs,path=paste0(params$dmrs_folder,"full_dmrs_summary",data_names,".csv")),
             error = "continue"),
  
  tar_target(pathways, gopath(dmrs,all.cpg=rownames(betas),n=40,ann=ann)),
  
  # # 
  # tar_target(save_gopath, 
  #            writexl::write_xlsx(pathways[FDR<1,],
  #                                paste0(params$pathway_folder,data_names,"_Pathways.xlsx"))),
  NULL
)
# list(
#   tar_target(samplesheet_path, "data/ss_H358.rds", format = "file"),
#   tar_target(samplesheet, readRDS(samplesheet_path)),
#   # tar_target(ss,samplesheet),
#   tar_target(ss,samplesheet[-1,]),
#   tar_target(category,samplesheet[1,]),
#   tar_target(rgSet, cnv.methyl::read.metharray.exp.par(ss,folder="ana",extended=T)),
#   # Qc report:
#   tar_target(QC_plots, cnv.methyl::qc(rgSet,sampGroups = "condition")),
#
#   # Calculate purity:
#   tar_target(purity, cnv.methyl::purify(myLoad=rgSet)),
#   tar_target(filtered, filter(targets=ss, rgSet=rgSet,sampGroups="condition")),
#   targets,
#   # tar_target( alldmps,
#   #             list(apply(tidyr::expand_grid(c("dmps_ANA","dmps_SANDRA","dmps_WT"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
#   # )
#
#   #tar_combine(combined_gopath_ANA,
#   #list(apply(tidyr::expand_grid(c("gopath_ANA"),c("noob", "pq","funn","noob_pq","Em")),1,function(x)paste(x,collapse = "_")))
#   NULL
#   )
#
# subset(ss,!is.na(ss$ANA_dom1))
# technical_features <- c("Sample_Name", "organism", "Basename", "barcode")
#
#
#  mval<- getM(npq[,!is.na(ss$ANA_dom1)])
#  pheno<-ss[!is.na(ANA_dom1)]
#  mod <- model.matrix(
#    formula(paste(" ~ ANA_dom1 +", paste(batch,sep="+",collapse="+")))
#    ,data = pheno
#  )
#  mod0 <- model.matrix(
#    formula(paste(" ~", paste(batch,sep="+",collapse="+")))
#    ,data = pheno
#  )
#  sva.results <- sva(mval, mod, mod0)
# #
# design <- model.matrix(
#   formula(paste(" ~", paste(covs,sep="+",collapse="+")))
#   ,data = ss
#   )
# n.sv = sva::num.sv(betas,design,method="leek")
#
# svobj = sva(rgSet,mod,mod0,n.sv=n.sv)
#
# 
# # subset ss
# BiocManager::install(c("Biobase", "conumee", "GenomicRanges", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "impute", "IRanges", "limma", "maxprobes", "minfi", "SummarizedExperiment"))

combined <- tarchetypes::tar_combine(
  combined_summary,
  targets[["sumaries"]],
  command = dplyr::bind_rows(!!!.x,.id = "CL")
)
combined_ss <- tarchetypes::tar_combine(
  ss,
  targets[["ss_clean"]],
  command = dplyr::bind_rows(!!!.x,.id = "CL")
)

save_combined<-tar_target(
  save_combined_summary,{
    data.table::setorder(combined_summary,CL,Contrast)
    data.table::fwrite(combined_summary,paste0(results_folder,"/full_summary.csv"))
  }
)


list(targets,combined,save_combined,combined_ss)