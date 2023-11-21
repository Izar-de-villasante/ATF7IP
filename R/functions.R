# global opts:
options(future.globals.maxSize= 1891289600)

# amke results folders:
make_results_dirs <- function(results_folder = "./results/", analysis_folder = "./analysis/",subf){
  params<-list(
    results_folder = results_folder,
    qc_folder  = paste(analysis_folder,subf,"intermediate/QC/",sep= .Platform$file.sep ),
    ss_clean_path = paste(analysis_folder,subf,sep= .Platform$file.sep),
    bplots_folder = paste(results_folder,subf,"plots/pca/bplots/",sep= .Platform$file.sep),
    corrplot_folder = paste(results_folder,subf,"plots/pca/corrplot/",sep= .Platform$file.sep),
    dmp_folder = paste(results_folder,subf,"dmps/",sep= .Platform$file.sep),
    dmpplots_folder = paste(results_folder,subf,"dmps/",sep= .Platform$file.sep),
    dmrs_folder = paste(results_folder,subf,"dmrs/",sep= .Platform$file.sep),
    pathway_folder = paste(results_folder,subf,"gopath/",sep= .Platform$file.sep)
    
  )
  sapply(params,function(x)  dir.create(x,recursive=T,showWarnings = F))
  
  return(params)
}


name_rgset<-function(res,targets,newname=NULL,exclude=NULL,idcol="barcode"){
  require("Biobase")
  require(SummarizedExperiment)
  require(data.table)
  targets<-data.table::as.data.table(targets)
  
  cn<-targets[[idcol]]
  colnames(res)<-cn
  colnames(res@assays@data$Green)<-cn
  colnames(res@assays@data$Red)<-cn
  data.table::setkey(targets,"barcode")
  
  pheno <- methods::as(targets, "DataFrame")
  rownames(pheno)<-pheno$barcode
  stopifnot(rownames(pheno)==colnames(res))
  # res@colData[[idcol]] <- sapply(res@colData$barcode, function(x) substr(x,nchar(x)-18,nchar(x)))
  # remove bad samples
  res@colData<-pheno
  
  if(!is.null(newname))colnames(res)<-res@colData[[newname]]
  res[,!colnames(res) %in%exclude]
  return(res)
}



#Normalization functions:

noob <- function(rgSet){
  minfi::preprocessNoob(rgSet)
}

funn <- function(rgSet){
  minfi::preprocessFunnorm(rgSet)
}
noob_pq <- function(rgSet){
  minfi::preprocessNoob(rgSet)%>% minfi::preprocessQuantile()
}
pq <- function(rgSet){
  minfi::preprocessQuantile(rgSet)
}

# Em2 <- function(rgSet, arraytype = NULL){
#   pd <- rgSet@colData
#   qc <- ENmix::QCinfo(rgSet,distplot = F)
#   mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
#   mSetSqn <- ENmix::mpreprocess(rgSet = rgSet,impute = T)
#   if(is.null(arraytype)){
#     arraytype <- ifelse(500000 < nrow(mSetSqn), "EPIC", "450K" )
#   }   
#   
#   if(arraytype=="EPIC"){
#     mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(
#       mat = mSetSqn,
#       array = "IlluminaHumanMethylationEPIC", 
#       annotation = "ilm10b4.hg19"
#     )
#   }else mSetSqn <- minfi::makeGenomicRatioSetFromMatrix(mSetSqn)
#   
#   mSetSqn@colData <- pd
#   return(mSetSqn)
# }

# Em <- function(rgSet, arraytype = NULL){
#   pd <- rgSet@colData
#   qc <- ENmix::QCinfo(rgSet,distplot = F)
#   mdat <- ENmix::preprocessENmix(rgSet, QCinfo = qc)
#   mSetSqn <- minfi::mapToGenome(mdat)
#   mSetSqn@colData <- pd
#   return(mSetSqn)
# }

filter<-function(targets, rgSet,sampGroups=NULL,sampNames="Sample_Name",frac=0.1,pval=0.01,remove_sex=TRUE,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
  # requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  # requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
  n <- ncol(rgSet)
  #o <- rev(order(sampNames))
  #rgSet <- rgSet[, o]
  #sampNames <- sampNames[o]
  if (is.null(sampGroups)) g <- rep(1, n) else g <- targets[[sampGroups]]
  if (is.null(g)) g<-1:n
  g <- factor(g)
  
  pal =  c(
    "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
             "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
             "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
             "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
             
  
  # Check quality of the combined probe signal (testing against negative controls)
  # 1. Remove low quality samples. Defaults more than 10% we throw samples away.
  detP <- minfi::detectionP(rgSet, type = "m+u")
  grDevices::jpeg(file = paste0(qc_folder,"mean_detection_pvalues.jpeg"))#,width = 480,height = 480)
  ylabels<-colnames(detP)
  par(mar=c(max(4.1,max(nchar(ylabels))/2.2) ,4.1 , 4.1, 2.1))
  
  barplot(colMeans(detP), col=pal[g], las=2,
          cex.names=0.8, ylim=c(0,max(0.002,max(colMeans(detP))*2)), main ="Mean detection p-values")
  graphics::abline(h=0.05,col="red")
  graphics::legend("topleft", legend=levels(g), fill=pal[1:length(levels(g))],
                   bg="white")
  grDevices::dev.off()
  
  # 2. Removing low-quality samples (with fraction of probes not passing pval)
  bad_samples <- colnames(detP)[colSums(detP >=pval)/nrow(detP) > frac]
  if(length(bad_samples)>0){
    warning("The following samples will be discarded since they fail to pass the p-value filter ( ",
            frac*100,"% of the probes with p-val >", pval, "): \n ", paste(bad_samples,collapse = ", " ))
    rgSet <- rgSet[,setdiff(colnames(detP),bad_samples)]
  }else{
    cat("All samples passed detection P-value filter")
  }
  # 3. Removing low-quality probes (with p-value below pval)
  bad_probes<-which(rowSums(detP < pval) < ncol(rgSet)*(1-frac))
  rgSet <- rgSet[-c(bad_probes),]
  if(length(bad_samples)>0){
    warning("The following probes will be discarded since more than", frac*100,
            "% of the samples have detection p-values > ", pval, "): \n ", paste(bad_samples,collapse = ", " ))
  }else{
    cat("All samples passed detection P-value filter")
  }
  return(rgSet)
}

#' Generate qc plots
#' @title generate qc plots for signal distribution prior to filtering
#' @param rgSet rgSet object containing channel intenisty values
#' @param sampNames variable containing barcodes or ids matching colnames of the rgsetdata
#' @param sampGroups variables to use for coloring groups
#' @param qc_folder path to the folder where plots will be saved
#' @return plots 
#' @author izar de Villasante
#' @export
#'
qc <- function(rgSet,sampGroups=NULL, sampNames= "Sample_Name",qc_folder="analysis/intermediate/QC/"){
  if (!(sampNames %in% names(rgSet@colData))) sampNames <- colnames(rgSet)
  # n <- ncol(rgSet)
  # o <- rev(order(sampNames))
  # rgSet <- rgSet[, o]
  # sampNames <- sampNames[o]
  # if (is.null(sampGroups))
  #   sampGroups <- rep(1, n)
  # sampGroups <- sampGroups[o]
  dir.create(qc_folder,recursive=T,showWarnings = F)
  minfi::qcReport(rgSet = rgSet,
                  pdf = paste0(qc_folder,"Report.pdf"),
                  sampGroups = rgSet@colData[[sampGroups]],
                  sampNames = rgSet@colData[[sampNames]])
  if(length(unique(rgSet@colData[[sampGroups]]))>1){
    grDevices::png(file = paste0(qc_folder,"density_plot.png"),
                   width = 480, # The width of the plot in inches
                   height = 620) # The height of the plot in inches
    pal =  c(
      "#191919", "#F0A0FF", "#0075DC", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
               "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
               "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
               "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
               minfi::densityPlot(
                 rgSet, sampGroups = rgSet@colData[[sampGroups]],main = "Beta",
                 pal =pal
               )
               grDevices::dev.off()
               # grDevices::png(file = paste0(qc_folder,"bean_plot.png"),
               #                width = 480, # The width of the plot in inches
               #                height = 620) # The height of the plot in inches
               # minfi::densityBeanPlot(
               #   rgSet, sampGroups = rgSet@colData[[sampGroups]],
               #   pal=pal
               # )
               # 
               # grDevices::dev.off()
               
  }
  mSet <- minfi::preprocessRaw(rgSet)
  qc   <- minfi::getQC(mSet)
  grDevices::png(file = paste0(qc_folder,"mean_qc.png"),   # The directory you want to save the file in
                 width = 480, # The width of the plot in inches
                 height = 480) # The height of the plot in inches
  minfi::plotQC(qc)
  
  grDevices::dev.off()
}

prep<-function(mSetSqn,remove_sex=TRUE,pval=0.01,arraytype=NULL,qc_folder= "analysis/intermediate/QC"){
  
  
  # 4. Removing probes with known SNPs at CpG site
  mSetSqn <-  minfi::mapToGenome(mSetSqn)
  mSetSqn <- minfi::dropLociWithSnps(mSetSqn)
  
  # 5. Removing cross reactive probes
  if(!is.null(mSetSqn@colData$arraytype)){
    if(unique(mSetSqn@colData$arraytype)%in%"EPICv2"){
      xreactive <- readRDS("data/mask.rds")
      mSetSqn <-  mSetSqn[setdiff(rownames(mSetSqn),xreactive),]
    }else{mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)}
  }else{
    mSetSqn <-  maxprobes::dropXreactiveLoci(mSetSqn)
  }
  # 6. Sex. prediction & removal
  mSetSqn$predictedSex <- minfi::getSex(mSetSqn, cutoff = -2)$predictedSex
  if(remove_sex){
    if(!is.null(arraytype)){anno<-minfi::getAnnotation(mSetSqn)
    }else{
      anno<-minfi::getAnnotation(mSetSqn)
      anno<-anno[!(anno$chr %in% c("chrX","chrY")),]
    }
  }
  return(mSetSqn)
}


# Add info:

# Purity:
#cnv.methyl::purify()

# Celltype:
#1. cellCounts <- FlowSorted.Blood.450k::estimateCellCounts(rgSet)
#2. FlowSorted.Blood

# Copy Number Variation:
#cnv.methyl::Kc_get(ss )

# Sex:
# minfi::get_sex()

# Age:
# 
# Enmix<-function(rgSet2){
#   qcE<-ENmix::QCinfo(rgSet2)
#   mdat<-ENmix::preprocessENmix(rgSet2, bgParaEst="oob", dyeCorr="RELIC",
#                                QCinfo=qc, nCores=6)
# }

# # obtaining the beta values
# beta_values <- getBeta(gmSet)
# colnames(beta_values) <- metadata$sample
# 
# saveRDS(beta_values, file = "results/beta_values.rds")
# 
# # PRINCIPAL COMPONENT ANALYSIS
# 
# # selecting the top 100 most variable CpG sites
top_beta <- function(beta_values, n=1000){
  sdv <- apply(beta_values, 1, sd)
  top100 <- names(head(sort(sdv,decreasing=T), n))
  beta_top100 <- beta_values[top100,]
  return(beta_top100)
}
pca_res <- function(beta_top100,scale=T, center=T){
  prcomp(t(beta_top100), scale=scale, center=center)
}



corpca <- function(beta_top100,metadata,vars=NULL,title='PC1-6 clinical correlations'){
  requireNamespace("PCAtools")
  p<-PCAtools::pca(beta_top100,metadata = metadata, removeVar = 0.1)
  if(is.null(vars)){
    vars<-names(p$metadata)[sapply(p$metadata,function(x){
      !any(is.na(x))& length(unique(x))>1
    })
    ]
    vars[sapply(p$metadata,function(x)length(unique(x)))>1]
  }
  PCAtools::eigencorplot(p,
                         components = PCAtools::getComponents(p, 1:6),
                         metavars = vars,
                         col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
                         cexCorval = 0.7,
                         colCorval = 'white',
                         fontCorval = 2,
                         posLab = 'bottomleft',
                         rotLabX = 45,
                         posColKey = 'top',
                         cexLabColKey = 1.5,
                         scale = TRUE,
                         main = title,
                         colFrame = 'white',
                         plotRsquared = FALSE)
  
}


bplot<-function(pca,ss,colgroup,s,combs=NULL, tit= NULL,folder = "analysis/pca/bplots/"){
  library(ggplot2)
  library(gplots)
  library(ggrepel)
  library(ggfortify)
  ss<-droplevels.data.frame(ss)
  pal=c("#1b9e77","#d95f02", "#191919", "#0075DC", "#F0A0FF", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
                 "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
                 "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
                 "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
                 # pal =  c(
  #   "#191919", "#0075DC", "#F0A0FF", "#993F00", "#005C31", "#5EF1F2", "#FF0010",
  #   "#FFCC99", "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", "#003380",
  #   "#FFA405", "#FFA8BB", "#426600", "#2BCE48", "#4C005C", "#00998F", "#E0FF66",
  #   "#740AFF", "#990000", "#FFFF80", "#FFE100", "#FF5005")
  lapply(colgroup, function(f) {
    f_len <- length(unique(with(ss,get(f))))
    cols <- pal[1:f_len]
    names(cols)<- unique(with(ss,get(f)))
    if(is.null(combs))combs<-combn(4,2)
    ss<-as.data.frame(ss)
    #rownames(ss)<-ss$Sample_Name
    for(i in 1:dim(combs)[2]){
      if(is.null(tit))tit <- paste0("colored.by.",f, "_shape.", s)
      ap<-ggplot2::autoplot(pca, x=combs[1,i], y=combs[2,i], data = ss, colour=f,shape=s,alpha=0.7,size=1)+
        geom_text_repel(aes(label = Sample_Name, color = with(ss,get(f))),
                        show.legend = FALSE, size = 1.5,max.overlaps = Inf,segment.size=0.2,min.segment.length = 0.8,point.size = 0.5)+
        #scale_color_brewer(palette = "Paired")+
        scale_color_manual(values = cols)+
        #geom_point(aes(size=0.2))+
        labs(colour=f)+ 
        
        ggtitle(tit)+
        theme_bw(base_size = 7)+
        theme(legend.key=element_blank(), legend.key.size=unit(1,"point"))
      
      
      dir.create(folder)
      #plot(ap)
      ggsave(paste0(folder,tit,i,".png"),plot=ap,width = 5.56, height = 2.80,units="in" )
    }
  })
  return(folder)
}

# Surrogate analysis:
surrogate<-function(grset,pheno,condition){
  mval<- getM(grset)
  pheno<-pData(grset)
  mod <- model.matrix(~as.factor(condition), data=pheno)
  mod0 <- model.matrix(~1, data=pheno)
  sva.results <- sva::sva(mval, mod, mod0)
}


# # PCA on the top 100 sites
# pca_res <- prcomp(t(beta_top100), scale=T, center=T)
# pca_all <- prcomp(t(beta_values), scale=T, center=T)
# 
# 
# 
# 
# 
# ## plotting PC1 and PC2 by condition
# autoplot(pca_res, x=1, y=2, data=metadata, colour="vascular_type", shape="type")+
#   geom_text_repel(aes(label=sample, color=vascular_type),hjust=-0.2, vjust=0, show.legend=F, size=3.5)+
#   labs(colour="Tissue", shape="Type")+
#   xlim(c(-0.5,0.3))+
#   theme_bw()+
#   ggtitle("PCA by tissue")
# 

#' Generate models
#' @title construct models and contrasts with limma
#' @param object your object containing beta values
#' @param group_var the variable used as independent variable
#' @param covs the set of variables to use as confounders
#' @param metadata the metadata or sample sheet
#' @param set a boolean vector to subset the observations
#' @param gr the group
#' @return fit2 ebayes model 
#' @author izar de Villasante
#' @export
#'
mod <- function(object, group_var, covs=NULL, metadata,set = TRUE,gr=NULL,pairwise = T,
                singular=F){
  data.table::setDT(as.data.frame(metadata))
  cont_sing=cont_pair=gr_cont_sing=gr_cont_pair=NULL
  metadata<-subset(metadata,set)
  metadata<-droplevels(metadata)
  object <- object[,metadata$barcode]
  covs_formula<-NULL
  if (!is.null(covs)&length(covs)>0)covs_formula<-paste0("+",paste0(covs,collapse=" + ",sep=""))
  design <- model.matrix( 
    formula(
      paste("~ 0 +" , paste0(group_var),covs_formula,sep= " " )
    ),
    data = metadata
  )
  fit <- limma::lmFit(object,design)
  cols <- with(metadata,paste0(group_var, unique(get(group_var))))
  if(pairwise == T){
    cont_pair <- apply(combn(cols,2),2,function(x) paste(x,collapse = "-"))
    
    if(!is.null(gr)) { 
      gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
      gr_cont_pair <- apply(combn(gr_cols,2),2,function(x) paste(x,collapse = "-"))
    }
  }
  
  
  if(singular == T){
    cont_sing<-apply(combn(cols,length(cols)-1),2,function(x){
      var <- setdiff(cols,x)
      group <- contgroup(group_var,levels=x)
      contrast <- paste0(var,"-", group)
      return(contrast)
    } )
    if(!is.null(gr)) { 
      gr_cols <- sapply(gr,function(x)contgroup(x,colnames(design)))
      gr_cont_sing <- apply(combn(gr_cols,length(gr_cols)-1),2,function(x){
        var <- setdiff(gr_cols,x)
        group <- contgroup(group_var,levels=x)
        contrast <- paste0(var,"-", group)
        return(contrast)
      } )
    }
  }
  
  cont <- c(cont_sing,cont_pair)
  gr_cont <- c(gr_cont_sing,gr_cont_pair)
  contMatrix <- limma::makeContrasts(
    contrasts=c(cont,gr_cont),
    levels=colnames(design)
  )
  # rename contrasts:
  # GR: 
  if(!is.null(gr)) colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), gr_cols, gr, vectorize_all = FALSE)
  
  # Singular 1vsMean.
  large <- colnames(contMatrix) %in% c(cont_sing,gr_cont_sing)
  colnames(contMatrix)[large] <- sapply(
    colnames(contMatrix)[large], function(x) paste0("sing_",strsplit(x,"-")[[1]][1]))
  
  # remove group_var prefix:
  colnames(contMatrix) <- stringi::stri_replace_all_fixed(colnames(contMatrix), group_var, "", vectorize_all = FALSE)
  fit2 <- limma::contrasts.fit(fit, contMatrix)
  fit2 <- limma::eBayes(fit2)
  return(fit2)
}


contgroup<-function(name,levels){
  cols<-levels[grepl(name, levels, fixed = TRUE)]
  l<-length(cols)
  paste0("(",paste(cols,collapse = "+"),")/",l)
}


betasdmps <- function(betas,dmps,rgSet,anno_cols=c("Name","chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group")){
  require(data.table)
  dmps<-data.table::setDT(as.data.frame(dmps))
  ss<-rgSet@colData
  pvals <- minfi::detectionP(rgSet)
  p<-pvals[dmps$Name,]
  colnames(p)<-paste0(ss[colnames(pvals),"Sample_Name"],".Detection_Pval")
  b<-betas[dmps$Name,]
  colnames(b)<-paste0(ss[colnames(betas),"Sample_Name"],".AVG_betas")
  d<-dmps[,.SD,.SDcols=anno_cols]
  df<-cbind(d,b,p)
  df<-data.table::setDT(as.data.frame(df))
  data.table::setorder(df,chr,pos)
  # df[order(chr,pos)]
  return(df)
}

plotDMP <- function(DMPann,names,path=NULL){
  if(NROW(DMPann) >0){
    library(ggplot2)
    library(data.table)
    data.table::setDT(DMPann)
    DMPresults <- data.frame(table(DMPann[ ,c("Contrast","Type")]))
    # plot DMPs (hypo/hyper)
    g1<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
      geom_bar(position="dodge", stat= "identity")+
      theme_bw()+
      scale_fill_manual(values=c("red", "skyblue")) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))+
      labs(x = "", y = "count", fill='Methylation')+
      ggtitle('Differently methylated probes')
    
    # ggplot2::ggsave(g1,paste0("analysis/DMPplots/",names,".png"))
    
    # plot with facets
    g2<-ggplot2::ggplot(DMPresults, aes(Contrast, Freq, fill = Type)) +
      geom_bar(position="dodge", stat= "identity")+
      facet_wrap(.~Type, scales = "free_x") +
      theme_bw()+
      scale_fill_manual(values=c("red", "skyblue")) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))+
      labs(x = "", y = "count", fill='Methylation')+
      ggtitle('Differently methylated probes')
    # ggplot2::ggsave("analysis/DNA_methylation/DMP_p0.01_m0.3.facet.png")
    
    # plot proportion of DMPs in CGI
    DMP_annCGI <- data.frame(DMPann[ ,c("Contrast","Type", "Relation_to_Island")])
    g3<-ggplot2::ggplot(DMP_annCGI, aes(Contrast, fill = Relation_to_Island)) +
      facet_wrap(.~Type, scales = "free_x") +
      geom_bar(position ="fill", width = 0.8) +
      theme_bw() +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      ylab("DMPs") +
      xlab("")
    # ggplot2::ggsave("analysis/DNA_methylation/DMP_annCGI.png")
    
    #table(DMPann[ ,c("Relation_to_Island","Type", "Contrast")])
    
    # plot proportion of DMPs in genomic elements
    DMPann$UCSC_RefGene_Group[which(DMPann$UCSC_RefGene_Group == "")] <- "."
    DMP_annGenomic<-DMPann[,.(
      UCSC_RefGene_Group_short = unlist(lapply(strsplit(UCSC_RefGene_Group, ";"),'['))),
      by = c("Contrast","Type")]
    g4<-ggplot2::ggplot(DMP_annGenomic, aes(Contrast, fill = UCSC_RefGene_Group_short)) +
      facet_wrap(.~Type, scales = "free_x") +
      geom_bar(position = "fill", width = 0.8) +
      theme_bw() +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(fill = "RefGene") +
      ylab("DMPs") +
      xlab("")
    # ggplot2::ggsave("analysis/DNA_methylation/DMP_annGenomic.png")
    
    # table(DMPann[ ,c("UCSC_RefGene_Group_short","Type", "Contrast")])
    plt_list<-list(g1,g2,g3,g4)
    n <- c( "DMP_count.png","DMP_count_facet.png","DMP_annCGI.png", "DMP_annGenomic.png")
    
    if(!is.null(path)){
      sapply(1:length(plt_list),function(x){
        ggplot2::ggsave(
          filename = n[x],
          plot = plt_list[[x]],
          device = NULL,
          path = path,
          scale = 1,
          width = 4,
          height = 5,
          units = c("in"),
          dpi = 600,
          limitsize = TRUE,
          bg = NULL
        )
      })
    } 
    return(plt_list)
  }else{
    warning("empty data.frame please, try again modifying filter params.")
    
  }
  
}
save_plot <-function(object,filename,path){
  grDevices::png(file = paste(path,filename,sep="/"),
                 width = 480, # The width of the plot in inches
                 height = 620) # The height of the plot in inches
  object  
  dev.off()
  
} 


#gopath(dmrs_ANA,all.cpg=rownames(betas[,!is.na(ss_clean$ANA_dom)]),n=10,ann=ann)->pat
gopath <- function(object,all.cpg=NULL,n=20,ann=NULL){
  library(future)
  library(future.apply)
  # object[[1]]
  #require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  requireNamespace(c("S4Vectors","future","future.apply","missMethyl","data.table"))
  # future::plan("multisession")
  plan(tweak(multisession, workers = 6))
  cont<-unique(object$Contrast)
  pathways <- future.apply::future_lapply(cont, function(x){
    
    results.ranges<- object[object$Contrast == x,]
    gst_go <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
                                                    all.cpg = all.cpg,
                                                    collection = "GO", 
                                                    array.type = "EPIC",
                                                    anno = ann,
                                                    sig.genes = T)
    g<-missMethyl::topGSA(gst, n=n)
    g$Contrast <- x
    g$method <- "GO"
    g
    },
    error = function(e) {
      data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO")
    }
    )
    gst_prom <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
                                                      all.cpg = all.cpg,
                                                      collection = "GO", 
                                                      array.type = "EPIC",
                                                      genomic.features=c("TSS200","TSS1500","1stExon"),
                                                      anno = ann,
                                                      sig.genes = T)
    g<-missMethyl::topGSA(gst, n=n)
    g$Contrast <- x
    g$method <- "GO_prom"
    g
    },
    
    error = function(e) {
      data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="GO_prom")
    }
    )
    gst_kegg <- tryCatch({gst <- missMethyl::goregion(results.ranges, 
                                                      all.cpg = all.cpg,
                                                      collection = "KEGG", 
                                                      array.type = "EPIC",
                                                      anno = ann,
                                                      sig.genes = T)
    g<-missMethyl::topGSA(gst, n=n)
    g$Contrast <- x
    g$method <- "KEGG"
    g$ONTOLOGY <- NA
    g$TERM<-g$Description
    g$Description<-NULL
    g
    },
    
    error = function(e) {
      data.frame(ONTOLOGY=NA,TERM=NA,N=0,DE=0,P.DE=1,FDR=1,SigGenesInSet=NA,Contrast=x,method="KEGG")
    }
    )
    result <- data.table::rbindlist(list(gst_go,gst_prom,gst_kegg),fill=T)
    return(result)
  },future.packages = c("S4Vectors","GenomicRanges"))
  
  
  return(data.table::rbindlist(pathways))
  
}

#' For each contrast extract a set of DMPs and add gene annotation and methylation values
#'
#'
#' @title Extract DMPs, annotation and methylation difference for each contrast
#'
#' @return data.table
#' @author Izar de Villasante
#' @export
#' @import minfi
#' @import data.table
#' @import limma
#' @param beta_normalized normalized betavalues, as produce by minfi::getBeta(grSet_noob)),
#'  where colnames(beta_normalized) == metadata$sample_Name
#' @param ContrastsDM list of contrasts as returned by limma::makeContrasts()
#' which will pass to limma topTable as input
#' @param mDiff absolute mean methylation difference between groups to filter by
#' @param ann annotation dataset from manifest with metadata such as gene info,
#' CGI, RefGene, etc. see topTable genelist arg.
#' @param writeOut save result as .csv default = TRUE.
#' @param writedir
#'
#' @inheritParams limma::topTable
#' @examples
#'
#' betas<-readRDS("data/beta_noob.rds")
#' fit<-readRDS("data/fit2.rds")
#' ann<-readRDS("data/ann.rds")
#' DMPann <- DMPextr(fit = fit,                       # linear contrast model
#'                   ContrastsDM = ContrastsDM,          # contrasts
#'                   p.value = 0.01,                      # filter significantly different probes
#'                   beta_normalized = beta_noob,        # extract mean group betas
#'                   mDiff = 0.5,                        # select mean methylation differences
#'                   ann = ann,                          # annotate positions (CGI, RefGene, etc)
#'                   writeOut = FALSE                    # write output to file


DMPextr <- function(
    fit, ContrastsDM=colnames(fit$contrasts), p.value, beta_normalized, mDiff, ann=NULL,
    writedir = "analysis/DMP_", writeOut = TRUE,ncores=NULL){
  require(foreach)
  require(bigstatsr)
  require(data.table)
  if(is.null(ann)){
    dt_epic<-as.data.table(ann<-IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other,keep.rownames="ProbeID")
    dt_450k<-as.data.table(ann<-IlluminaHumanMethylation450kanno.ilmn12.hg19::Other,keep.rownames="ProbeID")   
    ann<- merge(dt_epic,dt_450k,all=T)  
  }
  ann<-data.table::as.data.table(ann,keep.rownames = "ProbeID")
  
  data.table::setkey(ann,"ProbeID") 
  ann<-ann[rownames(beta_normalized),]
  ann[,rn:=1:.N]
  betas <- bigstatsr::FBM(NROW(beta_normalized),NCOL(beta_normalized));betas[] <- as.matrix(beta_normalized)
  
  if(is.null(ncores))ncores=length(ContrastsDM)
  
  cl<- parallel::makeCluster(ncores,outfile="",useXDR=F,type = "FORK")
  parallel::clusterEvalQ(cl,{
    requireNamespace(c("limma","data.table"))
    data.table::setDTthreads(0)
  })
  
  
  doParallel::registerDoParallel(cl)
  
  message("Processing ", length(ContrastsDM)," contrasts. Using ",ncores," cores.")
  res<-foreach::foreach(i=itertools::isplitIndices(length(ContrastsDM), chunks=ncores),#isplitIndices(1400,chunks=ncores),
                        .combine='rbind',
                        # .multicombine = F,
                        .packages = c("limma","data.table"),
                        # .export = c("guessArrayTypes",".default.epic.annotation"),
                        .inorder=F,
                        .errorhandling = "pass"
  )%dopar%{
    DMP_1 <- limma::topTable(fit,
                             num = Inf,
                             coef = i ,
                             genelist = ann,
                             p.value = p.value  # = p.adj
    )
    
    if(nrow(DMP_1)<2){
      warning(paste("No DMP found for contrast:", ContrastsDM[i], sep=" "))
      dt<-ann[0,]
      dt[,c("logFC","AveExpr","t","P.Value","adj.P.Val", "B"):=numeric()]
      dt$Type=character(length = 0L)
      dt$Contrast=character(length=0L)
      dt$diff_meanMeth=numeric(length = 0L)
      
      return(dt)
      
    }
    else{
      DMP_1$Type <- "Hyper"
      DMP_1$Type[which(DMP_1$logFC < 0)] <- "Hypo"  # invert direction (see above0)
      #
      # DMP_1 <- DMP_1[ , c("chr", "pos", "strand", "Name",
      #                     "Type", "P.Value" , "adj.P.Val" ,
      #                     "Islands_Name", "Relation_to_Island",
      #                     "UCSC_RefGene_Name", "UCSC_RefGene_Accession","UCSC_RefGene_Group",
      #                     "Phantom4_Enhancers","Phantom5_Enhancers","X450k_Enhancer",
      #                     "Regulatory_Feature_Name", "Regulatory_Feature_Group",
      #                     "GencodeBasicV12_NAME","GencodeBasicV12_Accession", "GencodeBasicV12_Group",
      #                     "GencodeCompV12_NAME", "GencodeCompV12_Accession", "GencodeCompV12_Group"
      # )]
      DMP_1$Contrast  <- ContrastsDM[i]
      
      # Extract the methylation values for the respective DMP and samples
      # Get contrast i:
      c1 <-fit$contrasts[,ContrastsDM[i]]
      c1 <-c1[c1!=0]
      # Variable names in contrast i:
      vars1 <-names(c1)
      # design matrix for contrasts:
      design <- fit$design[,vars1]
      design <- t(design) * c1
      
      # betas
      
      # b <-betas[rownames(beta_normalized) %in% DMP_1$ProbeID,]
      b <- betas[DMP_1$rn,]
      # diff  _mean:
      DMP_1$diff_meanMeth<-rowSums(apply(design,1,function(x) apply(b%*%diag(x),1,function(y)mean(y[y!=0]))))
      
      # filter for absolute methylation difference
      DMP_1 <- DMP_1[which(abs(DMP_1$diff_meanMeth)>= mDiff), ]
      
      # save DMPs
      
      
      # write output file
      if(writeOut == TRUE){
        
        cat(paste("writing analysis/DMP_", ContrastsDM[i], ".csv\n",sep =""))
        data.table::fwrite(DMP_1, file = paste(writedir, ContrastsDM[i], ".csv", sep =""))
      }
      
      return(DMP_1)
    }
  }
  
  parallel::stopCluster(cl)
  return(res)
}

find_dmrs<-function(object=NULL, betas=NULL, model, fdr = 0.05, p.value = "fdr", bcutoff = 0.3, min.cpg=5, ncores=NULL){
  # FDR threshold used to define DMRS is indexed at the rate of that of DMPs. This rate is defined at fdr. Should only use fdr not p.value
  require(DMRcate)
  require(S4Vectors)
  require(GenomicRanges)
  contrasts <- colnames(model$contrasts)
  require(foreach)
  require(bigstatsr)
  require(data.table)
  if(!is.null(betas)& is.null(object)){
    object <- DMPextr(
      fit = model,                         # Toptable & stats
      ContrastsDM = colnames(model$contrasts),
      beta_normalized = betas,
      p.value = 0.95,
      mDiff = 0.01,
      # ann = ann,
      writeOut = F)
    # }else if(!is.null(object)) {
    #   
    # 
  }
  
  
  conts <- colnames(model$contrasts)
  if(is.null(ncores))ncores<-  min(RcppParallel::defaultNumThreads(),124,22*length(conts))
  
  cl<- parallel::makeCluster(ncores,outfile="",useXDR=F,type = "FORK")
  parallel::clusterEvalQ(cl,{
    requireNamespace(c("limma","data.table","DMRcate","S4Vectors"))
    
  })
  
  dm_threads <- max(ncores%/%length(conts),1)
  doParallel::registerDoParallel(cl)
  
  message("Processing ", length(conts)," contrasts. Using ",ncores," cores.")
  results<-foreach::foreach(i=conts,#isplitIndices(1400,chunks=ncores),
                            .combine='rbind',
                            # .multicombine = F,
                            # .packages = c("limma","data.table","DMRcate"),
                            # .export = c("guessArrayTypes",".default.epic.annotation"),
                            .inorder=F,
                            .errorhandling = "pass"
  )%dopar%{
    # i}
    # dmps<-object
    #   
    object<-object[object$Contrast==i,]
    # remove chromosomes with 1 or less dmps:
    chromosomes <- names(which(table(object$chr)> 1))
    object<-object[object$chr %in% chromosomes,]
    if(nrow(object)<1){
      return(data.table( seqnames=character(),start=numeric(),end=numeric(),width=numeric(),strand=character(), no.cpgs=integer(), min_smoothed_fdr=numeric(),
                         Stouffer=numeric(), HMFDR=numeric(), Fisher=numeric(), maxdiff=numeric(), meandiff=numeric(),
                         overlapping.genes=character(), Contrast=character()))
    }else{
      annotated <- GenomicRanges::GRanges(as.character(object$chr), IRanges(object$pos, object$pos), stat = object$t,
                                          diff = object$logFC, ind.fdr = object$adj.P.Val,
                                          is.sig = object$adj.P.Val < 0.05,Contrasts=i)
    }
    names(annotated) <-  rownames(object)
    annotated <- sort(annotated)
    myAnnotation <- new("CpGannotated", ranges=annotated)
    
    out <- tryCatch(
      {
        # Just to highlight: if you want to use more than one
        # R expression in the "try" part then you'll have to
        # use curly brackets.
        # 'tryCatch()' will return the last evaluated expression
        # in case the "try" part was completed successfully
        
        # message("This is the 'try' part")
        DMRs <- DMRcate::dmrcate(myAnnotation,
                                 pcutoff = "fdr",
                                 betacutoff= bcutoff,
                                 min.cpgs =min.cpg,
                                 # mc.cores=chunks,
                                 C=2
        )
        
      },
      error=function(cond) {
        
        message(cond)
        # Choose a return value in case of error
        new("DMResults", coord=character(), no.cpgs=integer(),
            min_smoothed_fdr=numeric(),Stouffer=numeric(), HMFDR=numeric(),
            Fisher=numeric(), maxdiff=numeric(), meandiff=numeric())->out
      },
      finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>'
        message(paste("Processed contrast:", i))
        message("Next.")
        
      }
    )
    
    if(!identical(out@no.cpgs, integer(0))){
      results.ranges <- DMRcate::extractRanges(out)
      results.ranges$Contrast = i
      # data.frame(Contrast = x, results.ranges)
    }else{
      return(data.table( seqnames=character(),start=numeric(),end=numeric(),width=numeric(),strand=character(), no.cpgs=integer(), min_smoothed_fdr=numeric(),
                         Stouffer=numeric(), HMFDR=numeric(), Fisher=numeric(), maxdiff=numeric(), meandiff=numeric(),
                         overlapping.genes=character(), Contrast=character()))
    }
    
    return(data.table::as.data.table(results.ranges))
    
  }
  
  # 
  # while (is.list(results)) {
  #   results<-suppressWarnings(do.call("c",results[sapply(results,function(x)class(x)=="GRanges")]))
  # }
  # if (is.null(results)){
  #   opts=c("no.cpgs", "min_smoothed_fdr","Stouffer", "HMFDR", "Fisher", "maxdiff","meandiff","overlapping.genes", "Contrast")
  #   tp=list(integer(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),character(),character())   
  #   functionXXX<-function(opts,tp){
  #     library(GenomicRanges)
  #     names(tp)=opts
  #     df=do.call(data.frame,tp)
  #     gr=GenomicRanges::GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL)) 
  #     GenomicRanges::mcols(gr)<-df
  #     gr
  #   }
  #   results <- functionXXX(opts,tp)
  # } 
  return(results)
}

### DMPS:

summary_dmps <- function(DMPextr_object,dir="./results/dmps/",name="raw",write=F){
  require(data.table)
  dt<-data.table::as.data.table(DMPextr_object)
  dt_summary<-dt[,.(
    "Hyper" = sum(Type=="Hyper"),Hypo=sum(Type=="Hypo")
    # "min/max" = paste(range(round(diff_meanMeth,2)),collapse=" / "),
    # "p<0.05" = sum(adj.P.Val<0.05),
    # "min5" = paste(head(round(sort(diff_meanMeth),2),5),collapse=";"),
    # "max5" = paste(tail(round(sort(diff_meanMeth),2),5),collapse=";"),
    
  ),by=c("Contrast")
  ]
  if(write){
    data.table::fwrite(dt,paste0(dir,name,"_dmp_raw.csv.gz"))
    data.table::fwrite(dt_summary,paste0(dir,name,"_dmp_summary.txt"))
  }
  
  return(dt_summary)
}


filter_dmps <- function(dmps, p.value = 0.01, mDiff = 0.3, s=F){
  require(data.table)
  dmps<-data.table::as.data.table(dmps)
  dmps_f <- dmps[adj.P.Val <= p.value & 
                   abs(diff_meanMeth) >= mDiff, ]
  dmps_s<-summary_dmps(dmps_f)
  
  if(s==T){out<-dmps_s[,.SD,.SDcols=c("Hypo","Hyper","Contrast")]}else{out<-dmps_f}
  return(out)
}

apply_filter_dmps<-function(dmps,dev = "png",p.value=seq(0.00,0.1,.01),
                            mDiff=seq(0.15,0.5,.05),path="analysis/intermediate/dmps/"){
  if(length(dmps)<1|is.null(dmps)){warning("no dmps available...")}else{
    
    require(ggplot2)
    dir.create(path)
    params<-expand.grid(p.value,mDiff,T)
    # names(params)<-c("p.val","mDiff","s")
    # params|>purrr::map(function(x)filter_dmps(dmps,x[1],x[2],x[3]))
    # p2<-split(params, seq(nrow(params)))
    # res2<-p2 |> purrr::map(\(x)filter_dmps(dmps,unlist(x[1]),unlist(x[2]),unlist(x[3])))
    
    res1<-with(params, Map(function(a,b,c) {
      dt<-filter_dmps(dmps, a,b,c)
      dt$p.val<-a
      dt$mDiff<-b
      dt
    }, Var1, Var2, Var3))
    pdata<-rbindlist(res1)
    pdata[,All:=Hypo+Hyper]
    pd<-melt(pdata,measure.vars=c("Hyper","Hypo","All"))
    
    pd$variable<-factor(pd$variable)
    if(NROW(pd)>0 & length(levels(pd$variable))>1){
      plt_list<-list()
      plt_list[["dmp_params"]] <-ggplot2::ggplot(data=pd,aes(x=mDiff,y=value,group=p.val,color=p.val))+#,color=Contrast,group=varaible))+
        geom_line(aes(linetype=factor(p.val)))+
        facet_grid(variable~Contrast)#,margin="variable")
      lapply(1:length(plt_list),function(x)
        ggsave(plot = plt_list[[x]],
               filename = paste0(names(plt_list[x]),".",dev),
               path = path,
               device = dev))
    }
  }
}

#DMRs
summary_dmrs <- function(dmrs,path="/results/dmrs/",write=T){
  dmrs[,Type:=ifelse(meandiff>0,"Hyper","Hypo")]
  dmrs.l<-dmrs[,list(Hyper.DMRS=sum(Type=="Hyper"),Hypo.DMRS=sum(Type=="Hypo")),by=c("Contrast")]
  genes.l<-dmrs[,list(Hyper.Genes=length(unique(unlist(strsplit(overlapping.genes[Type=="Hyper"],",")))),Hypo.Genes=length(unique(unlist(strsplit(overlapping.genes[Type=="Hypo"],","))))),by=c("Contrast")]
  summary<-merge(dmrs.l,genes.l)
  if(write)data.table::fwrite(summary,path)
  return(summary)
}

summarize<-function(dmps,dmrs,path="results/"){
  dir.create(path)
  sdmps<-summary_dmps(dmps,write = F)
  sdmrs<-summary_dmrs(dmrs,write=F)
  s<-merge(sdmps,sdmrs,all=T)
  data.table::fwrite(s,paste0(path,"/summary.csv"))
  return(s)
}

get_pathways <- function(dmrs, res.folder="results/pathways/", cols=c("term_size","query_size","intersection_size"), pval=0.05, topN=50,savefile=FALSE){
  require(data.table)
  data.table::setDT(dmrs)
  if(nrow(dmrs)<1){
    warning("no genes supplied")
    varnames=c("query", "significant", "p_value", cols, "precision", "recall", "term_id", "source", "term_name", "effective_domain_size", "source_order", "parents", "FDR", "TERM", "Contrast")
    results <- setNames(data.table(matrix(nrow = 0, ncol = length(varnames))), varnames)
    
  }else{
    
    if(savefile==TRUE)suppressWarnings(dir.create(res.folder))
    full_pathways <- pathway(dmrs[,.SD,.SDcols=c("overlapping.genes","Contrast")], cols=cols, pval = pval, topN=topN ) 
    if(savefile==TRUE)data.table::fwrite(full_pathways,paste0(res.folder,"/full_pathway.csv"))
    hyper_pathways <- pathway(dmrs[meandiff>0,.SD,.SDcols=c("overlapping.genes","Contrast")], cols=cols, pval = pval, topN=topN)
    if(savefile==TRUE)data.table::fwrite(hyper_pathways,paste0(res.folder,"/hyper_pathway.csv"))
    hypo_pathways <- pathway(dmrs[meandiff<0,.SD,.SDcols=c("overlapping.genes","Contrast")], cols=cols, pval = pval, topN=topN)
    if(savefile==TRUE)data.table::fwrite(hypo_pathways,paste0(res.folder,"/hypo_pathway.csv"))
    hypo_pathways$status="hypo"
    hyper_pathways$status="hyper"
    full_pathways$status="both"
    results <- rbind(hypo_pathways,hyper_pathways,full_pathways)
  }
  return(results)
}

pathway <- function(dmrs, path="results/pathways.csv", cols=c("term_size","query_size","intersection_size"), pval=0.05, topN=50,savefile=FALSE){
  require(gprofiler2)
  require(data.table)
  pathways <- lapply(unique(dmrs$Contrast),function(cont){
    p<-gprofiler2::gost(signif = T ,unique(dmrs[Contrast==cont,overlapping.genes],user_threshold=pval))[[1]]
    dth <- data.table::as.data.table(p)
    pat <- data.table()
    # Adapt output naming to convention similar to missmethyl 
    if(length(dth) > 0){
      
      dth[,FDR:= p_value]
      dth[,TERM:= term_name]
      dth[,source:= factor(source)]
      dth[,Contrast:= cont]
      pat <- path_results(
        pathway = dth,
        topN = topN,
        method = "source",
        path = path,
        cols = cols,
        pval=pval,
        savefile=savefile
      )
    }  
    return(pat)
  })
  result<-do.call("rbind",pathways)
  return(result)
}



#' Generate pathway results in an excel sheet. Filters for 
#' at least topN pathways for each group and all terms with FDR <= FDR. 
#' @title generate results excel sheet for pathway analysis
#' @param pathway data.table with pathway analysis results compatible
#'with gopath function output values
#' @param topN numeric value. Minimum number of terms for each group
#' @param method .SDcols argument to group by. Could be a single value or
#' a vector of columns in pathway object to group by. contrast, 
#' @param FDR FDR threshols to filter out.
#' @return plots 
#' @author izar de Villasante
#' @export
#'
path_results<-function(pathway,topN=50,method="method",pval=0.05,path="results/pathways.csv",cols=NULL,savefile=FALSE){
  require(data.table)
  pathway$method<-pathway[[method]]
  data.table::setorder(pathway,method,FDR)
  sig_idx <- pathway[,.I[FDR < pval]  ,by=method]$V1
  head_idx<-pathway[,.I[1:min(..topN,.N)],by=c(method,"Contrast")]$V1
  res<-pathway[base::union(sig_idx,head_idx),]
  # res[,TERM:=ifelse(FDR<pval,paste0("*** ",TERM," ***"),TERM)]
  results<-res[,.SD,.SDcols=c("Contrast","FDR",cols,"TERM","method")]
  data.table::setorder(results,Contrast,method,FDR)
  if(savefile==TRUE) data.table::fwrite(results,path)
  return(results)
}


