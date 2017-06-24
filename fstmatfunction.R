
#Load packages


pkg_list = c("XML","corrplot","magrittr","dplyr")
lapply(pkg_list,require,character.only=TRUE)


ordered_fstmat <- function(filepath,remove0=TRUE,order=TRUE,pvalue=FALSE){
  
  #----Read the XML and extract contents----#
  arqxml <- xmlParse(filepath)
  arq_fstmat <- xpathSApply(arqxml,"//PairFstMat",xmlValue) %>% as.vector()
  arq_fstmat <-strsplit(arq_fstmat,"\n")
  
  #--Clean and transform FstMat to matrix--#
  arq_fstmat_mt <- do.call(cbind,arq_fstmat[1])
  arq_fstmat_mt %<>% subset(arq_fstmat_mt[,1] != "")
  arq_fstmat_mt <- gsub(" + "," ",arq_fstmat_mt)
  arq_fstmat_mt <- strsplit(arq_fstmat_mt," ")
  #--Warning: Matrix with non-consistence columns number--##
  try(arq_fstmat_mt <- do.call(rbind,arq_fstmat_mt))
  
  #--Clean the matrix again--#
  arq_fstmat_mt <- arq_fstmat_mt[3:nrow(arq_fstmat_mt),3:ncol(arq_fstmat_mt)]
  arq_fstmat_mt[upper.tri(arq_fstmat_mt)] <-NA
  
  #The matrix looks good but each cell is in character class and no col and rownames
  arq_fstmat_mt %<>% apply(c(1,2),as.numeric)
  diag(arq_fstmat_mt) <- NA
  arq_fstmat_mt %<>% as.data.frame()
  
  #Let find out the rowname in pairDistPopLabels tag
  arq_rowname <- xpathSApply(arqxml,"//pairDistPopLabels",xmlValue) %>% strsplit("\n")
  arq_rowname <- do.call(cbind,arq_rowname) 
  arq_rowname <- gsub("+\\d+:\t","",arq_rowname)
  arq_rowname <- arq_rowname[5:nrow(arq_rowname)] %>% as.vector()
  
  #Assign name back to the fstmat df
  rownames(arq_fstmat_mt) <- arq_rowname
  colnames(arq_fstmat_mt) <- arq_rowname
  
  #Prepare a full matrix form for corrplot
  arq_fstmat_mtmt <- as.matrix(arq_fstmat_mt)
  arq_fstmat_mtmt[upper.tri(arq_fstmat_mtmt)] <- t(arq_fstmat_mtmt)[upper.tri(arq_fstmat_mtmt)]
  diag(arq_fstmat_mtmt) <- 0
  
  #Remove negative value
  if(remove0 == TRUE){
    arq_fstmat_mtmt[arq_fstmat_mtmt<0] <-0
  }

  #Create pvalue matrix
  if(pvalue == TRUE){

    arq_fstp <- xpathSApply(arqxml,"//PairFstPvalMat",xmlValue) %>% as.vector()
    arq_fstp<-strsplit(arq_fstp,"\n")

    #arq_fstp now is a nested list

    #--Clean and transform fstp to matrix--#
    arq_fstp <- do.call(cbind,arq_fstp[1])
    arq_fstp %<>% subset(arq_fstp[,1] != "")
    arq_fstp <- gsub(" + "," ",arq_fstp)
    arq_fstp <- strsplit(arq_fstp," ")
    try(arq_fstp <- do.call(rbind,arq_fstp))
    arq_fstp <- arq_fstp[2:nrow(arq_fstp),3:ncol(arq_fstp)]
    arq_fstp[upper.tri(arq_fstp)] <-NA
    diag(arq_fstp) <- 0

    #Remove the SD of pvalue
    arq_fstp <- gsub(" ","",arq_fstp)
    arq_fstp <- gsub("\\+.*","",arq_fstp)

    #Full matrix
    arq_fstp[upper.tri(arq_fstp)] <- t(arq_fstp)[upper.tri(arq_fstp)]
    rownames(arq_fstp) <- arq_rowname
    colnames(arq_fstp) <- arq_rowname
    arq_fstp %<>% apply(c(1,2),as.numeric)
      }
        

  #Corrplot
  
  if(order==TRUE){

    if(remove0 == TRUE){
      fstmat_plot <- corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black",tl.srt=10,order="FPC",cl.lim = c(0,1))
    }
    else{
      fstmat_plot <- corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black", tl.srt=10,order="FPC")
    }

    if(pvalue==TRUE){
      fstmat_plot <-  corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black",tl.srt=10,order="FPC",cl.lim = c(0,1),
      p.mat = arq_fstp,sig.level = 0.05,insig="pch",pch.cex = 1,pch.col = "#6b717a")
    }
  }
  
  if(order==FALSE){
    if(remove0 == TRUE){
      fstmat_plot <- corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black",tl.srt=10,order="original",cl.lim = c(0,1))
    }
    else{
      fstmat_plot <- corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black", tl.srt=10,order="original")
    }

    if(pvalue==TRUE){
      fstmat_plot <-  corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black",tl.srt=10,order="FPC",cl.lim = c(0,1),
      p.mat = arq_fstp,sig.level = 0.05,insig="pch",pch.cex = 1,pch.col = "#6b717a")
    }
  }
  
  if(order=="alphabet"){
    if(remove0 == TRUE){
      fstmat_plot <- corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black",tl.srt=10,order="alphabet",cl.lim = c(0,1))
    }
    else{
      fstmat_plot <- corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black", tl.srt=10,order="alphabet")
    }
    if(pvalue==TRUE){
      fstmat_plot <-  corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black",tl.srt=10,order="FPC",cl.lim = c(0,1),
      p.mat = arq_fstp,sig.level = 0.05,insig="pch",pch.cex = 1,pch.col = "#6b717a")
    }
  }
   return(list("plot"=fstmat_plot,"matrix"=arq_fstmat_mtmt))
}

#remove0: turn all negative value to zero
#order: group the groups using first principal component
#pvalue: default is FALSE, if pvalue is TRUE, cross the value with p-value < 0.05