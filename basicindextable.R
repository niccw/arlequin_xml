#This function extract these data from the Arliquin XML output:
#	Group Name: groupname
#	No. of gene copies: genecopiesno (equals to group size if it is a mtDNA marker(default))
#	No. of sequence : sequenceno
#	No. of polymorphic site: polysiteno
#	Genetic diversity : genediversity
#	Nucleotide diversity : nucleodiversity

# groupfullname argument takes a list of full name and add it to the dataframe, default with no full group name

# Retrun a DT table and a dataframe object for further analysis


#Load packages
pkg_list =c("XML","tidyverse","magrittr","corrplot","tidyr","DT","tibble","stringr")

lapply(pkg_list,require,character.only=TRUE)

basictable <- function(filepath,groupfullname=FALSE,mtDNA=TRUE){
  
  #Get all data nodes
  arqxml <- xmlParse(filepath)
  alldata <- xpathSApply(arqxml,"//data",xmlValue) 
  alldatamt <- as.matrix(alldata)
  
  #--Extract information block (no. of gene copies/ no. of sequence/ no. of polymorphic sites/ gene diversity)--#
  #Comment attached with the single string version for reference
  
  datablk1 <- alldatamt[grep(".*No. of gene copies.*:.*",alldatamt)]  #grep(".*No. of gene copies.*:.*\n",alldatamt)
  
  #No. of gene copies
  genecopiesno <- sapply(datablk1,str_extract,".*No. of gene copies.*\n")    #str_extract(teststring,".*No. of gene copies.*\n")
  genecopiesno <-strsplit(genecopiesno,"\n")
  genecopiesno <- do.call(rbind,genecopiesno) #Yeah! matrix!
  genecopiesno <- str_extract(genecopiesno,"\\d+") %>% as.vector() 

  if(mtDNA==FALSE){
  	genecopiesno <- genecopiesno/2
  }
  
  #No. of sequence
  sequenceno <- sapply(datablk1,str_extract,".*No. of sequences.*\n")    
  sequenceno <-strsplit(sequenceno,"\n")
  sequenceno <- do.call(rbind,sequenceno) 
  sequenceno <- str_extract(sequenceno,"\\d+") %>% as.vector() 
  
  #No. of polymorphic sites
  polysiteno <- sapply(datablk1,str_extract,".*No. of polymorphic sites.*\n")    
  polysiteno <-strsplit(polysiteno,"\n")
  polysiteno <- do.call(rbind,polysiteno) 
  polysiteno <- str_extract(polysiteno,"\\d+") %>% as.vector() 
  
  #Gene (Haplotype) Diversity
  genediversity <- sapply(datablk1,str_extract,".*Gene diversity.*\n")    
  genediversity <-strsplit(genediversity,"\n")
  genediversity <- do.call(rbind,genediversity) 
  genediversity <- str_extract(genediversity,"\\d\\.\\d+.*\\d") %>% as.vector() 
  genediversity <- gsub(" ","",genediversity)
  genediversity <- gsub("%+-%","+/-",genediversity) #Convert the plus minus to plotmath format
  
  
  #--Extract information block(Nucleotide diversity)--#
  datablk2 <- alldatamt[grep(".*Mean number of pairwise differences .*:.*",alldatamt)]
  #Nucleotide Diversity
  nucleodiversity <- sapply(datablk2,str_extract,".*Nucleotide diversity.*\n")    
  nucleodiversity <-strsplit(nucleodiversity,"\n")
  nucleodiversity <- do.call(rbind,nucleodiversity) 
  nucleodiversity <- str_extract(nucleodiversity,"\\d\\.\\d+.*\\d") %>% as.vector() 
  nucleodiversity <- gsub(" ","",nucleodiversity)
  nucleodiversity <- gsub("%+-%","+/-",nucleodiversity)
  
  #--group list--#
  groupname <- xpathSApply(arqxml,"//pairDistPopLabels",xmlValue) %>% strsplit("\n")
  groupname <- do.call(cbind,groupname) 
  groupname <- gsub("+\\d+:\t","",groupname)
  groupname <- groupname[5:nrow(groupname)] %>% as.vector()
  
  
  #Checkk if all length are equal--!!#
  if(!all.equal(length(groupname),length(genecopiesno),length(sequenceno),length(polysiteno),length(genediversity),length(nucleodiversity))){
    warning("Dfferent list length.")
  }
  
  

  #----Combine all together to df----#
  #Function control inpit groupfullname
  if(groupfullname==FALSE){
    combineddf <- data.frame(groupname,genecopiesno,sequenceno,polysiteno,genediversity,nucleodiversity)
    colnames(combineddf) <- c("group","n","no. of allele","no. of polymorphic sites","haplotype diversity","nucleotide diversity")
  }
  else{
    if(length(groupname)==length(groupfullname)){
      combineddf <- data.frame(groupname,groupfullname,genecopiesno,sequenceno,polysiteno,genediversity,nucleodiversity)
      colnames(combineddf) <- c("group","location","n","no. of allele","no. of polymorphic sites","haplotype diversity","nucleotide diversity")
    }
    else{
      stop("Length of groupfullname is not identical with the no of group name.")
    }
  }



  #Clean the df for futher use in R
  hdtempmt <- do.call(rbind,strsplit(as.character(combineddf$`haplotype diversity`),"+/-",fixed = TRUE))
  ndtempmt <- do.call(rbind,strsplit(as.character(combineddf$`nucleotide diversity`),"+/-",fixed = TRUE))

  combineddf %<>% mutate(`haplotype diversity` = hdtempmt[,1], `haplotype diversity sd` = hdtempmt[,2],
  `nucleotide diversity` = ndtempmt[,1], `nucleotide diversity sd` = ndtempmt[,2])

  combineddf$`haplotype diversity` <- as.numeric(combineddf$`haplotype diversity`)
  combineddf$`nucleotide diversity` <- as.numeric(combineddf$`nucleotide diversity`)
  combineddf$`haplotype diversity sd` <- as.numeric(combineddf$`haplotype diversity sd`)
  combineddf$`nucleotide diversity sd` <- as.numeric(combineddf$`nucleotide diversity sd`)

  combineddf %<>% tbl_df()
  #--------#
  

  DT <- datatable(combineddf,rownames=FALSE, extensions = 'Buttons',options=list(dom="Bftrip",buttons=c('copy','csv','pdf')))
  
  return(list("DTtable"=DT,"dataframe"=combineddf))
  
}  


    