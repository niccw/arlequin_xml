pkg_list <- c("tidyverse","magrittr","XML","igraph","stringr","tibble","visNetwork")
lapply(pkg_list,require,character.only=TRUE)



visnetworkMSN <- function(arqxmlpath,dnaspoutpath,name){
  
  
  arqxml <- xmlParse(arqxmlpath)
  alldata <- xpathApply(arqxml,"//data",xmlValue) %>% as.matrix()
  
  
  ####====edge table====####
  allMSNlist <- alldata[grep("Connection length",alldata),]
  if (length(allMSNlist)==1){
    lastMSN <- allMSNlist %>% as.character()
  } else{
    lastMSN <- allMSNlist[grep("ALL_MST|All_MST",allMSNlist)] %>% as.character() #The last MSN contain all info
  }
  #Clean up the matrix
  lastMSN <- strsplit(lastMSN,"\n")
  lastMSN <- do.call(cbind,lastMSN)
  #Record down the relative position of the data wanted
  startindex <- grep("OTU 1",lastMSN)+2
  endindex <- grep("NEXUS notation for MST",lastMSN)-3
  lastMSN <- lastMSN[startindex:endindex,] %>% as.matrix()
  lastMSN %<>% apply(c(1,2),as.character)
  lastMSN <- strsplit(lastMSN,"\\s+")
  edgemt <- do.call(rbind,lastMSN)
  edgemt <- edgemt[,-1]
  colnames(edgemt) <-c("from","to","length")
  edgemt[,1] <- gsub("Hap_","",edgemt[,1])
  edgemt[,2] <- gsub("Hap_","",edgemt[,2])
  edgemt %<>% tbl_df()
  edgemt[,1:3] <- sapply(edgemt[,1:3],as.numeric) #edge df (from,to,length)
  edgemt %<>% mutate(value=1/length)
  edgemt %<>% mutate(label=length)
  edgedf <- edgemt
  
  #add value for the edges too!!!
  
  ####====node table (from DNAsp out)====####
  
  dnaspout <- read_lines(dnaspoutpath)
  dnaspout %<>% as.matrix()
  
  seqno <- dnaspout[grep("Number of haplotypes",dnaspout)]
  seqno <- str_extract(seqno,"\\d.*") %>% as.numeric()
  startline <- 20 #haplotype 1 number strart from line20 in DNASP output
  endline <- 20 + seqno -1
  
  dnaspout <- dnaspout[startline:endline,] %>% as.matrix()
  nodesize <- str_extract(dnaspout,"\\:.*\\[")
  nodesize <- gsub("\\D","",nodesize)
  nodeID <- c(seq_along(nodesize))
  nodelabel <- paste0("Hap",nodeID)
  
  nodedf <- cbind(nodeID,nodesize,nodelabel)
  
  nodedf %<>% tbl_df()
  nodedf[,1:2] <- sapply(nodedf[,1:2],as.numeric)
  colnames(nodedf) <- c("id","value","label")
  
  vnetwork <- visNetwork(nodedf,edgedf) %>%   visOptions(highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T)
  vnetwork
  visSave(vnetwork,paste0(name,"visnetwork.html"))
  
}

