arqxml <- xmlParse(filepath)
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


corrplot(arq_fstmat_mtmt,method="color",type="lower",na.label = " ",tl.col="black",tl.srt=20,order="FPC",cl.lim = c(0,1),p.mat = arq_fstp,sig.level = 0.05,insig="pch",pch.cex = 1,pch.col = "#6b717a")
