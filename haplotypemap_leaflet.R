pkg_list = c("readr","stringr","plyr","magrittr")

#----Clean the output from DNASP for the df----#

haptable <- function(filepath,name){

	dnaspout <- read_lines(filepath)
	dnaspout %<>% as.matrix()

	seqno <- dnaspout[grep("Number of haplotypes",dnaspout)]
	seqno <- str_extract(seqno,"\\d.*") %>% as.numeric()

	startline <- 20 + seqno + 1  #haplotype 1 number strart from line20 in DNASP output
	dnaspout <- dnaspout[startline:nrow(dnaspout),] %>% as.matrix()
	dnaspout <- dnaspout[-nrow(dnaspout),] %>% as.matrix()
	#Seperate to three columns
	haplist <- str_extract(dnaspout,"Hap_\\d*") %>% as.vector()
	haplist <- haplist[!is.na(haplist)]
	hapno <- str_extract(dnaspout,":\\d*.*\\[") %>% as.matrix() %>% str_extract("\\d+") %>% as.vector()
	hapno <- hapno[is.na(hapno)]

	hapdist <- str_extract(dnaspout,"\\[.*\\]")
	hapdist <- gsub("\\[","",hapdist)
	hapdist <- gsub("\\]","",hapdist)
	hapdist <- gsub("\\d+","",hapdist)
	hapdist <- hapdist[!is.na(hapdist)]
	hapdist %<>% strsplit(" ")

	#Count number of group in each haplotype 
	#hapdist is now a list

	df_list <- list()
	hapnocount <- 1
	for (group in hapdist){
		hapnocount_name <- paste0("Hap_",hapnocount)

		haptable <- table(group) %>% as.data.frame()
		try(colnames(haptable) <- c("group",hapnocount_name))
		df_list[[hapnocount]] <- haptable

		hapnocount = hapnocount + 1
	}

	#Bind all the df for each haplotype to one df
	finalhapdist <- join_all(df_list,type="full")
	assign(paste0(name,"haptable"),finalhapdist,envir = .GlobalEnv)
}

#----trying----#



rownames(finalhapdist) <-  finalhapdist$group
finalhapdist %<>% subset(select=-group) #drop the group columns
finalhapdist %<>% apply(c(1,2),as.numeric)
finalhapdist %<>% t()
finalhapdist %<>% as.data.frame()

#Pie chart
pie_HTT <- ggplot(na.omit(finalhapdist),aes(x=factor(1),fill=factor(HTT))) + geom_bar(width=1) +
coord_polar(theta='y')


# Ideas: Generate a platte for Hap, check if there's any function for this
# For each group, plot the pie using a function which specify the color used
# Also need GPS list which the order need to be identical to the pie charts

# And then let see how to plot each of the pie chart on leaflet