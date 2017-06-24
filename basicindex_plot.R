#Load packages
pkg_list = c("ggplot2","cowplot")
lapply(pkg_list,require,character.only=TRUE)


#Take the output from basicindextable function
basicindexplot <- function(input){
	#This function take the output of basicindextable function
	combineddf <- input$dataframe

	#Haplotype diversity bar chart
	combineddf$group <- factor(combineddf$group,
		levels = combineddf$group[order(combineddf$`haplotype diversity`)])

	hd_p <- ggplot(combineddf, aes(x = group, y = `haplotype diversity`))
	haplotypeDiversityPlot <- hd_p + geom_bar(stat = 'identity',fill="#1470e0",color="black",size=1,alpha=0.8) +  geom_errorbar(width=0.3,size=1,aes(ymin=`haplotype diversity`-`haplotype diversity sd`,ymax=`haplotype diversity`+`haplotype diversity sd`)) + 
	coord_flip() + theme_minimal() 

	#Nucleotide diversity bar chart
	combineddf$group <- factor(combineddf$group,
		levels = combineddf$group[order(combineddf$`nucleotide diversity`)])
	nd_p <- ggplot(combineddf, aes(x = group, y = `nucleotide diversity`))
	nucleotideDiversityPlot <- nd_p + geom_bar(stat = 'identity',fill="#f77533",color="black",size=1,alpha=0.8) +  geom_errorbar(width=0.3,size=1,aes(ymin=`nucleotide diversity`-`nucleotide diversity sd`,ymax=`nucleotide diversity`+`nucleotide diversity sd`)) + 
	coord_flip() + theme_minimal() 

	#Scatter plot (hd vs nd)
	scatter_p <- ggplot(combineddf,aes(x = `haplotype diversity`, y = `nucleotide diversity`, color=`haplotype diversity` + `nucleotide diversity`))
	hnScatterPlot <- scatter_p + geom_point(size=3) + geom_text(aes(label=group),vjust=1.5,check_overlap = TRUE) + 
	scale_color_gradient(limits = c(0, 0.001), low = "gray75", high = "#229399") + theme_minimal() + theme(legend.position="none")

	plot_grid(haplotypeDiversityPlot,nucleotideDiversityPlot,hnScatterPlot,labels=c("A","B","C"), nrow=2, ncol=2)

	}

	