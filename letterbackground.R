#####----Check the update of the html----####

#Get the text from html using XML package

library(XML)
doc.html = htmlTreeParse("http://xxx") #check if need to use internal node

#check the tag insde the table
doc.txt = unlist(xpathApply(doc.html,'tag'),xmlValue)

#Process the string to an organized df

#(Optional) Extract the gist contain (PDF), process it to get the zoning and proposed use 

#Compare the df with the last version, extract the new item for analysis

#Save the processed df with date


####----Check with the db----####

#Check if the site name in any ecological important list

#Check the previous related letter in the same site

#(If have the zoning and proposed info) Check similar cases

#Compare the lists and find out which are the most similar cases as ref

#?Antway to get the gov data in list/db form?#

####----(Additional)Update the db from different csv----####

####----(Ultimate)Analysis the content of related lettes (n>3), extract common sentences, construct a backbone of letter----####

