#########################
###Search PubMed for specific genes and their associations with disease e.g. Alzheimer's disease###
##########################
###based on###
### see https://cran.r-project.org/web/packages/reutils/reutils.pdf
#############
library(reutils)
genes <- read.table("list-of-genes.txt", header=FALSE)
for ( i in 1:length(genes[,1])){
   comp <- tryCatch({ ###skip if no papers found
   gene.name <- as.character(genes[i,])
   query <- paste(gene.name,"[all] and  alzheimer disease [mesh]", sep="")
   pmids <- esearch(query, "pubmed", usehistory = TRUE)
   articles <- efetch(pmids)
   titles <- articles$xmlValue("//ArticleTitle")
   write.table(rbind (gene.name,titles),file="list-of-genes-pubmed-Alzheimer.txt", col.names=FALSE, quote=FALSE, append=TRUE)
   }, error=function(e) NULL)
}
