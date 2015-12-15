###########################################################################
###############Pathwave for GEO Datasets: 11/23/2015#########################
setwd("~/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/PATHWAVE")
load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/1021sPTBvsPET.Rdata")
colnames(DATA)<-DATA[1,]
#RAM: Read in the csv file for expression matrix setting first column as rownames and header equals true
#First Column is rownames,
Exp.Mat<-as.data.frame(DATA[-c(1:2),])

#data = read.csv("C:/Users/rharihar/Desktop/CHDI_PROJECT/Multi_tissue_PW/Striatum_final.csv",header = T,row.names =1)

#Write out the tsv file for expression matrix

write.table(Exp.Mat, file='PTBvsPETExpressionMatrix.tsv', quote=FALSE, sep='\t', col.names = NA)

#read in dataclass file
#dataclass = read.csv("C:/Users/rharihar/Desktop/CHDI_PROJECT/Multi_tissue_PW/Striatum_final_class.csv", header = F)
dataclass=as.factor(DATA[2,])

#write out tsv file for classes

write.table(dataclass, file='C:/Users/rharihar/Desktop/CHDI_PROJECT/Multi_tissue_PW/Striatum_final_class.tsv', quote=FALSE, sep='\t')

#Running PathWave

library('PathWave')

##param.value: for threwshold, 0.5
##10,000 3 minutes
#preprocesssed.tag=check out, probably run user

#pwres = pathWave.run(preprocessed.tag="KEGG.hsa", input.exprdata=Exp.Mat, input.sampleclasses=dataclass,
                    # param.numperm=10000, param.pvalue.threshold=0.1,
                     #param.kegg.only_metabolism=TRUE, param.filter.size=3, verbose=TRUE,
                     #output.file.prefix="~/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/PATHWAVE/output")

pwres = pathWave.run(preprocessed.tag="KEGG.hsa", input.exprdata=Exp.Mat, input.sampleclasses=dataclass,
                     param.numperm=10000, param.pvalue.threshold=0.1,
                     param.kegg.only_metabolism=TRUE, param.filter.size=3, verbose=TRUE)

pwres$results.table

write.csv(pwres$results.table,"C:/Users/rharihar/Desktop/CHDI_PROJECT/Multi_tissue_PW/PW_res/Striatum_table.csv")

keggurls <- pathWave.getColorKEGGMapURLs(pwres$results.filtered,
                                         preprocessed.tag="KEGG.hsa",
                                         col=c("green","grey","red"),
                                         col.sign.pattern=c("red","red","green"))

write(keggurls,"results.txt")
