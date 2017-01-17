###########################################################################
###############Pathwave for GEO Datasets: 11/23/2015#########################
setwd("~/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/PATHWAVE")

##Step 1: Load in data file (Same as DIRAC)
load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/1021sPTBvsPET.Rdata")
colnames(DATA)<-DATA[1,]

#Generate Expression Matrix (remove identifying information)
Exp.Mat<-as.data.frame(DATA[-c(1:2),])

#Crete dataclass file (Case/Control status)
dataclass=as.factor(DATA[2,])

#Run pathwave

library('PathWave')

##param.value: for threshold, 0.5
##10,000 3 minutes

pwres = pathWave.run(preprocessed.tag="KEGG.hsa", input.exprdata=Exp.Mat, input.sampleclasses=dataclass,
                     param.numperm=10000, param.pvalue.threshold=0.2,
                     param.kegg.only_metabolism=FALSE, param.filter.size=3, verbose=TRUE)

pwres$results.table
