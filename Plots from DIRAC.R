###DIRAC: DIFFERENTIALLY REGULATAED GENES
library(biomaRt)
library(vioplot)

setwd("~/Documents/Project 1. Preterm Birth/AnalysisJan2016")
metacovar<-read.csv("IlluminaMetaCovariate.csv")
rownames(metacovar)<-metacovar$SampleIdentifier

#Track 1: PET vs. PTB####
load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/ProcessingJan2016/Option2ComBat.Rdata")
data2<-t(exprs(eset_COMBAT))

###Genes from biocarta pathway of interest####
#path.File<-"/Users/alisonpaquette/Documents/General Price Lab notes/GSEA Files/GSEAentrez2biocarta.gmt" #need to get from msig database, this is name of pathway and all the genes.  Can change to whatever you want
path.File<-"/Users/alisonpaquette/Documents/General Price Lab notes/GSEA Files/GSEAentrez2kegg.gmt"
readPathway<-function(path.File, exp.data){
  pathways<-readLines(path.File)
  path<-strsplit(pathways, "\t")

  return.pathways = list()
  for ( i in 1:length(path)){
    namesDIR<-as.matrix(path[[i]])

    if (length(path[[i]]) > 5 && length(which(rownames(exp.data) %in% namesDIR[,1])) > 5){
      return.pathways <-append(return.pathways, list(path[[i]]))
    }
  }

  names(return.pathways)<-sapply(return.pathways,'[',1)
  final.pathways = lapply(return.pathways,'[',-(1:2))
  return(final.pathways)
}
Tpathway.list<-readPathway(path.File, exprs(eset_COMBAT))
GENES1<-Tpathway.list$KEGG_AMINO_SUGAR_AND_NUCLEOTIDE_SUGAR_METABOLISM
GENES<-intersect(colnames(data2),GENES1)
data1<-data2[,GENES]

data<-merge(data1,metacovar,by='row.names',all=F)

##Entrez Ids to Gnes####
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "entrezgene", attributes= c("hgnc_symbol","entrezgene"),values=GENES,mart= mart)
G_list$entrezgene<-as.character(G_list$entrezgene)
rownames(G_list)<-as.character(G_list$entrezgene)


###PLOT####
PET<-subset(data,PET.status=="Y")
rownames(PET)<-PET[,1]
PET<-PET[,c(2:42)]

NOPET<-subset(data,PET.status=="N")
rownames(NOPET)<-NOPET[,1]
NOPET<-NOPET[,c(2:42)]

GENES==colnames(PET)
colnames(PET)==colnames(NOPET)


TT<-NULL
for(i in 1:42)
{
  x<- t.test(PET[,i],NOPET[,i])
  TT[[i]]<-x$p.value
}
TT

PET2<-PET[,c(1,6,7,13,14,16,19,20,21,23,31,33,35,38)]
NOPET2<-NOPET[,c(1,6,7,13,14,16,19,20,21,23,31,33,35,38)]
GENES2<-GENES[c(1,6,7,13,14,16,19,20,21,23,31,33,35,38)]
par(mfrow=c(3,5),mai = c(0.4, 0.4, 0.4, 0.4))
for(i in 1:14)
{
vioplot(PET2[,i],NOPET2[,i],col="azure3",names=c("PE"," No PE"))
title(G_list$hgnc_symbol[G_list$entrezgene==GENES2[i]])
}

TT<-NULL
for(i in 1:42)
{
 x<- t.test(PET[,i],NOPET[,i])
 TT[[i]]<-x$p.value
}
TT

#Track 2: PT vs T####
load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/ProcessingJan2016/Option3ComBat.Rdata")
data2<-t(exprs(eset_COMBAT))

###Genes from biocarta pathway of interest####
path.File<-"/Users/alisonpaquette/Documents/General Price Lab notes/GSEA Files/GSEAentrez2biocarta.gmt" #need to get from msig database, this is name of pathway and all the genes.  Can change to whatever you want
readPathway<-function(path.File, exp.data){
  pathways<-readLines(path.File)
  path<-strsplit(pathways, "\t")

  return.pathways = list()
  for ( i in 1:length(path)){
    namesDIR<-as.matrix(path[[i]])

    if (length(path[[i]]) > 5 && length(which(rownames(exp.data) %in% namesDIR[,1])) > 5){
      return.pathways <-append(return.pathways, list(path[[i]]))
    }
  }

  names(return.pathways)<-sapply(return.pathways,'[',1)
  final.pathways = lapply(return.pathways,'[',-(1:2))
  return(final.pathways)
}
Tpathway.list<-readPathway(path.File, exprs(eset_COMBAT))
GENES<-Tpathway.list$BIOCARTA_EGFR_SMRTE_PATHWAY
GENES<-GENES[-9]
data1<-data2[,GENES]

data<-merge(data1,metacovar,by='row.names',all=F)

##Entrez Ids to Gnes####
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "entrezgene", attributes= c("hgnc_symbol","entrezgene"),values=GENES,mart= mart)
G_list$entrezgene<-as.character(G_list$entrezgene)
rownames(G_list)<-as.character(G_list$entrezgene)


###PLOT####
PT<-subset(data,PretermGroup=="PT")
rownames(PT)<-PT[,1]
PT<-PT[,c(2:11)]

TERM<-subset(data,PretermGroup=="T")
rownames(TERM)<-TERM[,1]
TERM<-TERM[,c(2:11)]

GENES==colnames(PT)
colnames(PT)==colnames(TERM)

par(mfrow=c(2,5),mai = c(0.4, 0.4, 0.4, 0.4))
for(i in 1:10)
{
  vioplot(PT[,i],TERM[,i],col="azure3",names=c("PT","TERM"))
  title(G_list$hgnc_symbol[G_list$entrezgene==GENES[i]])
}

TT<-NULL
for(i in 1:10)
{
  x<- t.test(PT[,i],TERM[,i])
  TT[[i]]<-x$p.value
}
TT

###PRETERM VS TERM W/PREECLAMPSIA####
load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/ProcessingJan2016/Option4ComBat.Rdata")
data2<-t(exprs(eset_COMBAT))

###Genes from biocarta pathway of interest####
path.File<-"/Users/alisonpaquette/Documents/General Price Lab notes/GSEA Files/GSEAentrez2biocarta.gmt" #need to get from msig database, this is name of pathway and all the genes.  Can change to whatever you want
readPathway<-function(path.File, exp.data){
  pathways<-readLines(path.File)
  path<-strsplit(pathways, "\t")

  return.pathways = list()
  for ( i in 1:length(path)){
    namesDIR<-as.matrix(path[[i]])

    if (length(path[[i]]) > 5 && length(which(rownames(exp.data) %in% namesDIR[,1])) > 5){
      return.pathways <-append(return.pathways, list(path[[i]]))
    }
  }

  names(return.pathways)<-sapply(return.pathways,'[',1)
  final.pathways = lapply(return.pathways,'[',-(1:2))
  return(final.pathways)
}
Tpathway.list<-readPathway(path.File, exprs(eset_COMBAT))
GENES<-Tpathway.list$BIOCARTA_MCALPAIN_PATHWAY
GENES<-GENES[-c(4,15,18)]
data1<-data2[,GENES]

data<-merge(data1,metacovar,by='row.names',all=F)

##Entrez Ids to Gnes####
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "entrezgene", attributes= c("hgnc_symbol","entrezgene"),values=GENES,mart= mart)
G_list$entrezgene<-as.character(G_list$entrezgene)
rownames(G_list)<-as.character(G_list$entrezgene)


###PLOT####
PT<-subset(data,PretermGroup=="PT")
rownames(PT)<-PT[,1]
PT<-PT[,c(2:23)]

TERM<-subset(data,PretermGroup=="T")
rownames(TERM)<-TERM[,1]
TERM<-TERM[,c(2:23)]

GENES==colnames(PT)
colnames(PT)==colnames(TERM)

par(mfrow=c(2,5),mai = c(0.4, 0.4, 0.4, 0.4))
for(i in 13:22)
{
  vioplot(PT[,i],TERM[,i],col="azure3",names=c("PT","TERM"))
  title(G_list$hgnc_symbol[G_list$entrezgene==GENES[i]])
}

TT<-NULL
for(i in 13:22)
{
  x<- t.test(PT[,i],TERM[,i])
  TT[[i]]<-x$p.value
}
TT
