#This is an analysis script for detecting conditionally rare taxa in a temporal microbial community dataset.  
#Written by A. Shade 30 May 2013/02 Dec 2013, to accompany the manuscript: "Conditionally rare taxa disproportionately contribute to temporal changes in microbial diversity."  
#This script comes with no warranty.
#Questions?  shade.ashley@gmail.com

#####
#16 Oct 2014 bug fix.  ALS.  MaxRel filter was updated.  Also added option:  can discover of CRT based on MaxRel calculated from dataset with all OTUs OR dataset with only non-singleton OTUs.
####

#####
#What does the script do?
#This script will print the proportion of conditionally rare taxa detected in the dataset in the R console.  
#It will also output a file of the OTU IDs, and, if provided, the taxonomic assignments of those OTUs, for the conditionally rare taxa.

#The script allows the user to define thresholds of the coefficient of bimodality (b_thresh, default = 0.90), and the relative abundance maximum (abund_thresh, default = 0.005). 


#####
#What are the input files?
#The input file for this script is: An OTU (taxa) table, with samples in columns and taxa in rows.  
#The first row should include column names.  The first column should have taxa (OTU) IDs.  
#The first cell (row 1, col 1) should be empty.  
#It is optional that the last column contains taxonomic assignments of each OTU.

#####
#How do I use the script? 
#Step 1.
#load required R packages: vegan, TSA. 
library(vegan)
library(TSA)

#Step 2.
#Place the input file and script in the same working directory to run this script.  Change the working directory in R to match where the files have been placed.

#Step 3.
#Load the necessary functions into your R workspace, contained in a separate file, "rare_fncs.R" 
#source("../rare_fncs.R")

#####
#16 Oct 2014 bug fix.  ALS.  MaxRel filter was updated.  Also added option:  can discover of CRT based on MaxRel calculated from dataset with all OTUs OR dataset with only non-singleton OTUs.
####

#function to make relative abundance table - load into workspace
makeRFtable.f=function(data){
  cSum1<-colSums(data)
  
  #define an empty matrix to put your RF values into
  newdata<-matrix(0,dim(data)[1], dim(data)[2])
  
  #Assign the same column and row names to the new matrix.
  colnames(newdata)<-colnames(data)
  rownames(newdata)<-rownames(data)
  
  #Each cell will be divided by the column sum.
  for (i in 1:length(data)){
    newdata[,i] <- data[,i]/cSum1[i]
  }
  
  return(newdata)
}
###


#Rare to Prevalent OTUs
SimpleRareToPrev.f=function(otu_fp,abund_thresh = 0.005, abund_thresh_ALL=FALSE,b_thresh = 0.90, rdp_lastcol=TRUE){
  
  #Read in files
  otu=read.table(otu_fp, header=TRUE, check.names=FALSE, row.names=1, sep="\t")
  
  #If provided, remove rdp ids and save
  if(rdp_lastcol==TRUE){
    rdp=otu[,ncol(otu)]
    otu=otu[,-ncol(otu)]
  }
  
  
  #remove empty rows
  tmp=otu[rowSums(otu)>0,]
  no.otus=nrow(tmp)
  
  if(rdp_lastcol==TRUE){
    rdp2=rdp[rowSums(otu)>0]
  }
  
  #Remove singletons
  otu.nosigs=tmp[rowSums(tmp)>1,]
  
  if(rdp_lastcol==TRUE){
    rdp3=rdp2[rowSums(tmp)>1]
  }else{
    rdp3=NULL
  }
  
  
  #how many are left after singletons?
  no.sigs=nrow(otu.nosigs)
  
  
  #Make a rel abundance table - with the full dataset
  otu.rel=makeRFtable.f(tmp)
  #reduce the rel. abundance table to omit singletons
  otu.rel2=otu.rel[is.element(row.names(otu.rel), row.names(otu.nosigs)),]
  
  #loop for each OTU to calculate Coefficent of bimodality
  #For efficiency, loops through singleton-omitted dataset (singletons will not have rare-to-prevalent dynamics)
  out=NULL
  for(j in 1:nrow(otu.nosigs)){
    
    x=as.numeric(otu.nosigs[j,])
    k=kurtosis(x)
    s=skewness(x)
    
    #calculate the coefficient of bimodality for each OTU
    b=(1+(s^2))/(k+3)
    
    #determine whether OTU max and median relative abundance (based on full dataset, otu.rel)
    x2=as.numeric(otu.rel[j,])
    mx.rel.all=max(x2)
    med.rel.all=median(x2)
    
    x3=as.numeric(otu.rel2[j,])
    mx.rel.ns=max(x3)
    med.rel.ns=median(x3)
    
    out=rbind(out,c(row.names(otu.nosigs)[j],b,mx.rel.all, med.rel.all, mx.rel.ns, med.rel.ns))
  }
  
  #print(dim(out))
  #print(head(out))
  out=as.data.frame(out)
  colnames(out)=c("OTUID","CoefficientOfBimodality","MaxRel_All", "MedianRel_All", "MaxRel_NoSingletons", "MedianRel_NoSingletons")
  
  if(rdp_lastcol==TRUE){
    out=cbind(out, rdp3)
    colnames(out)[7]="TaxonomicAssignment"
  }
  
  #Filter 1: for at least one rel. abundance greater than abund_thresh (default = 0.005).  The default uses the whole dataset MaxRel (abund_thresh_ALL=TRUE), another option is the singleton-removed dataset.
  if(abund_thresh_ALL==TRUE){
    at="ALL"
    out.filter=out[as.numeric(as.vector(out[,"MaxRel_All"])) >= abund_thresh,]
    print(dim(out.filter))
  }else{
    at="NOSIG"
    out.filter=out[as.numeric(as.vector(out[,"MaxRel_NoSingletons"])) >= abund_thresh,]
    print(dim(out.filter))
  }
  
  #Filter 2: for coefficient of bimodality greater than b_thresh (default = 0.90)
  out.filter=out.filter[as.numeric(as.vector(out.filter[,"CoefficientOfBimodality"])) >= b_thresh,]
  print(dim(out.filter))
  
  write.table(out.filter, paste("ResultsFile_ConditionallyRareOTUID_", abund_thresh, "_", b_thresh, "_", at, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
  
  
  print("No. conditionally rare OTUs")
  print(nrow(out.filter))
  
  print("No. total OTUs")
  print(no.otus)
  
  print("Proportion conditional rare / total OTUs")
  print(nrow(out.filter)/no.otus)
  
  print("No singleton OTUs")
  print(no.sigs)
  
  print("Proportion conditionally rare / non-singletonOTUs")
  print(nrow(out.filter)/no.sigs)
  return(out.filter)
  
} 


#Step 4.  
#Change the options below to match your dataset.  The options are:  
#otu_fp - type the the full name of your dataset file, including the extension
#abund_thresh -  Change the maximum abundance threshold, if desired. Defaults to 0.005
#abund_thresh_ALL - Use TRUE if you want to use the full dataset (ALL OTUs) to calculate relative abundances. Use FALSE if you want to use the non-singleton (filtered) dataset to calculate relative abundances.  Default is FALSE.
#b_thresh - Change the coefficient of bimodality threshold, if desired.  Defaults to 0.90
#rdp_lastcol - Use TRUE if the last column of the dataset contains the taxonomic assignments of OTUs, use FALSE if not
#Then,to run the script, copy and paste the command into the R console:

(cond_rare <- SimpleRareToPrev.f(otu_fp="/Users/kiristern/Documents/GitHub/PBIN/data/ASVs_counts_copy.tsv",abund_thresh=0.005, abund_thresh_ALL=FALSE,b_thresh=0.90, rdp_lastcol=FALSE))
length(low_abundance(ASV_count, detection=0.5/100))

#select all the conditionally rare ASVs from phyloseq object
crare_virps <- subset_taxa(viral_physeq, (rownames(tax_table(viral_physeq)) %in% cond_rare$OTUID))
nonrare_virps <- subset_taxa(viral_physeq, !(rownames(tax_table(viral_physeq)) %in% cond_rare$OTUID))
tcrare_virps <- t(crare_virps)

#store data in timeseries object
ts_rare <- ts(tcrare_virps)

ts_rare$sample <- row.names(ts_rare)

ts_condrare <- reshape2::melt(ts_rare, id="sample")
ggplot(ts_condrare) +
  geom_line(aes(x=Var1, y=value, group=Var2, color=Var2))


#plot ASV over time
ASV1732 <- t(as.data.frame(virps["ASV_1732",]))
ASV2813 <- t(as.data.frame(virps["ASV_2813",]))
ASV841 <- t(as.data.frame(virps["ASV_841",]))

plot(ASV1732, pch="o", col="blue", lty=1, type="o", ylim=c(0,120))
#overlay on initial plot
points(ASV2813, col="red", pch="*")
lines(ASV2813, col="red", lty=2)
points(ASV841, col="green", pch="+")
lines(ASV841, col="green", lty=3)


