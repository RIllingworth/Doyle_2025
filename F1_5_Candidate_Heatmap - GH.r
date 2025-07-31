###########################################################################
#####Make heatmap of selected candidate genes##############################

#Parameters and libraries
library(limma)
library(RColorBrewer)
library(colorRamps)
library(Biobase)
library(ggplot2)
library(rgl)
options(scipen=20)
##Determine package versions
packageVersion("limma") # ‘3.60.2’
packageVersion("RColorBrewer") # ‘1.1.3’
packageVersion("colorRamps") # ‘2.3.4’
packageVersion("Biobase") # ‘2.64.0’
packageVersion("ggplot2") #‘3.5.1’
packageVersion("rgl") #‘1.3.17’

####Heatmap_Paramaters
#Heatmap name
name<-"Marker_NSC_HM"
#colours
HMcols<-c("blue","white","red")
#z score range for heatmap
zrange<-c(-3,0,3) 
##marker genes to plot
markers<-list(
c("Sox2","Sox9","Pax6","Six3","Emx2","Olig2","Foxg1","Epha5","Zic3","Nr2e1","Dmrta2"),
c("Twist1","Tbx15","Eya1","Hlx","Fgf1","Tbx18","Mylk"))           
##marker list groups - must be one per gene group
names(markers)<-c("Neurodevelopment","Other_Lineage")
##columns containing expresison data to plot in the correct order
Columns<-c(54,56,57,53,55,58)

#Set rootdirectory
Root_Directory<-"Your_Root_Dir"



# # # # # # # # # # # # Make heatmap and associated gene lists 

#Set root directory
setwd(Root_Directory)
setwd("R_Images")
load("1_Summary_Expr_Data.RData")

#Summary Expresison Data from R_image following script 2
Data<-Summary_Expression

###Create and Set output directory
setwd(Root_Directory)
setwd("Figure_1/Analaysis_Output/Heatmap&GeneLists")


##check that gene names are in table and output coordinates for them
for(i in 1:length(markers)){
    tmp<-which(match(markers[[i]],Data$Gene_Symbol)>0)
    markers[[i]]<-markers[[i]][tmp]}

####number of gene sets and coordinates within the list of unique markers
nclust<-length(markers)
genes<-unlist(markers)

##identify subset of data that includes all of your genes of interest
data<-which(match(Data$Gene_Symbol,genes)>0)

##Subset only the columns containing the log2 expression values         
Subset<-Data[data,Columns]
SSnames<-Data$Gene_Symbol[data]

##collapse to a single gene per subset - e.g. average transcript counts - requires to un-log the data (NB log2)
len<-length(genes)
SS<-matrix(ncol=ncol(Subset),nrow=len)

for(i in 1:len){
      c<-which(match(SSnames,genes[i])>0)
      SS[i,]<-(colMeans(2^Subset[c,]))
      }

colnames(SS)<-colnames(Subset)

#make trimmed and ordered gene lists and names 
Subset<-SS
subsetNames<-genes

###### Covert to z-scores
mean_data <- rowMeans(Subset)
sd_data <- apply(Subset, 1, sd)
Subset <- (Subset - mean_data) / sd_data

##Make table of final gene list ordered as the final heatmap will be
GenesOut<-cbind(subsetNames[length(subsetNames):1],Subset[length(subsetNames):1,])

#determine the cooridnates of each in the data table to plat the summary profiles 
Coordinates<-list()
for(i in 1:length(markers)){
      if(i == 1){
      Coordinates[[i]]<-c(1:length(markers[[i]]))}else{
      Coordinates[[i]]<-(1:length(markers[[i]])) + (
      Coordinates[[(i-1)]][length(Coordinates[[(i-1)]])]) 
      }}

nclust<-length(markers)

###write out a raw table file ordered as the heatmap
write.table(GenesOut,paste(name,"_OrderedGenes.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


#Write out cluster information
##write gene lists for each cluster
##N.B. writes in reverse order to orientate them as in the heatmap
for(k in 1:nclust){
    write.table(subsetNames[Coordinates[[k]]],paste("Cluster_",nclust-(k-1),".txt",sep=""),quote=FALSE,row.names=TRUE,sep="\t")}


##plot gene profiles for each cluster
##numbered from top of the heatmap to the bottom
for(k in 1:nclust){
    set<-unlist(Coordinates[k])
    pdf(paste("Gene Profiles for Cluster ",nclust-(k-1),".pdf",sep=""),height = 5, width = 7,pointsize=12,paper="special")
        plot(as.matrix(Subset[set[1],]),ylim=c(-5,5),ylab="log2 FC",type="l",axes=FALSE,xlim=c(0,ncol(Subset)+1),col="white",xlab="Population",main=paste("Gene Profiles for Cluster ",nclust-(k-1),sep=""))
        axis(1,at=seq(1,ncol(Subset),1),labels=colnames(Subset))
        axis(2,at=c(-5,0,5))
        SummaryStats<-matrix(ncol=ncol(Subset),nrow=5)
          for(i in 1:ncol(Subset)){
          SummaryStats[,i]<-quantile(Subset[set,i], c(.10, .25, .50, .75, .90))}
        #max to min
        poly1<-cbind(c((1:(ncol(SummaryStats))), ((ncol(SummaryStats)):1)),c(SummaryStats[1,],SummaryStats[5,ncol(SummaryStats):1]))
        poly2<-cbind(c((1:(ncol(SummaryStats))), ((ncol(SummaryStats)):1)),c(SummaryStats[2,],SummaryStats[4,ncol(SummaryStats):1]))
        polygon(poly1,col="grey90",border=NA)
        polygon(poly2,col="grey70",border=NA)
        #mean
        points(colMeans(Subset[set,]),type="l",lwd=1,lty=2,col="grey30")
        abline(h=0,lty=2,lwd=1,col="black")
        legend("topleft",legend=paste("n =",length(set),sep=" "),bty="n")
        dev.off()
}


###Plot heatmap of z score
#HM colours 
colfunc <- colorRampPalette(HMcols)
#data range - N.B. values beyond these limits will becapped to the limits
range<-zrange


##add in massd out values 
SS<-Subset
tmp<-which(Subset < range[1], arr.ind = TRUE)
SS[tmp]<-range[1]
tmp<-which(Subset > range[3], arr.ind = TRUE)
SS[tmp]<-range[3]

resolution<-0.1
pdf(paste(name,"_OrderedGenes.pdf",sep=""),width = 9, height = 15)
    par(mar=c(5,1,0.5,4.95))
    image(x=c(1:ncol(SS)), y=c(1:nrow(SS)),
    z=t(SS[,]), axes=FALSE, xlab="",ylab="",
    col=colfunc((length(seq(range[1],range[3],resolution))-1)),  breaks=seq(range[1],range[3],resolution))
    box()
    dev.off()
