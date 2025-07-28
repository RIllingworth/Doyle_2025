##Expression heatmap and scatter plot for Figure 1 (ex vivo NSCs)
#########################################################################
#parametrs and libraries
options(scipen=20)
library(limma)
library(Biobase)


#set root directory
Root_Directory<-"Your_Root_Dir"

setwd(Root_Directory)
setwd("R_Images")
load("1_Summary_Expr_Data.RData")

#Summary Expresison Data from R_image following script 2
Data<-Summary_Expression

## DEG expresison paramaters
co<-0.01
diff<-1

###Create and Set output directory
setwd(Root_Directory)
setwd("Figure_1/Analaysis_Output")
dir.create("Heatmap&GeneLists")
setwd("Heatmap&GeneLists")


#########################################################################
#Identify significantly differntilly expressed gene sets
#Up in NSC7_1(NSC_1) and NSC7_6(NSC_2) vs respective WTs
up1<-which(
     (Data$logFC_NSC7_1vNSC7_2>= diff)&(Data$adj.P.Val_NSC7_1vNSC7_2<=co)&
     (Data$logFC_NSC7_1vNSC7_9>= diff)&(Data$adj.P.Val_NSC7_1vNSC7_9<=co)&
     (Data$logFC_NSC7_6vNSC7_2> 0)&
     (Data$logFC_NSC7_6vNSC7_9> 0))
 
up6<-which(
     (Data$logFC_NSC7_6vNSC7_2>= diff)&(Data$adj.P.Val_NSC7_6vNSC7_2<=co)&
     (Data$logFC_NSC7_6vNSC7_9>= diff)&(Data$adj.P.Val_NSC7_6vNSC7_9<=co)&
     (Data$logFC_NSC7_1vNSC7_2> 0)&
     (Data$logFC_NSC7_1vNSC7_9> 0))

up_both<-which(
    (Data$logFC_NSC7_1vNSC7_2>= diff)&(Data$adj.P.Val_NSC7_1vNSC7_2<=co)&
    (Data$logFC_NSC7_1vNSC7_9>= diff)&(Data$adj.P.Val_NSC7_1vNSC7_9<=co)&
    (Data$logFC_NSC7_6vNSC7_2>= diff)&(Data$adj.P.Val_NSC7_6vNSC7_2<=co)&
    (Data$logFC_NSC7_6vNSC7_9>= diff)&(Data$adj.P.Val_NSC7_6vNSC7_9<=co))

up<-(unique(c(up1,up6,up_both)))
up1only<-setdiff(up1,up6)
up6only<-setdiff(up6,up1)

#Down in NSC7_1(NSC_1) and NSC7_6(NSC_2) vs respective WTs
down1<-which(
     (Data$logFC_NSC7_1vNSC7_2<= -diff)&(Data$adj.P.Val_NSC7_1vNSC7_2<=co)&
     (Data$logFC_NSC7_1vNSC7_9<= -diff)&(Data$adj.P.Val_NSC7_1vNSC7_9<=co)&
     (Data$logFC_NSC7_6vNSC7_2< 0)&
     (Data$logFC_NSC7_6vNSC7_9< 0))
 
down6<-which(
     (Data$logFC_NSC7_6vNSC7_2<= -diff)&(Data$adj.P.Val_NSC7_6vNSC7_2<=co)&
     (Data$logFC_NSC7_6vNSC7_9<= -diff)&(Data$adj.P.Val_NSC7_6vNSC7_9<=co)&
     (Data$logFC_NSC7_1vNSC7_2< 0)&
     (Data$logFC_NSC7_1vNSC7_9< 0))

down_both<-which(
    (Data$logFC_NSC7_1vNSC7_2<= -diff)&(Data$adj.P.Val_NSC7_1vNSC7_2<=co)&
    (Data$logFC_NSC7_1vNSC7_9<= -diff)&(Data$adj.P.Val_NSC7_1vNSC7_9<=co)&
    (Data$logFC_NSC7_6vNSC7_2<= -diff)&(Data$adj.P.Val_NSC7_6vNSC7_2<=co)&
    (Data$logFC_NSC7_6vNSC7_9<= -diff)&(Data$adj.P.Val_NSC7_6vNSC7_9<=co))

down<-(unique(c(down1,down6,down_both)))
down1only<-setdiff(down1,down6)
down6only<-setdiff(down6,down1)

#write out DEG lists 
write.table(unique(Data[up,1]),"Upreg_NSC7_Refseq.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(unique(Data[down,1]),"Downreg_NSC7_Refseq.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)




#########################################################################
#####write out an expression table for supplemental data
Data$DiffExpression       <- 0
Data$DiffExpression[up]   <- 1
Data$DiffExpression[down] <- -1

#####Add a column for NSC7 averaged expression (log2) 
Data$Log2_NSC7_I53A<-log2(((2^Data$NSC7_1)+(2^Data$NSC7_6))/2)
Data$Log2_NSC7_WT<-log2(((2^Data$NSC7_2)+(2^Data$NSC7_9))/2)

#####Write out an expression table for supplemental data
Sup<-Data[,c(1:7,8:10,14:16,29:31,32:34,23:25,44:46,53:59,60,61)]
write.table(Sup,"Master_NSC_Expr.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)



#########################################################################
#####Make a heatmap of DEGs 
#Add a unique gene and location ID to remove duplicate values
Data$ID<-paste(Data[,2],Data[,3],Data[,4],Data[,5],sep="_")

#Make subset to cluster ordered as indicated
cluster<-(Data[c(down1only,down6only,down_both,up_both,up6only,up1only),c(38,8,14,29,32,44,59)])

#remove duplicated genes
a<-which(duplicated(cluster$ID))
unq<-setdiff(1:nrow(cluster),a)
cluster<-cluster[unq,]
cluster<-cluster[,1:6]

#Name of heatmap
name<-"DiffExpr_LogFC"

#Set up heatmap colours
colfunc <- colorRampPalette(c("blue","white","red"))
#Set heatmap value range (log2 fold expression change)
range<-c(-8,0,8)

#Cap the colour scale to the maximum and minimum range values
a1<-which(cluster <= range[1],arr.ind = TRUE)
cluster[a1]<-range[1]
a2<-which(cluster >= range[3],arr.ind = TRUE)
cluster[a2]<-range[3]

#Set the resolution of the colour scale and make a postscript and tiff version of the heatmap
resolution<-(range[3]-range[1])/100
###run through 2x for plotting for ps and for tiff
for(p in 1:2){
      if( p==1 ){ 
          postscript(paste(name,".ps",sep=""),height = 5, width = 5,pointsize=12,paper="special",horizontal=FALSE)}else
          {tiff(paste(name,".tiff",sep=""),res=300,width = 1500, height = 1500, units = "px", pointsize = 12)}
          par(mar=c(5,0.5,3,3))
          image(x=c(1:ncol(cluster)), y=c(1:nrow(cluster)),
          z=t(cluster[]), axes=FALSE, xlab="",ylab="",
          col=colfunc((length(seq(range[1],range[3],resolution))-1)),  breaks=seq(range[1],range[3],resolution))
          box()
          dev.off()}

#Plot a colour bar
tiff("colbar_ratio.tiff",res=300,width = 1500, height = 1500, units = "px", pointsize = 12)
      plot(100,100,col="white",axes=FALSE,ylab="",xlab="",main="",ylim=c(0,100),xlim=c(0,100))
      rect(rep(20,40),c(20:60),rep(30,40),c(21:61),col=colfunc(41),border=NA)
      text(30,20,labels=range[1],adj=c(-0.5,0),pos=4)
      text(30,60,labels=range[3],adj=c(-0.5,0.5),pos=4)
      dev.off()





#########################################################################
#Make a scatter plot of averaged log2 expression values
#Covert to averaged NSCs 1 and 2 (6_1 and 6_6) expression values, covert back to Log2 and plot as a scatter plot. 
I53A<-log2(((2^Data$NSC7_1)+(2^Data$NSC7_6))/2)
WT<-log2(((2^Data$NSC7_2)+(2^Data$NSC7_9))/2)


###########################pdf
pdf("Expr.pdf",useDingbats=FALSE,width=7,height=7)
    par(mar=c(5.1, 5.1, 4.1, 2.1), bty="l")
    plot(WT,I53A,pch=16,cex=0.1,col="black",xlab="WT Expression (Log2)",ylab="I53A Expression (Log2)",main=name,cex.main=2,cex.lab=2,ylim=c(0,18),xlim=c(0,18),axes=FALSE)
    axis(1,at=c(0,9,18))
    axis(2,at=c(0,9,18))
    points(WT[up],I53A[up],pch=16, col="#d95f02",cex=2)
    points(WT[down],I53A[down],pch=16, col="#7570b3",cex=2)
    abline(0,1,lty=2)
    dev.off()
