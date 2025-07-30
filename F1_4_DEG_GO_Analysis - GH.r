###########################################################################
#####GO analysis of DEGs###################################################
##Note that the changes in annotation version s can impact on the reults

#parametrs and libraries
library(Biobase)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(limma)
options(scipen=20)

##Determine package versions
packageVersion("Biobase") # ‘2.58.0’
packageVersion("clusterProfiler") # ‘4.6.2’
packageVersion("AnnotationDbi") # ‘1.60.2’
packageVersion("org.Mm.eg.db") # ‘3.16.0’
packageVersion("limma") # ‘3.54.2’

#Set root directory
Root_Directory<-"Your_Root_Dir"
setwd(Root_Directory)
setwd("R_Images")
load("1_Summary_Expr_Data.RData")

#Summary Expresison Data from R_image following script 2
Data<-Summary_Expression

###Create and Set output directory
setwd(Root_Directory)
setwd("Figure_1/Analaysis_Output/Heatmap&GeneLists")


#####Identify DEGs
## DEG expresison paramaters
co<-0.01
diff<-1

#########################################################################
#Identify significantly differntilly expressed gene sets
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

up_c<-(unique(c(up1,up6,up_both)))

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

down_c<-(unique(c(down1,down6,down_both)))


#Gene Lists including background
up<-unique(Data$Gene_Symbol[up_c])
down<-unique(Data$Gene_Symbol[down_c])

#Background gene set
GS_All<-unique(Data$Gene_Symbol)

GO_results_up_BP<-enrichGO(gene=up, OrgDb="org.Mm.eg.db",keyType="SYMBOL",universe=GS_All,ont="BP",pvalueCutoff=0.05,pAdjustMethod = "BH")
GO_results_down_BP<-enrichGO(gene=down, OrgDb="org.Mm.eg.db",keyType="SYMBOL",universe=GS_All,ont="BP",pvalueCutoff=0.05,pAdjustMethod = "BH")


##write out plots and text outputs
pdf("GO_results_up_BP.pdf",width=7,height=10)
    print(dotplot(GO_results_up_BP,showCategory=2000,x="p.adjust",color="p.adjust"))
    dev.off()

write.table(as.data.frame(GO_results_up_BP),"GO_results_up_BP.txt",col.names=TRUE,quote=FALSE,sep="\t",row.names=FALSE)

pdf("GO_results_down_BP.pdf",width=7,height=10)
    print(dotplot(GO_results_down_BP,showCategory=2000,x="p.adjust",color="p.adjust"))
    dev.off()

write.table(as.data.frame(GO_results_down_BP),"GO_results_down_BP.txt",col.names=TRUE,quote=FALSE,sep="\t",row.names=FALSE)
