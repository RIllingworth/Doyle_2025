#########################################################################
###Generate the limma input files #######################################
#########################################################################

#Uses the un-normalised HOMER quantification of the raw NGS data for all samples 
#Quantifications were performed for exons on the same strand as the gene orientations
#Quant command:
#analyzeRepeats.pl /exports/cmvm/eddie/scs/groups/rilling2-lab/gene_anno/mm39/RefSeq_UCSC-refGene_mm39.gtf none -count exons -d \ [list of tag HOMER directories] -strand + -noadj > ./Quant/I53_NSCs_RefS_mm39_Exons_PosStr.txt

#### Set global paramaters
Root_Directory<-"Your_Root_Dir"
Quantitation_File<-"I53_NSCs_RefS_mm39_Exons_PosStr.txt"
Gene_Annotation_File<-"RefSeq_UCSC-refGene_mm39_inc_name.txt"
#Offset value to add to all genes / transcripts per RPK value
offset<-1




###Run Script#############################################################

setwd(Root_Directory)
setwd("Figure_1/Quant_Files")
Exon<-read.delim(Quantitation_File)

##Tidy table headers 
names<-colnames(Exon)
names[1]<-"Gene_ID"
names<-sub("..TagDir.","",names)
names<-sub("\\.\\..*","",names)
colnames(Exon)<-names

##add specific ID to remove exact duplicate regions 
Exon$MID<-paste(Exon[,1],Exon[,2],Exon[,3],Exon[,4],Exon[,5],sep="_")
Data<-Exon

#trim to required data columns
Data<-Data[,c(2:6,1,9:26)]

###convert NAs into 0s
for (i in c(7:ncol(Data))){
a<-match(Data[,i],NA)
b<-which(a>0)
if(length(b) > 0){Data[b,i]<-0}else{NULL}}

#make the value per kb + offset (RPK)
a<-c(7:ncol(Data))
Data[,a]<-(Data[,a]/(Data$Length/1000))+offset

###annotate including gene symbol
Anno<-read.delim(paste(Root_Directory,"/Annotations/",Gene_Annotation_File,sep=""))
##trim to unique Refseq IDs
a<-which(duplicated(Anno[,1]))
unq<-setdiff((1:nrow(Anno)),a)
Anno<-Anno[unq,]

#Annotate quant files with Gene Symbol From Annotation Files
Data<-merge(Data,Anno,by.x="Gene_ID",by.y="X.name")
Data<-Data[,c(1:6,29,7:24)]

##columns containing desired data
a<-c(8:ncol(Data))

#write out reads per kb exon data
setwd(Root_Directory)
setwd("Figure_1/Analaysis_Output")
dir.create("Files_For_LIMMA")
setwd("Files_For_LIMMA")

for(i in a[1]:a[length(a)]){
tmp<-Data[,c(1:7,i)]
tmp_names<-names(tmp)
tmp_names[length(tmp_names)]<-"Expression"
colnames(tmp)<-sub("\\.x","",tmp_names)
write.table(tmp,paste(sub("\\.x","",colnames(Data)[i]),"Exon.txt",sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")}
