library(ChAMP)
library(minfi)
library(FactoMineR)
library(factoextra)
####load the package needed

 
dir<-"C:/Users/fangyu/Desktop/DNAm/meth"
myLoad <- champ.load(dir,arraytype = "EPIC",detPcut = 0.05)
####read in the IDAT file in the meth directory
 
CpG.GUI(CpG = rownames(myLoad$beta),arraytype = "EPIC")
###check the CpG information
 
QC.GUI(beta=myLoad1$beta,arraytype = "EPIC")
###visualize the quality of the probes
 
myNorm <- champ.norm(arraytype="EPIC")
##normalization
   
   
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))
####combat comes from sva package.This is a way to remove potential biology variables that may have impact on high-through output experiment


   
myLoad$pd$Sample_Group[myLoad$pd$Sample_Group=="ME/CFS case"]="MECFS"      
myDMP <- champ.DMP(beta = myCombat,pheno=myLoad$pd$Sample_Group)
####change the group name to valid R group name and do the differential methylation probes analysis between case and control group

       
DMP_inter<-intersect(rownames(myDMP$MECFScase_to_Healthy),rownames(myDMP1$MECFScase_to_Healthy))
DMP_int<-myDMP$MECFScase_to_Healthy[DMP_inter,]
myDMP_1=myDMP
myDMP_1$MECFScase_to_Healthy=DMP_int
##get DMP intersext
 
##venn diagraom
library(VennDiagram)
venn.plot<-venn.diagram(x=list("Slides_only"=1:28507,"Slides_and_Sex"=25:(40583+24)),"SLIDES_SEX.png",fill=c("red","blue"))
CpG.GUI(CpG = rownames(myDMP$MECFScase_to_Healthy),arraytype = "EPIC")

myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="Bumphunter")
dmr_chr=myDMR$BumphunterDMR$seqnames
names<-table(dmr_chr)
num<-c(14,8,18,7,2,3,7,11,16,3,3,9,5,1,9,8,8,8,51,13,8)
chr<-names(names)
dt<-data.frame(num,chr)
data_L<-data.frame(Number=myDMR$BumphunterDMR$L,length=myDMR$BumphunterDMR$width)
data_L$DMR<-"All"
ggplot(data_L, aes(x=DMR, y=length,color=DMR)) + geom_boxplot()+scale_color_manual(values=c("#E69F00"))
ggplot(data_L, aes(x=DMR, y=Number,color=DMR)) + geom_boxplot()+scale_color_manual(values=c("#ff0000"))
###myDMR(differential methylation region) chrom distribution pie plot
myDMR_6<-myDMR$BumphunterDMR[which(myDMR$BumphunterDMR$seqnas=="chr6"),]
###DMR from chrom 6,extract information


myDMP_filter=myDMP
myDMP_filter$MECFS_to_Healthy<-myDMP$MECFS_to_Healthy[which(abs(myDMP$MECFS_to_Healthy$deltaBeta)>0.05),]   
filter<-myDMP_filter$MECFScase_to_Healthy
filter$gene<-as.character(filter$gene)
filter_2<-filter$P.Value
names(filter_2)=rownames(filter)
group_list=myLoad$pd


####functions for volcano and MA plot
visual_champ_DEM <- function(myLoad,myDMP,group_list,pro='test'){
  beta.m=myLoad$beta
  champDiff=myDMP[[1]] 
  head(champDiff)  
  colnames(champDiff)
  ## for volcano 
  if(T){
    nrDEG=champDiff
    head(nrDEG)
    attach(nrDEG)
    plot(logFC,-log10(P.Value))
    library(ggpubr)
    df=nrDEG
    df$logp= -log10(P.Value)
    ggscatter(df, x = "logFC", y = "logp",size=0.5)
    df$color='stable'
    for(i in 1:length(df$logFC))
    {
      if(df$P.Value[i]<=0.05)
      {
        if(df$logFC[i]> 0.2) {df$color[i]='up'}
        if(df$logFC[i]< -0.2) {df$color[i]='down'}
        if(df$logFC[i]> 0 & df$logFC[i]< 0.2) {df$color[i]='slight up'}
        if(df$logFC[i]< 0 & df$logFC[i]>= -0.2) {df$color[i]='slight down'}
      }
    }
#plot the down and up regulation according to logFC.        

df$name=rownames(df)
head(df)
ggscatter(df, x = "logFC", y = "logp",size=0.5,color = 'color')
ggscatter(df, x = "logFC", y = "logp", color = "color",size = 1,
          label = "name", repel = T,font.label = c(12,'black'),
          label.select =  head(rownames(df)[df$logp != 'stable']),
          palette = c("blue","red","red","blue") )
ggsave(paste0(pro,'_volcano.png'))
ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
table(df$p_c)
ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size = 1, 
          palette = c("green", "red", "black","blue") )
ggsave(paste0(pro,'_MA.png'))


}

###differential methylation probes heatmap  
  beta.m=myCombat            ##set the beta value,raw matrix is myload,norm is myNorm and 
  champDiff=myDMP_filter[[1]]

  ## for heatmap 
  dat=beta.m

  table(group_list)
  deg=champDiff
  x=deg$logFC 
  names(x)=rownames(deg)
  cg=c(names(head(sort(x),100)),
       names(tail(sort(x),100)))   ##sort the cg by logFC value       
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
  #n=t(scale(t(dat[cg,])))
  #n[n>2]=2
  #n[n< -2]= -2

  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(n)
  pheatmap(dat[cg,],show_colnames =F,
           show_rownames = F,
           cluster_cols = T, fontsize = 5,
           annotation_col=ac,
           filename = paste0("11",'_heatmap_top200_DEG_scale.png'))  
  pheatmap(dat[cg,],show_colnames =F,
           show_rownames = F,
           cluster_cols = T, fontsize = 5,
           annotation_col=ac,
           filename = paste0("11",'_heatmap_top200_DEG_raw.png'))  
  
}
###plot top and tail(most and least) 100 significant DMPs
 

###intersect
pone=readxl::read_xlsx("pone.xlsx")
a<-intersect(pone$`Probe Name`,rownames(myDMP$MECFS_to_Healthy))
pone<-data.frame(pone)
bb<-pone[which(pone$Probe.Name %in% a),]
bb1<-myDMP$MECFS_to_Healthy[a,]
bb=-bb$Mean.Beta.Difference
bb1=bb1$deltaBeta
t.test(bb,bb1,paired=TRUE)

DMP.GUI(DMP=myDMP[[1]],beta=myCombat,pheno=myLoad$pd$Sample_Group)
##visualize DMP 
DMR.GUI(DMR = myDMR,arraytype = "EPIC",pheno = myLoad$pd$Sample_Group)
##visualize DMR


pie.labels=names(table(myDMP_filter$MECFS_to_Healthy$feature))
pie.number=c(309,466,964,4467,3008,1781,625)
df <- data.frame("number"=pie.number, "labels"=pie.labels)
ggplot(df, aes(x = "", y = number, fill = labels)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  # geom_bar will generate a barplot
  coord_polar("y", start = 0)+
  # coord_polar will create a pie chart based on barplot
  scale_fill_brewer(palette="Dark2") +
  # change the color
  theme_void()
               

beta.m=myCombat
group_list=myLoad$pd$Sample_Group
dim(beta.m) 
###quality check plots:PCA,variable heatmap       
if(T){
  
  
  dat=t(beta.m)
  
  library("FactoMineR")#needed to plot PCA plot
  library("factoextra")  
  # this step will not be fast
  dat.pca <- PCA(dat , graph = FALSE) 
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave('all_samples_PCA.png')
  
  dat=beta.m
  
  cg=names(tail(sort(apply(dat,1,sd)),500) ##rank the beta matrix by sd
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
  #n=t(scale(t(dat[cg,]))) 
  n[n>2]=2 
  n[n< -2]= -2
  

  ac=data.frame(group=group_list)
  rownames(ac)=colnames(n)  
  pheatmap(n,show_colnames =F,show_rownames = F, fontsize = 5,
           annotation_col=ac,filename = 'heatmap_top500_sd.png')  ##heatmap plot showing top 500 variable probes
  dev.off()
  
  exprSet=beta.m
  pheatmap::pheatmap(cor(exprSet)) 
  
  colD=data.frame(group_list=group_list)
  rownames(colD)=colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet),
                     annotation_col = colD,
                     show_rownames = F,fontsize = 5,
                     filename = 'cor_all.png')   ###correlation plot
  dev.off() 
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  #M=cor(log2(exprSet+1)) 
  M=cor(exprSet)
  pheatmap::pheatmap(M,annotation_col = colD)
  pheatmap::pheatmap(M,
                     show_rownames = F,
                     annotation_col = colD,fontsize = 5, 
                     filename = 'cor_top500.png')
  dev.off() 
  
}
