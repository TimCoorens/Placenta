#Placenta figure 3

#Figure 3c
#tc16

options(stringsAsFactors = F)
patients=read.table("bulk_patients.txt")[,1]

library(gplots)
library(RColorBrewer)
col=colorRampPalette(brewer.pal(9,"Blues"))(100)
col_breaks = seq(0,1,0.01)

Genotype_all=read.csv("Extended_data_table5.csv") #read in supplementary table

asymmetry=c()
for(patient in patients){
 if(grepl("PD45",patient)) UC_sample=paste0(patient,"f")
    if(patient%in%c("PD45559","PD45590","PD45555")) UC_sample=paste0(patient,"g")
    if(grepl("PD42",patient)) UC_sample=paste0("PD",as.numeric(substr(patient,3,7))-1,"b")
    
    Genotype=Genotype_all[,grepl(patient,colnames(Genotype_all))|colnames(Genotype_all)==UC_sample]
    Genotype=Genotype[rowSums(Genotype=="-")==0,]
    first_mut=rownames(Genotype)[which.max(as.numeric(Genotype[,UC_sample]))]
    asymmetry=rbind(asymmetry,c(patient,first_mut,Genotype[first_mut,UC_sample],binom.test(x = NV[first_mut,UC_sample],n=NR[first_mut,UC_sample])$conf.int[1:2]))
  }
  print(patient)
}

asymmetry=as.data.frame(asymmetry)
asymmetry[asymmetry[,1]=="PD45560",3:5]=0.5*as.numeric(asymmetry[asymmetry[,1]=="PD45560",3:5])
asymmetry=asymmetry[order(asymmetry[,3],decreasing = T),]

colnames(asymmetry)=c("patient","mut","VAF","low","high")

plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(-0.4,1.2),ylim=c(-2,35))
for(n in 1:nrow(asymmetry)){
  rect(xleft=0,xright=min(2*as.numeric(asymmetry$VAF[n]),1),ytop=35-n+0.45,ybottom=35-n-0.45,col = rgb(18,90,168,maxColorValue = 255))
  rect(xright=1,xleft=min(2*as.numeric(asymmetry$VAF[n]),1),ytop=35-n+0.45,ybottom=35-n-0.45,col = rgb(201,201,201,maxColorValue = 255))
  text(label=asymmetry$patient[n],x=0,y=35-n,pos=2)
}

#Figures 3d-f
options(stringsAsFactors = F)
library(ape)
library(ggtree)
library(data.table)

patient="PD42142"
UC="PD42141b"
opt="snp"

Muts_on_branches=read.table(paste0(patient,"/",patient,"_snp_assigned_to_branches.txt"),sep="\t",header=T)
rownames(Muts_on_branches)=paste(Muts_on_branches$Chr,Muts_on_branches$Pos,Muts_on_branches$Ref,Muts_on_branches$Alt,sep="_")
tree=read.tree(paste0(patient,"/",patient,"_snp_tree_with_branch_length.tree"))
NV=read.table(paste0(patient,"/snp_NV_filtered_all.txt"))
NR=read.table(paste0(patient,"/snp_NR_filtered_all.txt"))
Genotype=NV/NR

tree_collapsed=tree
tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
tree_collapsed$edge.length[tree$edge.length==0]=0
tree_df=as.data.frame(fortify(tree_collapsed))

root=which(tree_df$node==tree_df$parent)
daughters=which(tree_df$parent==root&tree_df$node!=root)

if(length(daughters)==2){
  cont1=mean(Genotype[rownames(Muts_on_branches)[Muts_on_branches$Branch==daughters[1]],UC])
  cont2=mean(Genotype[rownames(Muts_on_branches)[Muts_on_branches$Branch==daughters[2]],UC])
  p=ggtree(tree_collapsed)
  dat1 = data.frame(node=daughters,
                    V1=c(cont1,cont2),
                    V2=c(cont2,cont1))
}
if(length(daughters)==3){
  cont1=mean(Genotype[rownames(Muts_on_branches)[Muts_on_branches$Branch==daughters[1]],UC])
  cont2=mean(Genotype[rownames(Muts_on_branches)[Muts_on_branches$Branch==daughters[2]],UC])
  cont3=mean(Genotype[rownames(Muts_on_branches)[Muts_on_branches$Branch==daughters[3]],UC])
  tot=cont1+cont2+cont3
  
  p=ggtree(tree_collapsed)
  dat1 = data.frame(node=daughters,
                    V1=c(cont1,cont2,cont3),
                    V2=tot-c(cont1,cont2,cont3))
}

pies = nodepie(dat1, cols=2:3,color=c('firebrick',"lightgrey"))

pdf(paste0(patient,"_pies.pdf"),width=5,height=3)
inset(p,pies,width=0.25,height=0.25)
dev.off()

#Figure 3g

options(stringsAsFactors = F)
library(data.table)

SNPs=fread("../chr10_1000G_loci.bed",header=F,data.table = F)
genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"
library("GenomicRanges")
library("Rsamtools")
library("MASS")
Ref=as.vector(scanFa(genomeFile, GRanges("10", IRanges(as.numeric(SNPs$V2),as.numeric(SNPs$V2)))))
samples=c("PD45581c","PD45581e","PD45581f","PD45665b")
coords=paste(SNPs$V1,SNPs$V2,sep="_")
rownames(SNPs)=coords
all_counts = array(0,dim=c(length(samples),length(coords),4),
                   dimnames=list(samples,coords,c("A","C","G","T")))


Nref=Ntot=matrix(0,ncol=length(samples),nrow=length(coords))
for (k in 1:length(samples)){
  data=fread(paste0(samples[k],".allelecounts.chr10.txt"),header=T,data.table=F)
  Ntot[,k]=data$Good_depth
  for(n in c("A","C","G","T")){
    Nref[Ref==n,k]=data[Ref==n,paste0("Count_",n)]
  }
  print(k)
}
rownames(Nref)=rownames(Ntot)=coords
colnames(Nref)=colnames(Ntot)=samples

BAF=1-Nref/Ntot
BAF=BAF[rowSums(is.infinite(BAF))==0,]
BAF=BAF[rowSums(is.nan(BAF))==0,]
BAF=BAF[rowMeans(BAF)>0,]

col=rep("grey40",length(BAF))
col[BAF[,"PD45665b"]<0.1&BAF[,"PD45581c"]>0.1]="dodgerblue" #select paternal SNPs

pos=SNPs[rownames(BAF),"V2"]
x_max=pos[length(pos)]/10^6

plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(-0.15*x_max,1.05*x_max),ylim=c(-0.5,3.5))
points(x=pos/10^6,y=2.3+BAF[,"PD45581c"],col=col,pch=16,cex=0.14)
points(x=pos/10^6,y=1.15+BAF[,"PD45581e"],col=col,pch=16,cex=0.14)
points(x=pos/10^6,y=BAF[,"PD45581f"],col=col,pch=16,cex=0.14)
segments(x0=-0.01*x_max,y0=-0.05,y1=3.3)
segments(x0=-0.01*x_max,y0=-0.05,x1=x_max)
segments(x0=seq(0,125,by=25),y0=-0.05,y1=-0.1)
segments(x0=-0.01*x_max,x1=-0.02*x_max,y0=c(c(0,0.5,1),1.15+c(0,0.5,1),2.3+c(0,0.5,1)))
text(y=c(c(0,0.5,1),1.15+c(0,0.5,1),2.3+c(0,0.5,1)),x=-0.02*x_max,labels=rep(c(0,0.5,1),3),pos=2)
text(x=seq(0,125,by=25),label=seq(0,125,by=25),y=-0.2)
text(x=x_max/2,y=-0.4,label="Position on Chr10 (Mb)")
text(x=-0.12*x_max,y=1.65,label="BAF",srt=90)
text(x=1.025*x_max,y=c(0.5,1.65,2.8),label=paste0("PD45581",c("f","e","c")),srt=270,font=2)
