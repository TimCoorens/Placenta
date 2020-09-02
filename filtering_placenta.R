options(stringsAsFactors = F)
patient=commandArgs(T)[1]
opt=commandArgs(T)[2] #either 'snp' or 'indel'

library(data.table)

source("/lustre/scratch117/casm/team274/tc16/Scripts/R_scripts/germline_exact_binom.R")
source("/lustre/scratch117/casm/team274/tc16/Scripts/R_scripts/beta_binom_flt.R")

samples_with_CNVs=read.table("samples_with_CNVs.txt")[,1]
data = fread(paste0(patient,"/",patient,".",opt,".tsv"),header=T,data.table=F)
Muts = paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_")
Genotype = data[,grepl("VAF",colnames(data))&colnames(data)!="PDv37is_VAF"]
NR = data[,grepl("DEP",colnames(data))&colnames(data)!="PDv37is_DEP"]
NV = data[,grepl("MTR",colnames(data))&colnames(data)!="PDv37is_MTR"]

rownames(Genotype)=rownames(NV)=rownames(NR)=Muts
samples=colnames(Genotype)=colnames(NR)=colnames(NV)=gsub("_VAF","",colnames(Genotype))
samples_patient=grepl(patient,samples)

XY_chromosomal = grepl("X|Y",Muts)
autosomal = !XY_chromosomal
xy_depth=mean(rowMeans(NR[XY_chromosomal,samples_patient]))
autosomal_depth=mean(rowMeans(NR[autosomal,samples_patient]))

gender='male'
if(xy_depth>0.8*autosomal_depth) gender='female'

noCNVs=!samples%in%samples_with_CNVs
select_samples=samples_patient&noCNVs

#filter out sites of consistently low and high depth
low=10
upper=100
if(gender=="male") depth_filter=(rowMeans(NR)>low&rowMeans(NR)<upper&autosomal)|
  (rowMeans(NR[,noCNVs])>low&rowMeans(NR[, noCNVs])<0.5*upperXY_chromosomal)
if(gender=="female") depth_filter=rowMeans(NR[,noCNVs])>low&rowMeans(NR[,noCNVs])<upper

germline=exact.binomial(gender=gender,NV=NV[,select_samples],NR=NR[,select_samples],cutoff = -5) #determine which variants are germline

write.table(Muts[germline],paste0(patient,"/germline_",opt,"_ids.txt"),row.names = F,col.names = F,quote=F)
write.table(Muts[!germline],paste0(patient,"/somatic_",opt,"_ids.txt"),row.names = F,col.names = F,quote=F)

NR_flt=NR[!germline&depth_filter,]
NV_flt=NV[!germline&depth_filter,]

NR_flt_nonzero=NR_flt
NR_flt_nonzero[NR_flt_nonzero==0]=1 # avoid creating NAs
shared_muts=rowSums(NV_flt>0)>1

rho_est = beta.binom.filter(NR=NR_flt_nonzero[shared_muts,select_samples],NV=NV_flt[shared_muts,select_samples])
flt_rho=log10(rho_est)<(-1)
rho_filtered_out = rownames(NR_flt_nonzero[shared_muts,])[flt_rho]
write.table(rho_filtered_out,paste0(patient,"/",opt,"_bbinom_filtered_out.txt"))
write.table(rho_est,paste0(patient,"/",opt,"_rho_est.txt"))

NR_flt_2 = NR_flt[!rownames(NR_flt)%in%rho_filtered_out,]
NV_flt_2 = NV_flt[!rownames(NV_flt)%in%rho_filtered_out,]

write.table(NR_flt_2,paste0(patient,"/",opt,"_NR_filtered_all.txt"))
write.table(NV_flt_2,paste0(patient,"/",opt,"_NV_filtered_all.txt"))

#If bulk sample, stop here
#If samples are used as input for a phylogeny builder, continue:

NR_flt_nonzero=NR_flt_2
NR_flt_nonzero[NR_flt_nonzero==0]=1

XY_chromosomal=grepl("X|Y",rownames(NR_flt_2))
autosomal=!XY_chromosomal
genotype_bin=as.matrix(NV_flt_2/NR_flt_nonzero)
if(gender=="male"){
  genotype_bin[autosomal,][genotype_bin[autosomal,]<0.1]=0
  genotype_bin[autosomal,][genotype_bin[autosomal,]>=0.3]=1
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]<0.2]=0
  genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]>=0.6]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
}
if(gender=="female"){
  genotype_bin[genotype_bin<0.1]=0
  genotype_bin[genotype_bin>=0.3]=1
  genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
}

muts=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(genotype_bin),split="_")),byrow = T))

#If LOH events, change presence to '?'
if (any(!noCNVs)){
  samples_cnv=samples[!noCNVs]
  for (s in samples_cnv){
    cnvs=read.csv(paste0("/nfs/cancer_ref01/nst_links/live/2058/",s,"/",s,".ascat_ngs.summary.csv"),header=F) #point to ASCAT or battenberg output
    loh=cnvs[cnvs$V8==0,]
    muts_affected=apply(muts,1,function(x) any(x[1]==loh$V1&x[2]>=loh$V2&x[2]<=loh$V3))
    print(paste0(sum(muts_affected&genotype_bin[,s]==0)," mutations in LOH regions in sample ",s))
    genotype_bin[muts_affected&genotype_bin[,s]==0,s]=0.5
  }
}

genotype_bin=genotype_bin[,select_samples&lcm_samples]
dna_strings = list()

if(opt=="snp"){
  Ref = muts[,3]
  Alt = muts[,4]
  dna_strings[1]=paste(Ref,sep="",collapse="")
  for (k in 1:ncol(genotype_bin)){
    Mutations = Ref
    Mutations[genotype_bin[,k]==0.5] = '?'
    Mutations[genotype_bin[,k]==1] = Alt[genotype_bin[,k]==1]
    dna_string = paste(Mutations,sep="",collapse="")
    dna_strings[k+1]=dna_string
  }
}
if(opt=="indel"){
  Ref = rep("A",nrow(genotype_bin))
  Alt = rep("T",nrow(genotype_bin))
  dna_strings[1]=paste(Ref,sep="",collapse="")
  for (n in 1:ncol(genotype_bin)){
    Mutations = Ref
    Mutations[genotype_bin[,n]==0.5] = '?'
    Mutations[genotype_bin[,n]==1] = Alt[genotype_bin[,n]==1]
    dna_string = paste(Mutations,sep="",collapse="")
    dna_strings[n+1]=dna_string
  }
}
bed_muts=cbind(muts[,1],as.numeric(muts[,2])-1,muts[,2])
bed_muts=bed_muts[order(bed_muts[,1],as.numeric(bed_muts[,2])),]
write.table(bed_muts,paste0(patient,"/",opt,"_good.bed"),quote=F,row.names=F,col.names=F,sep="\t")

names(dna_strings)=c("Ancestral",colnames(genotype_bin))
require(seqinr)
write.fasta(dna_strings, names=names(dna_strings),paste0(patient,"/",opt,"_for_MPBoot.fa"))
