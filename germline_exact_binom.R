#-------------------------------------------------
# Filtering germline substitutions using an 
# exact binomial test
# Tim Coorens - January 2019
#-------------------------------------------------

options(stringsAsFactors = F)

#-------------------------------------------------
# Functions
#-------------------------------------------------

#1. Exact binomial test

exact.binomial=function(gender,NV,NR,cutoff=-5){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting 
  # variants and NR the matrix with total depth (samples as columns, mutations rows, 
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  XY_chromosomal = grepl("X|Y",rownames(NR))
  autosomal = !XY_chromosomal
  
  if(gender=="female"){
    NV_vec = rowSums(NV)
    NR_vec = rowSums(NR)
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    NV_vec = rowSums(NV)[autosomal]
    NR_vec = rowSums(NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    NV_vec = rowSums(NV)[XY_chromosomal]
    NR_vec = rowSums(NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
      if (n%%1000==0){
        print(n)
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  qval = p.adjust(pval,method="BH")
  germline = log10(qval)>cutoff
  return(germline)
}


# Toy data set, coverage of 30X, supposed "clonal" biopsies

NV = rbind(rbinom(p=0.5,n=10,size=30),
           rbinom(p=0.5,n=10,size=30),
           rbinom(p=0.5,n=10,size=30),
           rbinom(p=0.95,n=10,size=30),
           c(rep(0,5),rbinom(p=0.5,n=5,size=30)),
           c(rep(0,3),rbinom(p=0.5,n=7,size=30)),
           c(rep(0,4),rbinom(p=0.95,n=6,size=30)),
           c(rbinom(p=0.5,n=3,size=30),rep(0,7)))
NR=matrix(30,ncol=10,nrow=nrow(NV))
rownames(NV)=rownames(NR)=c(paste0(c(1:3,"X"),"_germ"),
                            paste0(c(1,2,"X",4),"_somatic"))
colnames(NV)=colnames(NR)=paste0("sample",1:10)
germline=exact.binomial(NV=NV,NR=NR,gender="male")
rownames(NV)[germline] #Should be all variants labeled as 'germ'
rownames(NV)[!germline] #Should be all variants labeled as 'somatic'

# Toy data set, coverage of 30X, bulk samples

NV = rbind(rbinom(p=0.5,n=10,size=30),
           rbinom(p=0.5,n=10,size=30),
           rbinom(p=0.5,n=10,size=30),
           rbinom(p=0.95,n=10,size=30),
           rbinom(p=0.15,n=10,size=30),
           rbinom(p=0.25,n=10,size=30),
           rbinom(p=0.5,n=10,size=30),
           rbinom(p=0.1,n=10,size=30))
NR=matrix(30,ncol=10,nrow=nrow(NV))
rownames(NV)=rownames(NR)=c(paste0(c(1:3,"X"),"_germ"),
                            paste0(c(1,2,"X",4),"_somatic"))
colnames(NV)=colnames(NR)=paste0("sample",1:10)
germline=exact.binomial(NV=NV,NR=NR,gender="male")
rownames(NV)[germline] #Should be all variants labeled as 'germ'
rownames(NV)[!germline] #Should be all variants labeled as 'somatic'