
library(ggbeeswarm)
library(ggplot2)

manifest=read.csv("Extended_data_table2.csv") #Read in Extended data table 2

#Figure 2b

lcm_samples = manifest[manifest$Histo_desc != 'Bulk' & manifest$Patient != 'PD42146',]
lcm_samples$Histo_desc = factor(lcm_samples$Histo_desc, levels = c('Trophoblast', 'Mesenchymal core'))

lcm_samples[lcm_samples$Histo_desc == 'Mesenchymal core', ]$SNVs_medianVAF

#is the difference in VAF between the groups statistically significant?
wilcox.test(lcm_samples[lcm_samples$Histo_desc == 'Trophoblast', ]$SNVs_medianVAF, lcm_samples[lcm_samples$Histo_desc == 'Mesenchymal core', ]$SNVs_medianVAF, alternative = "two.sided")
#W = 1929, p-value = 9.941e-13

median(lcm_samples[lcm_samples$Histo_desc == 'Trophoblast', ]$SNVs_medianVAF) #0.387993
median(lcm_samples[lcm_samples$Histo_desc == 'Mesenchymal core', ]$SNVs_medianVAF) #0.2020832

pdf('/lustre/scratch119/casm/team274sb/to3/placenta/final_files/tb_v_mc_placenta_VAF_20200702.pdf', height = 5, width = 4, useDingbats = F)
ggplot(data = lcm_samples) + 
  geom_beeswarm(mapping = aes(x = Histo_desc, y = SNVs_medianVAF, col = Histo_desc), cex=4.5, shape=19, dodge.width = 0.8, size = 3) +
  geom_pointrange(mapping = aes(x = Histo_desc, y = SNVs_medianVAF),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black") + coord_cartesian(clip = 'off', ylim = c(0, 0.5), expand = T) +
  theme_classic() +
  theme(panel.background = element_rect(colour=NA, fill=NA),
        axis.text.x = element_text(size = 11, vjust = -2), panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 11), axis.text.y = element_text(size = 11), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), face = 'bold')) +
  labs(x = "Histological description", y = 'Median substitution VAF') + scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), expand = expansion(add = c(0.05, 0))) + scale_color_manual(values = c("#E6007EFF", "#999999"))
dev.off()

#Figure 2c

options(stringsAsFactors = F)
library(ggplot2)
library(ggtree)
library(ape)
library(data.table)
patients=read.table("patients.txt")[,1]

find_parent=function(sample,tree_df){
  child=which(tree_df$label==sample)
  parent=tree_df$parent[child]
  all_nodes=c(child,parent)
  while(child!=parent){
    child=parent
    parent=tree_df$parent[child]
    all_nodes=c(all_nodes,parent)
  }
  return(unique(all_nodes))
}

MCs=manifest$Sample[manifest$Histo_desc=="Mesenchymal core"]
biopsy_TB_all=biopsy_MC_all=c()
TB_dists_all=list()
MC_dists_all=list()
for(patient in patients){
  tree=read.tree(paste0(patient,"_snp_tree_with_branch_length.tree"))
  tree_df=as.data.frame(fortify(tree))
  samples=tree_df$label[tree_df$isTip]
  samples=samples[grepl(patient,samples)]
  MC_samples=samples[samples%in%MCs]
  TB_samples=samples[!samples%in%MCs]
  if(grepl("PD421",patient)){
    biopsies=paste0(patient,"b",c("",2:4))
  }else{
    biopsies=paste0(patient,letters[2:5])
  }
  TB_dists=list()
  MC_dists=list()
  for(b in biopsies){
    biopsy_TB=TB_samples[unlist(strsplit(TB_samples,split = "_"))[seq(1,(2*length(TB_samples)-1),by=2)]==b]
    biopsy_MC=MC_samples[unlist(strsplit(MC_samples,split = "_"))[seq(1,(2*length(MC_samples)-1),by=2)]==b]
    biopsy_MC=biopsy_MC[!biopsy_MC%in%MC_remove]
    biopsy_MC_all=c(biopsy_MC_all,biopsy_MC)
    biopsy_TB_all=c(biopsy_TB_all,biopsy_TB)
    TB_dists[[b]]=matrix(NA,nrow=length(biopsy_TB),ncol=length(biopsy_TB))
    for(i in 1:(length(biopsy_TB)-1)){
      for(j in (i+1):length(biopsy_TB)){
        nodes=c(find_parent(biopsy_TB[i],tree_df),find_parent(biopsy_TB[j],tree_df))
        ancestral_nodes=names(table(nodes))[table(nodes)==2]
        if(length(ancestral_nodes)){
          ancestral_length=sum(tree_df$branch.length[as.numeric(ancestral_nodes)])
          TB_dists[[b]][i,j]=ancestral_length/mean(tree_df$x[tree_df$label%in%c(biopsy_TB[i],biopsy_TB[j])])
        }else{
          TB_dists[[b]][i,j]=0
        }
      }
    }
    if(length(biopsy_MC)>1){
      MC_dists[[b]]=matrix(NA,nrow=length(biopsy_MC),ncol=length(biopsy_MC))
      for(i in 1:(length(biopsy_MC)-1)){
        for(j in (i+1):length(biopsy_MC)){
          nodes=c(find_parent(biopsy_MC[i],tree_df),find_parent(biopsy_MC[j],tree_df))
          ancestral_nodes=names(table(nodes))[table(nodes)==2]
          if(length(ancestral_nodes)){
            ancestral_length=sum(tree_df$branch.length[as.numeric(ancestral_nodes)])
            MC_dists[[b]][i,j]=ancestral_length/mean(tree_df$x[tree_df$label%in%c(biopsy_MC[i],biopsy_MC[j])])
          }else{
            MC_dists[[b]][i,j]=0
          }
        }
      }
    }
  }
  TB_dists_all[[patient]]=TB_dists
  MC_dists_all[[patient]]=MC_dists
} 

TB_dist_vec=unlist(TB_dists_all)
TB_dist_vec=as.numeric(TB_dist_vec[!is.na(TB_dist_vec)])
MC_dist_vec=unlist(MC_dists_all)
MC_dist_vec=as.numeric(MC_dist_vec[!is.na(MC_dist_vec)])

wilcox.test(TB_dist_vec,MC_dist_vec)
mean(MC_dist_vec)
mean(TB_dist_vec)

#Figure 2e
plot(0,type = "n", xaxt = "n", yaxt = "n", xlab = "", bty="n",
     ylab = "",xlim=c(-1,9),ylim=c(-0.2,1.2))
a=1
for(p in patients){
  points(x=jitter(x=rep(a-0.5,sum(grepl(p, names(TB_dist_vec)))),amount=0.3),y=as.numeric(TB_dist_vec[grepl(p, names(TB_dist_vec))]),pch=16,col="grey60")
  segments(x0=a-0.75,x1=a-0.25,y0=mean(as.numeric(TB_dist_vec[grepl(p, names(TB_dist_vec))])),lwd=5)
  a=a+1
}
points(x=jitter(x=rep(5.5,length(MC_dist_vec)),amount=0.3),y=MC_dist_vec,pch=16,col="grey60")
segments(x0=5.25,x1=5.75,y0=mean(MC_dist_vec),lwd=5))

segments(x0=0,y0=-0.015,y1=1,lwd=1)
segments(x0=0,x1=-0.1,y0=seq(0,1,0.2))
text(x=-0.1,y=seq(0,1,0.2),labels=seq(0,1,0.2),pos=2)
text(x=c(2.5,5.5),y=-0.1,labels=c("Trophoblasts","MC"))

segments(x0=0,x1=6,y0=-0.015,lwd=1)
segments(x0=5,y0=-0.015,y1=1,lwd=2,lty='dashed')
dev.off()



