library(ggplot2)
library(ggbeeswarm)
library(reshape2)

manifest=read.csv("Extended_data_table2.csv") #Read in Extended data table 2

#Figure 1b
#to3

bulk_mb = manifest[manifest$Histo_desc == 'Bulk' & manifest$Patient != 'PD45581',] #remove trisomic rescue sample
bulk_mb$status = ifelse(bulk_mb$Clinical_group == 8, 'Healthy', 'Abnormal') #divide by whether or not the placenta was classified as healthy
bulk_mb$status = factor(bulk_mb$status, levels = c('Healthy', 'Abnormal'))

mean(bulk_mb$adjusted_SNV_MB_for_sensitivity) #144.9859
min(bulk_mb$adjusted_SNV_MB_for_sensitivity) #38.11963
max(bulk_mb$adjusted_SNV_MB_for_sensitivity) #259.4976

ggplot(data = bulk_mb) + 
  geom_beeswarm(mapping = aes(x = status, y = adjusted_SNV_MB_for_sensitivity), cex=4.5, shape=19, dodge.width = 0.8, size = 3, col = "#999999") +
  geom_pointrange(mapping = aes(x = status, y = adjusted_SNV_MB_for_sensitivity),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black") + coord_cartesian(clip = 'off', ylim = c(0, 300), expand = T) +
  theme_classic() +
  theme(panel.background = element_rect(colour=NA, fill=NA),
        axis.text.x = element_text(size = 11, vjust = -2), panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.text.y = element_text(size = 11), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), face = 'bold')) +
  labs(x = "Bulk samples", y = 'Substitution burden (adjusted)') + scale_y_continuous(breaks = c(0, 100, 200, 300), labels = c(0, 100, 200, 300), expand = expansion(add = c(5, 0))) 

#Figure 1c
#to3

bulk_mb = manifest[manifest$Histo_desc == 'Bulk' & manifest$Patient != 'PD45581',] #remove trisomic rescue sample
bulk_mb$status = ifelse(bulk_mb$Clinical_group == 8, 'Healthy', 'Abnormal')
bulk_mb$status = factor(bulk_mb$status, levels = c('Healthy', 'Abnormal'))

mean(bulk_mb$SNVs_medianVAF) #0.2437741
min(bulk_mb$SNVs_medianVAF) #0.1463415
max(bulk_mb$SNVs_medianVAF) #0.4418605

ggplot(data = bulk_mb) + 
  geom_beeswarm(mapping = aes(x = status, y = SNVs_medianVAF), cex=4.5, shape=19, dodge.width = 0.8, size = 3, col = "#999999") +
  geom_pointrange(mapping = aes(x = status, y = SNVs_medianVAF),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black") + coord_cartesian(clip = 'off', ylim = c(0, 0.5), expand = T) +
  theme_classic() +
  theme(panel.background = element_rect(colour=NA, fill=NA),
        axis.text.x = element_text(size = 11, vjust = -2), panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), axis.text.y = element_text(size = 11), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), face = 'bold')) +
  labs(x = "Bulk samples", y = 'Median substitution VAF') + scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), expand = expansion(add = c(0.05, 0))) 

#Figure 1d
#to3

library(reshape2)
library(ggplot2)

bulk_sigs <- manifest[manifest$Histo_desc == 'Bulk', c('Sample', 'Patient', 'Clinical_group', 'prop_SBS1', 'prop_SBS5', 'prop_SBS18')]
row.names(bulk_sigs) = bulk_sigs$Sample
bulk_sigs_norm <- bulk_sigs[bulk_sigs$Patient != 'PD45581',]

bulk_sigs_norm$Patient <- as.factor(bulk_sigs_norm$Patient)
bulk_sigs_norm <- bulk_sigs_norm[order(bulk_sigs_norm$Patient),]
bulk_sigs_norm$not_fit <- 1 - (bulk_sigs_norm$prop_SBS1 + bulk_sigs_norm$prop_SBS5 + bulk_sigs_norm$prop_SBS18) 
bulk_sigs_norm[bulk_sigs_norm$not_fit < 0,]$not_fit = 0 #rounding errors create tiny minus values
bulk_sigs_norm$status = ifelse(bulk_sigs_norm$Clinical_group == 8, 'Healthy', 'Abnormal')
bulk_sigs_norm$status = factor(bulk_sigs_norm$status , levels = c('Healthy', 'Abnormal'))

mean(bulk_sigs_norm$prop_SBS18) #0.434562

dfm.bulk <- reshape2::melt(bulk_sigs_norm[, c(1,4:8)], id.vars = c("Sample", "status"))
dfm.bulk$variable <- factor(dfm.bulk$variable, levels = c('not_fit', 'prop_SBS18', 'prop_SBS5', 'prop_SBS1'))

ggplot(dfm.bulk, aes(x = Sample, y = value)) + 
  geom_bar(aes(fill = variable), stat = "identity", position = "fill", width = 1) + 
  facet_grid(.~ status, scales = 'free_x', space = 'free_x', switch = 'x') +
  theme(panel.spacing = unit(0, "mm"),
        strip.background = element_blank(),
        strip.placement = 'outside',
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(colour = 'black'), 
        strip.text.x = element_text(size = 11), 
        axis.text = element_text(size = 11), 
        axis.title.x = element_text(face="bold", vjust = -3), 
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        plot.title = element_text(size=16, face="bold"), 
        panel.background = element_rect(fill = NA, colour = 'black', size = 0.5),
        panel.ontop = T,
        legend.position = 'right') +
  theme(plot.margin=unit(c(0.5, 0.1, 0.5, 0.5),"cm")) +
  coord_cartesian(ylim = c(0, 1), expand = F) + scale_x_discrete(expand = c(0,0)) +
  labs(y = 'SBS signature contribution', x = 'Bulk samples') +
  scale_fill_manual(values = c("#656565", "white", "#009FE3FF", "black"), name = '', limits = rev(c('not_fit', 'prop_SBS18', 'prop_SBS5', 'prop_SBS1')), labels = c("SBS1", "SBS5", "SBS18",  "Not fit")) +
  guides(fill = guide_legend(override.aes = list(colour = "black")))

#Figure 1e
#to3

#Data from Lee-Six et al, Nature, 2019
intestine_wgs <- read.table('/lustre/scratch119/casm/team274sb/to3/placenta/reference_datasets/lee_six_intestine_sigs.txt', header = T, sep = '\t', stringsAsFactors = F)

intestine_wgs$intestine <- 'LargeIntestine' #split these anatomically distinct structures with distinct malignancy profiles
for(i in 1:nrow(intestine_wgs)){
  if(intestine_wgs$site[i] == 'Ileum') intestine_wgs$intestine[i] = 'Ileum' 
}

bowel_analysis = intestine_wgs[, c('crypt', 'intestine', 'sbs.SBS18', 'sbstotal')]
bowel_analysis$sbs18.prop = bowel_analysis$sbs.SBS18 / bowel_analysis$sbstotal

bowel_plot = bowel_analysis[, c('crypt', 'intestine', 'sbs18.prop')]
names(bowel_plot) = c('Sample', 'status', 'prop_SBS18')

#read in our bulk and microdissection samples data, exclude the trisomic rescue
bulk_mb = manifest[manifest$Histo_desc == 'Bulk' & manifest$Patient != 'PD45581',] #remove trisomic rescue sample
bulk_mb$status = ifelse(bulk_mb$Clinical_group == 8, 'Healthy', 'Abnormal')
bulk_mb$status = factor(bulk_mb$status, levels = c('Healthy', 'Abnormal'))
placenta = bulk_mb[bulk_mb$Histo_desc %in% c('Bulk', 'Trophoblast'), c('Sample', 'status', 'prop_SBS18')]

mean(placenta$prop_SBS18) #43% when combined across all clinical groups

#combine datasets and order them
sbs18_plotting_data = rbind(bowel_plot, placenta)
sbs18_plotting_data$status = factor(sbs18_plotting_data$status, levels = c('Healthy', 'Abnormal', 'LargeIntestine', 'Ileum'))

aggregate(prop_SBS18 ~ status, sbs18_plotting_data, summary) #use mean values of prop_SBS18
#Healthy 45%, Abnormal 42%, Colorectal epithelium 13%, ileum 14%

wilcox.test(sbs18_plotting_data[sbs18_plotting_data$status == 'Healthy',]$prop_SBS18, sbs18_plotting_data[sbs18_plotting_data$status == 'LargeIntestine',]$prop_SBS18, alternative = "two.sided") #W = 16102, p-value < 2.2e-16
wilcox.test(sbs18_plotting_data[sbs18_plotting_data$status == 'Abnormal',]$prop_SBS18, sbs18_plotting_data[sbs18_plotting_data$status == 'LargeIntestine',]$prop_SBS18, alternative = "two.sided") #W = 19462, p-value < 2.2e-16

ggplot(data = sbs18_plotting_data) +
  geom_quasirandom(mapping = aes(x = status, y = prop_SBS18, col = status), groupOnX = T, pch = 19, size = 3) +
  scale_color_manual(values = c("#E6007EFF", "#E09DC3", "#999999", "#999999")) +
  scale_x_discrete(labels = c('Normal', 'Abnormal\nparameter', 'Colon', 'Ileum')) + 
  geom_pointrange(mapping = aes(x = status, y = prop_SBS18),
                  stat = "summary", shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, fill="black", size = 1) + coord_cartesian(clip = 'off', ylim = c(0, 1), expand = T) +
  theme_classic() +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(add = c(0, 0))) +
  labs(y = 'Proportion of single base substitutions attributed to SBS18') +
  theme(panel.background = element_rect(colour=NA, fill=NA),
        axis.text.x = element_text(size = 11), panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(colour = 'black'), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 11), axis.text.y = element_text(size = 11), axis.title.x = element_blank()) +
  geom_vline(xintercept = 2.5, linetype = 'dashed')
