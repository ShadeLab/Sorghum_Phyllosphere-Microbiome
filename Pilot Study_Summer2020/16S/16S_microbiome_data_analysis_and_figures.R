library(RColorBrewer)
library(viridis)
library(microbiome)
library(reshape2)
library(ggpubr)
library(broom)
library(ggfortify)
library(dplyr)
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(indicspecies)
library(lemon)
library(ggalluvial)
library(pheatmap)
library(metagMisc)
library(metagenomeSeq)
library(qiime2R)
library(DESeq2)

install.packages("processx")
library(devtools)

devtools::install_github("vmikk/metagMisc")
library(metagMisc)

#Set working directory
setwd("/Volumes/ShadeLab/WorkingSpace/MarcoMechan_WorkingSpace/Sorghum_PilotStudy")

#Import files
otu=read.table("phyloseq/otu_table.txt", header = T, row.names = 1)
tax=read.delim("phyloseq/taxonomy.tsv",row.names = 1)
map=read.table('metadata.txt', header = T, row.names = 1)

tax=tax[-2]
tax_df <- colsplit(tax$Taxon, '; ', names =  c("Kingdom", "Phylum", "Class", 
                                               "Order", "Family", "Genus", "Species"))
tax_df[1:7] <- lapply(tax_df[1:7], function(x) gsub(".*__", "", x))
rownames(tax_df) <- rownames(tax)

OTU=otu_table(as.matrix(otu), taxa_are_rows = T)
TAX=tax_table(as.matrix(tax_df))
MAP=sample_data(map)

otuPhyloseq=phyloseq(OTU,TAX,MAP)
sample_sums(otuPhyloseq)

#Filtering mitochondria, Chloroplast and Unclassified taxa
otuPhyloseq_filt <- otuPhyloseq %>%
  subset_taxa(Kingdom != "Unassigned")

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Family != "Mitochondria") | is.na(Family))

otuPhyloseq_filt <- otuPhyloseq_filt %>%
  subset_taxa((Order != "Chloroplast") | is.na(Class))

sample_sums(otuPhyloseq_filt) 

# Subset Mucilage samples
Muc <- subset_samples(otuPhyloseq_filt, Compartment%in%c("Mucilage"))
Muc1 <- prune_taxa(taxa_sums(Muc) > 0, Muc)
summary(sample_data(Muc1)$Compartment)
print(Muc1)

Muc1.rarefied = rarefy_even_depth(Muc1, rngseed=1, sample.size=0.9*min(sample_sums(Muc1)), replace=F)
sample_sums(Muc1)
sample_sums(Muc1.rarefied)

# Relative abundance
Muc_phylum <- Muc1.rarefied %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  #filter(Abundance > 0.01) %>%
  arrange(Phylum)

write.csv(Muc_phylum, 'Mucilage_RA_phylum.csv') 

n <- dim(Muc_phylum)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(Muc_phylum,aes(x=Sample,y=Abundance,fill=Phylum)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance \n") +
  xlab("Sorghum mucilage") +
  theme(axis.text.x=element_text(size=16,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Phylum Composition")


Muc_family <- Muc1.rarefied %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Family)

write.csv(Muc_family, 'Mucilage_RA_family.csv')

n <- dim(Muc_family)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(Muc_family,aes(x=Sample,y=Abundance,fill=Family)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Family > 1%) \n") +
  xlab("Sorghum mucilage") +
  theme(axis.text.x=element_text(size=16,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Family Composition")


Muc_genus <- Muc1.rarefied %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(Genus)

write.csv(Muc_genus, 'Mucilage_RA_genus.csv')

n <- dim(Muc_genus)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(Muc_genus,aes(x=Sample,y=Abundance,fill=Genus)) +
  geom_bar(position="fill",stat="identity") + 
  scale_fill_manual(values = col_vector) + 
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) + 
  ylab("Relative Abundance (Genus > 1%) \n") +
  xlab("Sorghum mucilage") +
  theme(axis.text.x=element_text(size=16,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=18, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=18,face="bold", vjust = 10),
        legend.text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("Genus Composition")

rarecurve(t(otu_table(Muc1)), step=50, cex=0.5)

#alpha Diversity
alpha_meas = c("Observed", "Chao1", "Shannon", "Simpson")
(p <- plot_richness(Muc1.rarefied, "Site", measures=alpha_meas))
p + geom_boxplot(data=p$data, aes(x=Site, y=value, color=NULL), alpha=0.1)

rich = estimate_richness(Muc1.rarefied)
pairwise.wilcox.test(rich$Observed, sample_data(Muc1.rarefied)$Site)
pairwise.wilcox.test(rich$Shannon, sample_data(Muc1.rarefied)$Site)
pairwise.wilcox.test(rich$Simpson, sample_data(Muc1.rarefied)$Site)

# Bray-Curtis
bray_diss_Muc = phyloseq::distance(Muc1.rarefied, method="bray")
ordination_Muc = ordinate(Muc1.rarefied, method="PCoA", distance=bray_diss_Muc)

adonOut<- adonis(bray_diss_Muc ~ sample_data(Muc1.rarefied)$Site)
pVal=adonOut$aov.tab[6]$`Pr(>F)`[1]

plot_ordination(Muc1.rarefied, ordination_Muc, color = 'Site') + theme(aspect.ratio=1)+
  geom_point(size=3) + 
  theme(axis.text.x=element_text(size=14,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2)),
        axis.title=element_text(size=14,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  annotate(geom = 'text', label=paste("PERMANOVA\nP-val=",pVal, sep=''), x=.3, y=-.1)
  ggtitle("PCoA: Bray-Curtis")
  
  # Differential abundance of OTUs
  deseq_Muc1 = phyloseq_to_deseq2(Muc1, ~ Site)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(deseq_Muc1), 1, gm_mean)
deseq_Muc1 = estimateSizeFactors(deseq_Muc1, geoMeans = geoMeans)
deseq_Muc1 = DESeq(deseq_Muc1, fitType="local")

res = results(deseq_Muc1)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
res = results(deseq_Muc1, contrast=c("Site", "Fertilized", "Non_Fertilized"), alpha=alpha)
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Muc1)[rownames(sigtab), ], "matrix"))
head(sigtab)
# To write all OTUs that were significant different: positives and negatives
sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(sigtab, 'DEseq_all_values_Muc1.csv')

# To subset positives values
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

write.csv(posigtab, 'Differential_abundance_Muc1.csv')

theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x=element_text(size=16,angle=90,hjust =0.5),
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size=16,colour = "black", face = "bold"), 
        strip.text.y = element_text(size=16, face = 'bold'),
        plot.title = element_text(size = rel(2), face="bold"),
        axis.title=element_text(size=16,face="bold", vjust = 10),
        legend.text = element_text(size=14)) +
  ylab("TAXA \n") +
  xlab("Log2FoldChange") +
  ggtitle("OTU DESeq Mucilage: Fertilized vs Non-Fertilized")
