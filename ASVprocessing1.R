setwd("~/AmpliconTest/")
library("dplyr")
library("phyloseq")
library("rstatix")
library("vegan")
library("ggplot2")
library("dada2")
library("ggpubr")
library("reshape2")
library(gplots)
library(picante)

#Load in R files from ASV picking/taxonomy assignment
load("intermediate/dada2_seqtab_taxa.rda")
#Load fresh metadata
meta1<-sample_data(metadata_df)
sample_names(meta1)<-meta1$ID
sample_names(meta1)
OTU<-otu_table(seqtab.nochim, taxa_are_rows=FALSE)
sample_names(OTU)
sample_names(taxa)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(meta1),
               tax_table(taxa))
ps
tree <- read_tree("Output/dada2_nochim.tree")
ps <- merge_phyloseq(ps, tree)

#Remove unclassified bacteria, and human/plantreads
ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))
ps.clean
#Filter low abundance taxa - less than 5 reads total
ps.clean.p0 <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 5}, prune=TRUE)
ps.clean.p0

#Collapse to phylum level
sample_meta <- sample_data(ps.clean.p0)
ps.clean.phylum = phyloseq::tax_glom(ps.clean.p0, taxrank = rank_names(ps.clean.p0)[2])
phyloseq::taxa_names(ps.clean.phylum) = phyloseq::tax_table(ps.clean.phylum)[,"Phylum"]
otu_table(ps.clean.phylum)[1:5,1:5]

#Collapse to genus level
ps.clean.re <- transform_sample_counts(ps.clean.p0, function(x) x / sum(x))
ps.clean.genus = phyloseq::tax_glom(ps.clean.re, taxrank = rank_names(ps.clean.re)[6])
phyloseq::taxa_names(ps.clean.genus) = phyloseq::tax_table(ps.clean.genus)[,"Genus"]
otu_table(ps.clean.genus)[1:5,1:5]

ps.clean.phylum1<-subset_samples(ps.clean.phylum,Treatment!="E4")

#Total read cound
phyloseq::psmelt(ps.clean.phylum1) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free", nrow = 2) +
  theme_classic()
ggsave("TreatmentEffectAtPhylum.pdf",width=11)

# alpha diversity should be calculated before filtering on abundance and prevalence
tree = phyloseq::phy_tree(ps)
samp = data.frame(phyloseq::otu_table(ps))

#Generate a data.frame with adiv measures
adiv <- data.frame(
  phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps)), ID,Day, Subject, Treatment)) %>%
  dplyr::select(-se.chao1)

measures = colnames(adiv)[1:6]
#head(adiv)

#dataframes to look remove or focus on controls
adiv2<-subset(adiv,Treatment!="E4")
adiv3<-subset(adiv,Day!="13D")
adiv3<-subset(adiv3,Day!="3D")
adiv3<-subset(adiv3,Day!="6D")
adiv3<-subset(adiv3,Day!="5D")
adiv3<-subset(adiv3,Day!="7D")

#adiv$Diet<-factor(adiv$Diet,levels=c("FF","W","P","AE"))
#adiv2$Diet<-factor(adiv2$Diet,levels=c("FF","W","P","AE"))
#Plot adiv measures
adiv3 %>%
  gather(key = metric, value = value, measures) %>% #c("Observed", "Shannon", "Chao1", "Simpson", "PD")
  mutate(metric = factor(metric, levels = measures)) %>% #c("Observed", "Shannon", "Chao1", "Simpson", "PD")
  ggplot(aes(x = Day, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Treatment), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measures") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none") +theme_classic()+scale_color_manual(values=c("brown","forestgreen","yellow4","turquoise"))+stat_compare_means()

#ggsave("DayEffectAlphaDiv.pdf",width=8.5,height=6.8)

#Plot single Metric
ggplot(data=adiv,aes(x=Day,y=InvSimpson))+facet_wrap(~ Treatment, scales = "free")+geom_boxplot()+geom_jitter(aes(color=Subject), height = 0, width = .2)+
  theme_classic()#+scale_color_manual(values=c("blue3","yellowgreen","violet"))+stat_compare_means()

#Focus on last two timepoints by subsetting based on DSS2
adiv4<-subset(adiv,DSS2!="NA")
adiv4<-subset(adiv4,ID!="101C")

adiv4$DSS2<-factor(adiv4$DSS2,levels=c("Pre","Post"))
ggplot(data=adiv4,aes(x=DSS2,y=InvSimpson))+facet_wrap(~ Diet, scales = "free")+geom_boxplot()+geom_point(aes(color=Microbiome),size=1.6)+
  theme_classic()+scale_color_manual(values=c("blue3","yellowgreen","violet"))+stat_compare_means()+geom_line(aes(group=Mouse,color=Microbiome),linewidth=.7,alpha=.7)

res = adiv %>%
  gather(key = metric, value = value, measures) %>% 
  mutate(metric = factor(metric, levels = measures)) %>%
  group_by(metric) %>%
  rstatix::wilcox_test( value ~ DSS2, p.adjust.method = "none") %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "group", scales = "free_y")

#Calculate beta diversity
ps.beta <- ps.clean.re
ordBC <- ordinate(ps.beta, "PCoA", "bray")
ordJC <- ordinate(ps.beta, "PCoA", "jaccard")
ordUF <- ordinate(ps.beta, "PCoA", "unifrac")
ordwUF <- ordinate(ps.beta, "PCoA", "wunifrac")
smpID <- sample_data(ps.beta)$ID

#Combine 1st 2 vectors of every metric 
df <- rbind(data.frame(ordBC$vectors[,1:2], sample = smpID, method = 'BC'),
            data.frame(ordJC$vectors[,1:2], sample = smpID,method = 'Jaccard'),
            data.frame(ordUF$vectors[,1:2], sample = smpID,method = 'unifrac'),
            data.frame(ordwUF$vectors[,1:2], sample = smpID,method = 'wunifrac'))

# add sample_data info
df <- merge(df, data.frame(sample_data(ps.beta)),by.x="sample",by.y="ID")
#df<-subset(df,DSS2!="NA")

ggplot(data = df, aes(Axis.1,Axis.2, color = Treatment) )+geom_point() + stat_ellipse() +
  theme_classic()+ facet_wrap(~method,scales = 'free')

#Pull a single diversity metric
dfU<-subset(df,method=="Jaccard")
dfU<-subset(dfU,Treatment!="E4")


ggplot(data = dfU, aes(Axis.1,Axis.2, color = Treatment) )+geom_point() +
  theme_classic()+scale_color_manual(values=c("blue3","yellowgreen","pink2"))+ stat_ellipse()#+geom_line(aes(group=Subject),color="grey",linewidth=.6,alpha=.5)
ggsave("TreatmentEffectBetaDiv.pdf",width=6,height=5)

#Separate ps object by subject
ps.clean.genus2<-subset_samples(ps.clean.genus,DSS2!="NA")
ps.clean.genusS1<-subset_samples(ps.clean.genus2,Microbiome=="S1")
ps.clean.genusS7<-subset_samples(ps.clean.genus2,Microbiome=="S7")
ps.clean.genusS11<-subset_samples(ps.clean.genus2,Microbiome=="S11")

#Plots box plots for all factors
#S1$DSS2<-factor(S1$DSS2,levels = c("Pre","Post"))
#  ggplot(data = S1, aes(x = DSS2, y = Abundance)) +geom_boxplot(outlier.shape  = NA) +
#  geom_point(aes(color = Diet), size=.2) +labs(x = "", y = "Abundance\n") +
#  facet_wrap(~ OTU, scales = "free", nrow = 4) +theme_classic()+
#  geom_line(aes(group=Mouse,color=Diet))+stat_compare_means(label="p")

##############
#This function first filters taxa that are very low abundance (less than .01% in 25% of the samples), 
#tests for differences using wilcoxon test and calculates the log2 fold change

calcWilcox <- function(data) {
  result_list <- list()
  wh1 <- genefilter_sample(data, filterfun_sample(function(x) x > 0.0001), A = 0.25 * nsamples(data))
  data.1 <- prune_taxa(wh1, data)
  S1.1 <- as.data.frame(otu_table(data.1))
  colnames(S1.1) <- gsub("[[:space:]]+", "_", colnames(S1.1))  # Replace spaces with underscores
  colnames(S1.1) <- gsub("[^[:alnum:]_]", "", colnames(S1.1))  # Remove special characters
  taxa <- as.data.frame(names(S1.1))
  names(taxa) <- "taxa"
  S1.1$ID <- rownames(S1.1)
  S1.2 <- merge(metadata_df, S1.1)
  for (i in taxa$taxa) {
    if (length(unique(S1.2[[i]])) > 10) {
      result <- S1.2 %>%
        group_by(Treatment) %>%
        summarise(taxa =i, p_val = wilcox.test(.data[[i]] ~ Day)$p.val,log_fold_change = log2(
          (mean(.data[[i]][Day == "2D"]+ runif(1,.00005,.00009))) / (mean(.data[[i]][Day == "4D"])+ runif(1,.00005,.00009)))
        )
      result$p_adj <- p.adjust(result$p_val, method = "holm")
      result_list[[i]] <- result
    } else {
      cat("Not enough variability in column:", i, "\n")
      cat("Unique values:", unique(S1.2[[i]]), "\n")
    }
  }
  p <- do.call(rbind, result_list)
  return(p)
}

ps.clean.genus1<-subset_samples(ps.clean.genus,Day!="13D")
ps.clean.genus1<-subset_samples(ps.clean.genus1,Day!="3D")
ps.clean.genus1<-subset_samples(ps.clean.genus1,Day!="5D")
ps.clean.genus1<-subset_samples(ps.clean.genus1,Day!="6D")
ps.clean.genus1<-subset_samples(ps.clean.genus1,Day!="7D")
ps.clean.genus1<-subset_samples(ps.clean.genus1,Treatment!="E4")
WilcoxOut1<-calcWilcox(ps.clean.genus1)



#Make a list of taxa that are significantly different Pre/Post in at least one diet
WilcoxOut1.1<-subset(WilcoxOut1,p_val<=.2)
sigtaxaS1<-unique(WilcoxOut1.1$taxa)
WilcoxOut1.2<-WilcoxOut1[which(WilcoxOut1$taxa %in% sigtaxaS1),]

#Make the output into a matrix suitable for heatmapping:
WilcoxOut1.2_E1<-subset(WilcoxOut1.2,Treatment=="E1")
WilcoxOut1.2_E2<-subset(WilcoxOut1.2,Treatment=="E2")
WilcoxOut1.2_E3<-subset(WilcoxOut1.2,Treatment=="E3")


S1LFMat<-as.matrix(cbind(WilcoxOut1.2_E1$log_fold_change,WilcoxOut1.2_E2$log_fold_change,
                         WilcoxOut1.2_E3$log_fold_change))
S1WilMat<-as.matrix(cbind(WilcoxOut1.2_E1$p_val,WilcoxOut1.2_E2$p_val,
                          WilcoxOut1.2_E3$p_val))

colnames(S1WilMat)<-c("E1","E2","E3")
rownames(S1WilMat)<-WilcoxOut1.2_E1$taxa
colnames(S1LFMat)<-c("E1","E2","E3")
rownames(S1LFMat)<-WilcoxOut1.2_E1$taxa

add_stars <- function(value) {
  if (is.numeric(value) && !is.na(value)) {
    if (value < 0.001) {
      return("**")
    } else if (value < 0.05) {
      return("*")
    } else {
      return("")
    }
  } else {
    return("")
  }
}

S1LFSig<-apply(S1WilMat, MARGIN = c(1, 2), FUN = function(x) add_stars(x))


#Make those heat maps yo!
#heatmap.2(S1WilMat,scale="column", col = pal,trace = "none",margins = c(3,11),cellnote = Sig1[2:13],notecol = "black",notecex = 1.5)

pal<-colorpanel(50,"darkblue","white","violetred1")

#Pending how many cells formatting using heatmap2 can be kind of a pain
#But I love the legend+histogram that it makes, some other good packages are superheat, and base R heatmap
#This particular code is specific to subject 1 in my data:
#save it by uncommenting the pdf command below, running the graphing command then dev.off

#pdf("Subject1_LogfoldChangeHeatmap.pdf",width = 8,height=10)
heatmap.2(S1LFMat,scale="none", col = pal,trace = "none", margins=c(3.5,1),
          symkey=FALSE,symbreaks=TRUE, cellnote = S1LFSig,notecol = "black",notecex = 1.5,
          dendrogram='none',density.info='histogram', 
          denscol="black",keysize=6,cexRow = .9,
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(3.5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(0.65, 1.3), lwid=c(1, 1.5, 1))
dev.off()
#Save it:








