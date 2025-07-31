
rm(list=ls())



library(vegan)
library(ggplot2)
library(phyloseq)
library(ape)
library(picante)


## Data analyses for species-rich necromass experiment
## 1. PCoA and PERMANOVA (Fig.4a)
## 2. mantel tests and scatter plots (Extended Data Fig.8)
## 3. NDC richness analyses (Fig.4b; Supplementary Table S8)
## 4. Structure equation modeling (Fig.4c; Extended Data Fig. 9)





## 1. PCoA, distance matrix and PERMANOVA


TPY.meta <- read.table("TPY_metadata.txt",header = T)

TPY.rarefy <- t(read.table("TPY_rarefy_otu.txt",header = T))
dim(TPY.rarefy)
sum(TPY.rarefy)


TPY.KO.data <- t(read.table("TPY_KO_data.txt",header = T))

TPY.sepp <- read.tree(file = "TPY_sepp_tree.nwk")


TPY.phylo1 <- phyloseq(otu_table(TPY.rarefy,taxa_are_rows = FALSE),sample_data(TPY.meta),phy_tree(TPY.sepp))


NPC.phylo <- subset_samples(TPY.phylo1,Type=="NPC")
NDC.phylo <- subset_samples(TPY.phylo1,Type=="NDC")


KO.phylo <- phyloseq(otu_table(TPY.KO.data,taxa_are_rows = FALSE),sample_data(TPY.meta))
KO.NPC <- subset_samples(KO.phylo,Type=="NPC")
KO.NPC <- prune_taxa(taxa_sums(KO.NPC)>0,KO.NPC)




## PCoA

samples.phylo1 <- subset_samples(TPY.phylo1,Type=="NPC" | Type=="NDC")

ordinate1 <- ordinate(samples.phylo1, "MDS", "bray")

p1 <- plot_ordination(samples.phylo1,ordinate1,color="Medium",shape = "Type") +
  geom_point()+
  scale_shape_manual(values = c(17,15)) +
  xlab("PCoA1 (30.4%)") +
  ylab("PCoA2 (14.5%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  ggtitle("") 
p1



## calculate distance matrix for NPC samples


NDC.BC <- vegdist(otu_table(NDC.phylo),method = "bray",binary = F)  ## Bray-Curtis dissimilarities among NDC
NPC.BC <- vegdist(otu_table(NPC.phylo),method = "bray",binary = F)  ## Bray-Curtis dissimilarities among NPC
NPC.UF <- UniFrac(NPC.phylo,weighted = T)   ## UniFrac dissimilarities among NPC
NPC.KO <- vegdist(otu_table(KO.NPC),method = "bray",binary = F)   ## Bray-Curtis of predicted KO profiles among NPC


distance.data <- as.data.frame(cbind(as.vector(NDC.BC),as.vector(NPC.BC),as.vector(NPC.UF),as.vector(NPC.KO)))
colnames(distance.data) <- c("NDC.BC","NPC.BC","NPC.UF","NPC.KO")
write.csv(distance.data,"TPY distance data.csv")


## adonis tests, differences between NPC and NDC
adonis1 <- adonis(vegdist(otu_table(samples.phylo1),method="bray",binary=F)~Type,data=subset.data.frame(TPY.meta,Type=="NPC" | Type=="NDC"),permutations=9999)
adonis1$aov.tab


## adonis tests, differences among three groups of NDC
adonis1 <- adonis(NDC.BC~Medium,data=subset.data.frame(TPY.meta,Type=="NDC"),permutations=9999)
adonis1$aov.tab


## adonis tests, differences among three groups of NPC
adonis1 <- adonis(NPC.BC~Medium,data=subset.data.frame(TPY.meta,Type=="NPC"),permutations=9999)
adonis1$aov.tab





## 2. mantel tests and scatter plots

TPY.dis.data <- read.csv("TPY distance data.csv",row.names = 1,header = T)

mantel(NPC.BC,NDC.BC)
mantel(NPC.UF,NDC.BC)
mantel(NPC.KO,NDC.BC)


p4 <- ggplot(TPY.full.data,aes(x=residual.BC,y=community.BC)) + 
  geom_point(colour="darkgrey") + 
  geom_smooth(
    method = "lm",
    se = FALSE,          
    color = "darkgrey",
    linewidth = 1
  ) +
  xlab("Compositional distances among NPS")+
  ylab("Compositional distances among NDCs")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  ggtitle("a") 
p4


p5 <- ggplot(TPY.full.data,aes(x=residual.UF,y=community.BC)) + 
  geom_point(colour="darkgrey") + 
  geom_smooth(
    method = "lm",   
    se = FALSE,          
    color = "darkgrey",
    linewidth = 1
  ) +
  xlab("Phylogenetic distances among NPS")+
  ylab("Compositional distances among NDCs")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  ggtitle("b") 
p5



p6 <- ggplot(TPY.full.data,aes(x=residual.KO,y=community.BC)) + 
  geom_point(colour="darkgrey") + 
  geom_smooth(
    method = "lm",   
    se = FALSE,         
    color = "darkgrey",
    linewidth = 1
  ) +
  xlab("Functional distances among NPS")+
  ylab("Compositional distances among NDCs")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  ggtitle("c ") 
p6


## 3. NDC richness analyses 
## ANCOVA, Supplementary Table S8

TPY.data <- read.csv("TPY alpha.csv",header = T)
head(TPY.data)

anova(lm(NDC.richness ~ NPC.richness * Medium,data=TPY.data))
##  the interaction term was not significant

## Marginal and conditional effects of NPC.richness and Medium were assessed with two models
anova(lm(NDC.richness ~ NPC.richness + Medium,data=TPY.data))
anova(lm(NDC.richness ~ Medium + NPC.richness,data=TPY.data))

## linear model, Fig.4b
summary(lm(NDC.richness ~ NPC.richness,data=TPY.data))


## 4. SEM 

## TPY SEM

TPY.data <- read.csv("TPY alpha.csv",header = T)
head(TPY.data)

attach(TPY.data)


library(psych)
library(piecewiseSEM)

pairs.panels(
  TPY.data [, 3:6],
  method = "pearson",   
  density = TRUE,       
  ellipses = FALSE,     
  smooth = FALSE,
  lm=TRUE,
  hist.col = "darkgrey",
  rug = FALSE
)



sem_m1 <- psem(
  lm(NDC.richness~NPC.richness+NPC.mpd+NPC.KO),
  lm(NPC.mpd~NPC.richness),
  lm(NPC.KO~NPC.richness),
  data=TPY.data
)

basisSet(sem_m1)

dTable_m1 <- dSep(sem_m1)
dTable_m1
fisherC(dTable_m1)

coefs(sem_m1)
AIC(sem_m1)   # 23.88
rsquared(sem_m1)
summary(sem_m1)


## m2, lowest AIC
sem_m2 <- psem(
  lm(NDC.richness~NPC.richness+NPC.mpd),
  lm(NPC.mpd~NPC.richness),
  lm(NPC.KO~NPC.richness),
  data=TPY.data
)

basisSet(sem_m2)

dTable_m2 <- dSep(sem_m2)
dTable_m2
fisherC(dTable_m2)

coefs(sem_m2)
AIC(sem_m2)   # 24.486
rsquared(sem_m2)
summary(sem_m2)

# 


sem_m3 <- psem(
  lm(NDC.richness~NPC.richness),
  lm(NPC.mpd~NPC.richness),
  lm(NPC.KO~NPC.richness),
  data=TPY.data
)

basisSet(sem_m3)

dTable_m3 <- dSep(sem_m3)
dTable_m3
fisherC(dTable_m3)

coefs(sem_m3)
AIC(sem_m3)   # 29.322
rsquared(sem_m3)
summary(sem_m3)

## Given an AIC difference (¦¤AIC) of 0.6 between model1 vs. model2\
## we selected the more parsimonious model 2



