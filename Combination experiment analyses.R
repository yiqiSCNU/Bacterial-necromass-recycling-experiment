
rm(list=ls())


library(vegan)
library(ggplot2)
library(phyloseq)
library(ape)
library(picante)


## Combination necromass experiment
## 1. PCoA and PERMANOVA
## 2. mantel tests and scatter plots (Extended Data Fig. 3)
## 3. NDC richness analyses (Fig.2c)
## 4. Structure equation modeling (Fig.2d; Extended Data Fig. 5)


## 1. PCoA and PERMANOVA

ASV <- read.table("NPS_T6_otu.txt",header = T)
metadata <- read.csv("NPS metadata.csv",row.names = 1,header = T)

metadata$Diversity <- as.factor(metadata$Diversity)

NPS.phylo <- phyloseq(otu_table(ASV,taxa_are_rows = T),sample_data(metadata))
NPS.phylo


bray.A <- vegdist(t(otu_table(NPS.phylo)),method = "bray",binary=FALSE)

ordinate.A <- ordinate(NPS.phylo, "MDS", bray.A)

PCoA.data <- ordinate.A$vectors
dim(PCoA.data)


PCoA.data <- cbind(PCoA.data,metadata)


mycolors_2 <- c("#4DAF4A","#969696","brown","orange","#9400D3","#E41A1C",
                "blue","#FFD700","#EE00EE","black","#FB8072","#63B8FF")



## PCoA for all treatments, and shown in individual PNS level 
## Fig. 2a; Extended Data Fig. 2

p0 <- ggplot(PCoA.data,aes(x=PCoA1,y=PCoA2)) +
  geom_point(aes(color=as.factor(Diversity))) +
  scale_color_manual(values=c("red", "chartreuse4","#9400D3","orange","grey40","dodgerblue")) +
  scale_x_continuous(limits=c(-0.5, 0.38)) +
  scale_y_continuous(limits=c(-0.5, 0.34)) +
  xlab("PCoA1 (21.5%)") + ylab("PCoA2 (4.8%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  ggtitle("A Necromass combination experiment") 
p0




p1 <- ggplot(subset.data.frame(PCoA.data,Diversity==1),aes(x=PCoA1,y=PCoA2)) +
  geom_point(aes(color=Species)) +
  scale_color_manual("Species",values = mycolors_2)+
  scale_x_continuous(limits=c(-0.5, 0.38)) +
  scale_y_continuous(limits=c(-0.5, 0.34)) +
  xlab("PCoA1 (21.5%)") + ylab("PCoA2 (4.8%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  ggtitle("B 1 NPS") 
p1


p2 <- ggplot(subset.data.frame(PCoA.data,Diversity==2),aes(x=PCoA1,y=PCoA2)) +
  geom_point(aes(color=Species)) +
  scale_color_manual("Species",values = mycolors_2)+
  scale_x_continuous(limits=c(-0.5, 0.38)) +
  scale_y_continuous(limits=c(-0.5, 0.34)) +
  xlab("PCoA1 (21.5%)") + ylab("PCoA2 (4.8%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  ggtitle("C 2 NPS") 
p2


p3 <- ggplot(subset.data.frame(PCoA.data,Diversity==3),aes(x=PCoA1,y=PCoA2)) +
  geom_point(aes(color=Species)) +
  scale_color_manual("Species",values = mycolors_2)+
  scale_x_continuous(limits=c(-0.5, 0.38)) +
  scale_y_continuous(limits=c(-0.5, 0.34)) +
  xlab("PCoA1 (21.5%)") + ylab("PCoA2 (4.8%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  ggtitle("D 3 NPS") 
p3


p4 <- ggplot(subset.data.frame(PCoA.data,Diversity==4),aes(x=PCoA1,y=PCoA2)) +
  geom_point(aes(color=Species)) +
  scale_color_manual("Species",values = mycolors_2)+
  scale_x_continuous(limits=c(-0.5, 0.38)) +
  scale_y_continuous(limits=c(-0.5, 0.34)) +
  xlab("PCoA1 (21.5%)") + ylab("PCoA2 (4.8%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  ggtitle("E 4 NPS") 
p4

p5 <- ggplot(subset.data.frame(PCoA.data,Diversity==6 | Diversity==12),aes(x=PCoA1,y=PCoA2)) +
  geom_point(aes(color=Species)) +
  scale_color_manual("Species",values = mycolors_2)+
  scale_x_continuous(limits=c(-0.5, 0.38)) +
  scale_y_continuous(limits=c(-0.5, 0.34)) +
  xlab("PCoA1 (21.5%)") + ylab("PCoA2 (4.8%)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  +
  ggtitle("F 6 and 12 NPS") 
p5


## PERMANOVA, Supplementary Table S3


NPS.phylo <- phyloseq(otu_table(ASV,taxa_are_rows = T),sample_data(metadata))
NPS.phylo

adonis1 <- adonis(bray.A~Species,data=metadata,permutations=9999)
adonis1$aov.tab


T6.phylo.1 <- subset_samples(NPS.phylo,Diversity==1)
T6.phylo.1 <- prune_taxa(taxa_sums(T6.phylo.1)>0,T6.phylo.1)
bray.A <- vegdist(t(otu_table(T6.phylo.1)),method = "bray",binary=FALSE)
adonis1 <- adonis(bray.A~Species,data=subset.data.frame(metadata,Diversity==1),permutations=9999)
adonis1$aov.tab


T6.phylo.2 <- subset_samples(NPS.phylo,Diversity==2)
T6.phylo.2 <- prune_taxa(taxa_sums(T6.phylo.2)>0,T6.phylo.2)
bray.A <- vegdist(t(otu_table(T6.phylo.2)),method = "bray",binary=FALSE)
adonis1 <- adonis(bray.A~Species,data=subset.data.frame(metadata,Diversity==2),permutations=9999)
adonis1$aov.tab


T6.phylo.3 <- subset_samples(NPS.phylo,Diversity==3)
T6.phylo.3 <- prune_taxa(taxa_sums(T6.phylo.3)>0,T6.phylo.3)
bray.A <- vegdist(t(otu_table(T6.phylo.3)),method = "bray",binary=FALSE)
adonis1 <- adonis(bray.A~Species,data=subset.data.frame(metadata,Diversity==3),permutations=9999)
adonis1$aov.tab


T6.phylo.4 <- subset_samples(NPS.phylo,Diversity==4)
T6.phylo.4 <- prune_taxa(taxa_sums(T6.phylo.4)>0,T6.phylo.4)
bray.A <- vegdist(t(otu_table(T6.phylo.4)),method = "bray",binary=FALSE)
adonis1 <- adonis(bray.A~Species,data=subset.data.frame(metadata,Diversity==4),permutations=9999)
adonis1$aov.tab

T6.phylo.6 <- subset_samples(NPS.phylo,Diversity==6)
T6.phylo.6 <- prune_taxa(taxa_sums(T6.phylo.6)>0,T6.phylo.6)
bray.A <- vegdist(t(otu_table(T6.phylo.6)),method = "bray",binary=FALSE)
adonis1 <- adonis(bray.A~Species,data=subset.data.frame(metadata,Diversity==6),permutations=9999)
adonis1$aov.tab

## Mantel tests, scatter plots with moving average

NPS.data <- read.csv("NPS distance data.csv",row.names = 1,header = T)
attach(NPS.data)
mantel(NDC.BC,NPS.BC)
mantel(NDC.BC,NPS.UF)
mantel(NDC.BC,NPS.KO.BC)


## 3.3 * 3.3 inches

library(ggrastr)

p1 <- ggplot(NPS.data,aes(x=NPS.BC,y=NDC.BC)) + 
  geom_jitter_rast(colour="grey",alpha=0.1,width=0.03) + 
  geom_smooth( method = "gam", 
               formula = y ~ s(x, k = 8, bs = "cs"),  
               se = FALSE,       
               color = "darkgrey") +
  xlab("Compositional distances among NPS")+
  ylab("Compositional distances among NDCs")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")+
  ggtitle("a ") 
p1


p2 <- ggplot(NPS.data,aes(x=NPS.UF,y=NDC.BC)) + 
  geom_point_rast(colour="grey",alpha=0.1) + 
  geom_smooth(
    method = "loess",    
    span = 0.5,          
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
p2



p3 <- ggplot(NPS.data,aes(x=NPS.KO.BC,y=NDC.BC)) + 
  geom_point_rast(colour="grey",alpha=0.1) + 
  geom_smooth(
    method = "loess",    
    span = 0.5,          
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
  ggtitle("c") 
p3



## 3. NDC richness analyses

NPS.data <- read.csv("NPS SEM data.csv",header=T)
head(NPS.data)

m1<- lm(NDC.richness~NPS.diversity+I(NPS.diversity^2),data=NPS.data)
summary(m1)

m2 <-lm(NDC.richness~NPS.diversity,data=NPS.data)
summary(m2)

anova(m1,m2)
AIC(m1)
AIC(m2)


## 4. Structure equation modeling

## pairwise correlation

library(psych)
library(piecewiseSEM)

NPS.data <- read.csv("NPS SEM data.csv",header=T)
head(NPS.data)

pairs.panels(
  NPS.data [, 4:8],
  method = "pearson",   
  density = TRUE,       
  ellipses = FALSE,   
  smooth = TRUE,        
  hist.col = "darkgrey",
  rug = FALSE
)


attach(NPS.data)

sem_m0 <- psem(
  lm(KO.richness~NPS.diversity+NPS.mpd),
  lm(metabolic~KO.richness),
  lm(NDC.richness~metabolic),
  data=NPS.data
)


basisSet(sem_m0)

dTable_m0 <- dSep(sem_m0)
dTable_m0
fisherC(dTable_m0)


# potential significant paths are missing,
# model needs to be restructured


sem_m1 <- psem(
  lm(KO.richness~NPS.diversity+NPS.mpd),
  lm(metabolic~KO.richness+NPS.diversity+NPS.mpd),
  lm(NDC.richness~metabolic+NPS.diversity),
  data=NPS.data
)

basisSet(sem_m1)

dTable_m1 <- dSep(sem_m1)
dTable_m1
fisherC(dTable_m1)

coefs(sem_m1)
AIC(sem_m1)   # 30.081
rsquared(sem_m1)
summary(sem_m1)




