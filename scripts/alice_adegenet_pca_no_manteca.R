#Author: Alice
#Date: 12/14/22, edited 1/13/23

#NOTE THE THIS IS NOW WITH ONLY MAN1, MAN2, AND MAN3 REMOVED

library(adegenet)
#library(pegas) #can read vcf files; easy to interconvert with adegenet files. but takes forehhhver to convert and thus
library(vcfR) #this one is what I actually used for the reading/conversion

#convert to adegenet file (genind) with vcfR
sunflowers_vcf <- read.vcfR("./data/filtered_hw_no_manteca.recode.vcf")

#need to make a popfile so genind has pop info
popfile <- read.csv("data/popfile_no_manteca.csv") 
popvector <- popfile$pop

sunflowers <- vcfR2genind(sunflowers_vcf, pop=popvector)


#check for missing data
sum(is.na(sunflowers$tab))
#there are 89635 NAs

#replace NAs with mean allele frequency
sunflowers_data <- scaleGen(sunflowers, NA.method="mean")
#check
class(sunflowers_data) #matrix
dim(sunflowers_data) #78 (number of samples) by 11168

#generate PCA
sunflower_pca1 <- dudi.pca(sunflowers_data,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)

#Plotting PCA
pca_plot1 <- s.class(sunflower_pca1$li, pop(sunflowers))
title("PCA of H.winteri and H. annuus \naxes 1-2")

color_scheme <-virid(8) #adegenet version of viridis. 8 because there are 8 pops
nicer_pca_plot <- s.class(sunflower_pca1$li, pop(sunflowers), xax=1,yax=2, col=transp(color_scheme, .6), axesell = FALSE, cstar=0, cpoint=3, grid=FALSE, clabel = 0.65)


#allow PCA to handle missing data instead
x.sunflowers <-tab(sunflowers, freq=TRUE, NA.method="mean")

pca.sunflowers <- dudi.pca(x.sunflowers, center=TRUE, scale=FALSE)

#retained 4 axes, though retaining 3 or 5 could also make sense

pca.sunflowers

#start generating plots
s.class(pca.sunflowers$li, fac=popvector, col=virid(15))

#better plot
s.class(pca.sunflowers$li, fac=popvector, 
        col=transp(virid(15), .6),
        axesell = FALSE, cstar=0, cpoint=3,
        grid=FALSE,
        clabel=0.4)

#verify eigenvalues are variance of components
pca.sunflowers$eig[1]
#57.57624
pc1 <- pca.sunflowers$li[,1]
var(pc1)
#58.32399

#standardize eigenvalues to percentages
eig.perc <- 100*pca.sunflowers$eig/sum(pca.sunflowers$eig)
head(eig.perc)
#First axis=9.430756
#Second axis=4.413757
#third axis=4.323089
#fourth axis=3.205738
#fifth=2.680033
#sixth=2.339150


#even better plot
#better plot
s.class(pca.sunflowers$li, fac=popvector, 
        col=transp(virid(15), .6),
        axesell = FALSE, cstar=0, cpoint=2,
        grid=FALSE,
        clabel=0.5,
        cellipse = FALSE)

#looking at other axes
s.class(pca.sunflowers$li, fac=popvector, 
        col=transp(virid(15), .6),
        axesell = FALSE, cstar=0, cpoint=2,
        grid=FALSE,
        clabel=0.5,
        cellipse = FALSE,
        xax=2,
        yax=3)
s.class(pca.sunflowers$li, fac=popvector, 
        col=transp(virid(15), .6),
        axesell = FALSE, cstar=0, cpoint=2,
        grid=FALSE,
        clabel=0.5,
        cellipse = FALSE,
        xax=3,
        yax=4)

#plotting with ggplot
library(ggplot2)
library(dplyr)
pca_no_petiolaris <- cbind(popfile, pca.sunflowers$li)

#make table with center values for each population
pca_no_petiolaris_summary <- pca_no_petiolaris %>% 
        group_by(pop) %>% 
        summarize(mean_x = mean(Axis1),
        mean_y = mean(Axis2), 
        mean_3 = mean(Axis3), 
        mean_4 = mean(Axis4)) %>% 
        ungroup() 

#add species to summary table
species_vector <- c("HA", "HW", "HA", "HA", "HA", "HA", "HW", "HW")
pca_no_petiolaris_summary$sp <- species_vector

#plot
plot1 <- ggplot(data = pca_no_petiolaris, aes(x = Axis1, y = Axis2, color = sp)) +
        geom_point(alpha=0.8, size=0.2)+
        scale_color_manual(values = c("HW" = "#225ea8",
                                      "HA" ="#31a354")) +
        geom_label(data=pca_no_petiolaris_summary, label.padding = unit(0.05, "lines"), aes(x=mean_x, y=mean_y, label=pop), alpha=0.7, size=2)+
        theme_classic(base_size=5)+
        theme(legend.position="none")+
        xlab("PCA 1 (9.43%)")+
        ylab("PCA 2 (4.41%)")
plot1
#for Mol Ecol paper
ggsave("figures/paper_pca.pdf", plot=plot1, width=80, height=70, units="mm")

plot1_medium <- ggplot(data = pca_no_petiolaris, aes(x = Axis1, y = Axis2, color = sp)) +
        geom_point(alpha=0.8, size=0.8)+
        scale_color_manual(values = c("HW" = "#225ea8",
                                      "HA" ="#31a354")) +
        geom_label(data=pca_no_petiolaris_summary, label.padding = unit(0.05, "lines"), aes(x=mean_x, y=mean_y, label=pop), alpha=0.7, size=2)+
        theme_classic(base_size=8)+
        theme(legend.position="none")+
        xlab("PCA 1 (9.43%)")+
        ylab("PCA 2 (4.41%)")
ggsave("figures/paper_pca_medium.pdf", plot=plot1_medium, width=112, height=90, units="mm")

        
plot2 <-ggplot(data = pca_no_petiolaris, aes(x = Axis2, y = Axis3, color = sp)) +
        geom_point(alpha=0.8)+
        scale_color_manual(values = c("HW" = "#225ea8",
                                      "HA" ="#31a354")) +
        geom_label(data=pca_no_petiolaris_summary, label.padding=unit(0.10, "lines"), aes(x=mean_y, y=mean_3, label=pop), alpha=0.8, size=3)+
        theme_classic()+
        theme(legend.position="none")+
        xlab("PCA 2 (4.41%)")+
        ylab("PCA 3 (4.32%)")
plot2


plot3 <-ggplot(data = pca_no_petiolaris, aes(x = Axis3, y = Axis4, color = sp)) +
        geom_point(alpha=0.8)+
        scale_color_manual(values = c("HW" = "#225ea8",
                                      "HA" ="#31a354")) +
        geom_label(data=pca_no_petiolaris_summary, label.padding=unit(0.10, "lines"),aes(x=mean_3, y=mean_4, label=pop), alpha=0.8, size=3)+
        theme_classic()+
        theme(legend.position="none")+
        xlab("PCA 3 (4.32%)")+
        ylab("PCA 4 (3.21%)")
plot3
