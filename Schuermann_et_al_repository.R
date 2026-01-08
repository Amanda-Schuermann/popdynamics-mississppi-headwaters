#Authors: Amanda Schuermann, Andrew Hafs
#Date: Nov. 2024-March 2025, updated for vers. 4.5.2 in Dec. 2025
#Publication:  "Population Dynamics of Native Benthic Macroinvertebrates 
#Utilizing the Zebra Mussel (Dreissena polymorpha) as Habitat in the Mississippi River Headwater Region"

#This code contains all visualizations done to explore the data, as well as the code to
#create the figures in the paper. Code corresponding to a figure or table are denoted by four #, so you
#can jump directly to them in the code.Some annotations and explanations are provided for ease of use.
#Happy coding!


### Clear environment and free up RAM ###
rm(list=ls(all=TRUE))

### Call the librarian! ###

#install.packages("vegan")
#install.packages("ellipse")
#install.packages("dplyr")
#install.packages("mvnormtest")
#install.packages('devtools')
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#install.packages("DTK")
#install.packages("C:\\Users\\aschu\\OneDrive\\Documents\\Masters Program\\Pub\\DTK_3.5.tar.gz", repos = NULL)

#This package cannot be loaded directly, must download the .tar file to computer and replace path.
#Use this website: https://ftp.osuosl.org/pub/cran/contrib/main/Archive/DTK/

library(vegan)# metaMDS(), envfit()
library(ellipse)# ordiellipse()
library(dplyr)
library(mvnormtest)
library(devtools)
library(pairwiseAdonis) #pairwiseadonis2
library(DTK)

set.seed(944679)
#Setting seed to provide fixed starting point for NMDS iterations

### Data Loading and Cleaning ###

zebs<-read.csv("C:\\Users\\aschu\\OneDrive\\Documents\\Masters Program\\Research\\Data\\Data_key_taxa_NMDS.csv")
div <- read.csv("C:\\Users\\aschu\\OneDrive\\Documents\\Masters Program\\Research\\Data\\diversity_corr_updated_1.csv")
permanova_test <- read.csv("C:\\Users\\aschu\\OneDrive\\Documents\\Masters Program\\Research\\Data\\diversity_PERMANOVA_upt.csv")
AFDM <- read.csv("C:\\Users\\aschu\\OneDrive\\Documents\\Masters Program\\Research\\Data\\data_AFDM_csv.csv")

tot <- rowSums(zebs[,4:13])
zebs_2 <- zebs[tot > 0, ]
zebs_2 <- zebs_2[-8,]

str(div)
div_final <- div[-6, ]
str(div_final)

par(family = "sans")

names(div_final)[18]<-paste("Sphaeriidae")
names(zebs_2)[13]<-paste("Sphaeriidae")
names(permanova_test)[18]<-paste("Sphaeriidae")


## Prepparing Data for Visualization and Analysis ##

#extract each site
control_site <- div_final[div_final$Site==1,]
#control

site_two <- div_final[div_final$Site==2,]
#downstream Bemidji

site_three <- div_final[div_final$Site==3,]
#Power Dam

site_four <- div_final[div_final$Site==4,]
#Downstream Wolf

site_five <- div_final[div_final$Site==5,]
#Downstream Andrusia

treatment_sites <- div_final[div_final$Site!= 1 ,]
#all treatment sites only, no site 1

#Extract by subseason
sub_1 <- div_final[div_final$Subseason==1,]

sub_2 <- div_final[div_final$Subseason==2,]

sub_3 <- div_final[div_final$Subseason==3,]

#Want to compare if there is a difference between Simpson Diversity Index
#for the treatment sites between both SITE and SUBSEASON

#Visualize diversity index by site and subseason

macro_uninfested <- div_final[div_final$Site == 1, "Macro_Count"]

macro_infested <- div_final[div_final$Site > 1, "Macro_Count"]

shapiro.test(macro_infested) #p<0.05
shapiro.test(macro_uninfested)

wilcox.test(macro_infested,macro_uninfested)

#p<0.05, p=0.008781


## PERMANOVA ##

community_test <- permanova_test[,c(2,9:22)]
dist_community_test <- dist(community_test[,2:11], method="euclidean")
community_test$Site_new <- as.factor(community_test$Site_new)
community_test$Subseason_new <- as.factor(community_test$Subseason_new)
str(community_test)
nrow(community_test)

(community_site_new<-adonis2(dist_community_test~community_test$Subseason_new, data=div_final, permutations=9999, method='euclidean')) #p<0.05
(community_subseason_new<-adonis2(dist_community_test~community_test$Site_new, data=div_final, permutations=9999, method='euclidean')) #p<0.05

test_subseason_new <- permutest(betadisper(vegdist(dist_community_test, method='bray'), div_final$Subseason), permutations=9999)
test_site_new <- permutest(betadisper(vegdist(dist_community_test, method='bray'), div_final$Site), permutations=9999)#p<0.001

#Permutest is a post-PERMANOVA test
#If p<0.05, there may be a confounding factor.
#If p>0.05, there is likely not a confounding factor. 
(test_subseason_new) #p>0.05, not due to dispersion

pairwise_subseason_new <- pairwise.adonis2(dist_community_test~Subseason_new, data=community_test, perm=9999, sim.method='euclidean', p.adjust.m="holm")

(pairwise_subseason_new) 
#IS difference between June and July, and May and July. NO difference between May and June.

## NMDS ##

tot_test <- rowSums(community_test[,2:11])
NMDS_data <- community_test[tot_test > 0, ]

NMDS_test<-metaMDS(NMDS_data[,2:11], k=3, distance = "bray")

NMDS_test

#NMDS_3<-metaMDS(zebs_2[,4:13], k=3)

#NMDS_3

dim(NMDS_test$points)
nrow(NMDS_test$points)
length(zebs_2$Subseason)


par(mfrow=c(1,3))
#Axis 1v2
nmds_test_visual_1<-plot(NMDS_test$points, type="n",xlab="NMDS Axis 1", ylab="NMDS Axis 2")
text(NMDS_test$points[,1],NMDS_test$points[,2],labels=as.character(zebs$Sample_ID),cex=0.85)

take_1 <- c(1,2)

text(NMDS_test,display="species",cex=0.8,col='red', choices=take_1)

ordiellipse(NMDS_test,NMDS_data$Subseason_new,kind="se",conf=0.95,col="blue", label = TRUE, choices=take_1)

#Axis 2v3

take_2 <- c(2,3)

nmds_test_visual_2<-plot(NMDS_test$points, type="n",xlab="NMDS Axis 2", ylab="NMDS Axis 3")

text(NMDS_test$points[,2], NMDS_test$points[,3],labels=as.character(zebs_2$Sample_ID),cex=0.85)

text(NMDS_test,display="species",cex=0.8,col='red')

ordiellipse(NMDS_test,NMDS_data$Subseason_new,kind="se",conf=0.95,col="blue", label = TRUE,choices=take_2)

#Axis 1v3

take_3 <- c(1,3)
nmds_test_visual_3<-plot(NMDS_test$points, type="n",xlab="NMDS Axis 1", ylab="NMDS Axis 3")

text(NMDS_test$points[,1],NMDS_test$points[,3],labels=as.character(zebs$Sample_ID),cex=0.85)

text(NMDS_test,display="species",cex=0.8,col='red', choices=take_3)

ordiellipse(NMDS_test,NMDS_data$Subseason_new,kind="se",conf=0.95,col="blue", label = TRUE, choices=take_3)

#Bubble chart

NMDS_test<-metaMDS(NMDS_data[,2:11], k=3, distance = "bray")
NMDS_test

par(mfrow=c(1,1))

plot(NMDS_test, display="sites", choices = c(1,2))
plot(NMDS_test, display="sites", choices = c(1,3))
plot(NMDS_test, display="sites", choices = c(2,3))

(BugsOrd <- scores(NMDS_test, display="sites"))

plot(BugsOrd)
#Now it's easy to do species overlays... (using cex= to vary size of plotted circles)

library(RColorBrewer)
par(mfrow=c(1,1))
display.brewer.pal(name="Dark2", n=3)

par(mfrow=c(2,5))
size <- 5
par(mar = c(6, 4, 2, 2))
par(family = "sans")

Subseason_1_NMDS <- NMDS_data[NMDS_data$Subseason==1,]

Subseason_2_NMDS <- NMDS_data[NMDS_data$Subseason==2,]

Subseason_3_NMDS <- NMDS_data[NMDS_data$Subseason==3,]

#### Figure 2 ####
par(mfrow=c(2,5))

custom_colors=c("gray80", "gray48", "black")

splt <- "Diptera"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)

splt <- "Amphipoda"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Trichoptera"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Ephemeroptera"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Coleptera"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Odonata"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Plecoptera"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt,col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Gastropoda"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Platyhelminthes"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)


splt <- "Sphaeriidae"
plot(BugsOrd, cex=NMDS_data[,splt]/size,xlim=c(-1.2,1.8),ylim=c(-1,1.2),main=splt, col=custom_colors[as.numeric(NMDS_data$Subseason)], lwd=1.5)

par(mfrow=c(1,1))
plot(1,1)
legend("bottom", legend = c("Subseason 1","Subseason 2","Subseason 3"), col=c("gray80", "gray48", "black"), horiz = TRUE, pch=19)

legend("bottom", legend = c("Subseason 1","Subseason 2","Subseason 3"), fill = unique(NMDS_data$Subseason), horiz = TRUE,inset = c(-1, -0.6), xpd=TRUE)

#Jittered stripchart, do this for site and season

par(mfrow=c(1,3))
plot_sub_1 <- sub_1[,c(9:17)]
stripchart(plot_sub_1, method="jitter", vert=T, pch=1, cex=1, main="Subseason 1", cex.axis=0.6, xlab="Taxa", ylab="Frequency", ylim=c(0,80))

plot_sub_2 <- sub_2[,c(9:17)]
stripchart(plot_sub_2, method="jitter", vert=T, pch=1, cex=1, main="Subseason 2", cex.axis=0.6, xlab="Taxa", ylab="Frequency", ylim=c(0,80))

plot_sub_3 <- sub_3[,c(9:17)]
stripchart(plot_sub_3, method="jitter", vert=T, pch=1, cex=1, main="Subseason 3", cex.axis=0.6, xlab="Taxa", ylab="Frequency", ylim=c(0,80))

#Outlier Removed
par(mfrow=c(1,3))
plot_sub_1 <- sub_1[,c(9:17)]
stripchart(plot_sub_1, method="jitter", vert=T, pch=1, cex=1, main="Subseason 1", cex.axis=0.6, xlab="Taxa", ylab="Frequency", ylim=c(0,50))

plot_sub_2 <- sub_2[,c(9:17)]
stripchart(plot_sub_2, method="jitter", vert=T, pch=1, cex=1, main="Subseason 2", cex.axis=0.6, xlab="Taxa", ylab="Frequency", ylim=c(0,50))

plot_sub_3 <- sub_3[,c(9:17)]
stripchart(plot_sub_3, method="jitter", vert=T, pch=1, cex=1, main="Subseason 3", cex.axis=0.6, xlab="Taxa", ylab="Frequency", ylim=c(0,50))

#Each graph on own
par(mfrow=c(1,1))
stripchart(plot_sub_1, method="jitter", vert=T, pch=1, cex=1, main="Subseason 1", cex.axis=1, xlab="Taxa", ylab="Frequency")
stripchart(plot_sub_2, method="jitter", vert=T, pch=1, cex=1, main="Subseason 2", cex.axis=1, xlab="Taxa", ylab="Frequency")
stripchart(plot_sub_3, method="jitter", vert=T, pch=1, cex=1, main="Subseason 3", cex.axis=1, xlab="Taxa", ylab="Frequency")

#ANOVA
#Does Simpson diversity index vary by site? How about by subseason?

#Site
par(mfrow=c(1,1))
plot(div_final$use_simp~div_final$Site, xlim=c(1,5), ylim=c(0,1.5), xaxp = c(1,5,4), xlab="Site ID", ylab="Simpson Diversity Index", pch=div$Subseason)
legend("topright", legend = paste("Subseason", 1:3), pch = 1:3, bty = "n")

bartlett.test(div_final$use_simp~div_final$Site) #p<0.05, unequal variance

par(mfrow=c(1,5))
shapiro.test(div_final$use_simp[div$Site=="1"])
qqnorm(div_final$use_simp[div_final$Site=="1"], ylab="Test")
qqline(div_final$use_simp[div_final$Site=="1"])
shapiro.test(div_final$use_simp[div_final$Site=="2"])
qqnorm(div_final$use_simp[div_final$Site=="2"])
qqline(div_final$use_simp[div_final$Site=="2"])
shapiro.test(div_final$use_simp[div_final$Site=="3"])
qqnorm(div_final$use_simp[div_final$Site=="3"])
qqline(div_final$use_simp[div_final$Site=="3"])
shapiro.test(div_final$use_simp[div_final$Site=="4"])
qqnorm(div_final$use_simp[div_final$Site=="4"])
qqline(div_final$use_simp[div_final$Site=="4"])
shapiro.test(div_final$use_simp[div_final$Site=="5"])
qqnorm(div_final$use_simp[div_final$Site=="5"])
qqline(div_final$use_simp[div_final$Site=="5"])

#Data all has unequal variance and is non-normal, proceed with oneway test and DTK pairwise comparisons

oneway.test(div_final$use_simp~div_final$Site) #P<0.05, there is at least one group that is significantly different
oneway.test(div_final$use_simp~div_final$Subseason) #P>0.05, there is at least one group that is significantly different

library("DTK")

DTK_result_site <- DTK.test(div_final$use_simp,div_final$Site)
DTK_result_site
par(mfrow=c(1,1))
DTK.plot(DTK_result_site)


par(mfrow=c(1,1))
plot(div_final$use_simp~div_final$Site, xlim=c(1,5), ylim=c(0,1.5), xaxp = c(1,5,4), xlab="Site ID", ylab="Simpson Diversity Index", pch=div$Subseason)
legend("topright", legend = paste("Subseason", 1:3), pch = 1:3, bty = "n")
text(1, 1.1, "a") 
text(2, 1.1, "a") 
text(3, 1.1, "ad")
text(4, 1.1, "bc") 
text(5, 1.1, "bd") 
#subseason markers don't really say anything here?

#### Figure 3 ####
boxplot(control_site$use_simp, site_two$use_simp, site_three$use_simp, site_four$use_simp, site_five$use_simp, names=c("Uninfested Site","Downstream Bemidji", "Downstream Stump", "Downstream Wolf", "Downstream Andrusia"), ylim=c(0,1.3), ylab="Simpson Diversity Index", xlab="Site Names")
text(1, 1.1, "a") 
text(2, 1.1, "a") 
text(3, 1.1, "ad")
text(4, 1.1, "bc") 
text(5, 1.1, "bd") 

#List of sig diff sites:
# 4-1, 5-1, 4-2, 5-2, 4-3

#Subseason

par(mfrow=c(1,1))
plot(div_final$use_simp~div_final$Subseason, xlim=c(1,3), ylim=c(0,1.5), xaxp = c(1,3,2), xlab="Site ID", ylab="Simpson Diversity Index", pch=div$Site)
legend("topright", legend = paste("Site", 1:5), pch = 1:5, bty = "n")

bartlett.test(div_final$use_simp~div_final$Subseason) #p>0.05, equal variance

par(mfrow=c(1,3))
shapiro.test(div_final$use_simp[div$Subseason=="1"])
qqnorm(div_final$use_simp[div_final$Subseason=="1"], ylab="Test")
qqline(div_final$use_simp[div_final$Subseason=="1"])
shapiro.test(div_final$use_simp[div_final$Subseason=="2"])
qqnorm(div_final$use_simp[div_final$Subseason=="2"])
qqline(div_final$use_simp[div_final$Subseason=="2"])
shapiro.test(div_final$use_simp[div_final$Subseason=="3"])
qqnorm(div_final$use_simp[div_final$Subseason=="3"])
qqline(div_final$use_simp[div_final$Subseason=="3"])

#All data is non-normal and has equal variance, proceed with kruskal.test

kruskal.test(div_final$use_simp~div_final$Subseason) #p>0.05, no significant difference between subseasons.

par(mfrow=c(1,1))
plot(div_final$use_simp~div_final$Subseason, xlim=c(1,3), ylim=c(0,1.5), xaxp = c(1,3,2), xlab="Site ID", ylab="Simpson Diversity Index", pch=div$Site)
legend("topright", legend = paste("Site", 1:5), pch = 1:5, bty = "n")

#t-test

control_simpson <- control_site$use_simp

treatment_simpson <- treatment_sites$use_simp

shapiro.test(control_simpson) #p<0.05
shapiro.test(treatment_simpson) #p<0.05
#non normal, use wilcox.test

wilcox.test(control_simpson, treatment_simpson) #p<0.05, populations are not the same

hist(control_simpson, col=rgb(0,0,1,1/4), ylim=c(0,50), xlim=c(0,1.2), xlab="Simpson Diversity Index", main="")
hist(treatment_simpson, col=rgb(1,0,0,1/4), add=T)
legend("topright", legend=c("Uninfested", "Infested"), col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)), fill = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
#have unequal sample sizes, more infested data points that uninfested data points.


#### Figure 4 ####
par(mfrow=c(1,2))
boxplot(control_simpson, treatment_simpson, ylim=c(0,1.5), col=c("gray35","gray65"), ylab="Simpson Diversity Index", names=c("Uninfested Site", "Infested Sites"))
text(2, 1.45, "p<0.001", cex=1.5)
#AFDM

par(mfrow=c(1,1))

str(AFDM)

control_AFDM <- AFDM[AFDM$Site==1,]

treatment_AFDM <- AFDM[AFDM$Site!=1,]

shapiro.test(control_AFDM$percent_carbon) #p>0.05
shapiro.test(treatment_AFDM$percent_carbon) #p<0.05
#non normal, use wilcox.test

wilcox.test(control_AFDM$percent_carbon, treatment_AFDM$percent_carbon) #p>0.05, populations are not the same

boxplot(control_AFDM$percent_carbon, treatment_AFDM$percent_carbon, ylim=c(0,1.3),col=c("gray35","gray65"), ylab="Percent Carbon", names=c("Uninfested Site", "Infested Sites"))
text(2, 1.25, "p=0.9457", cex=1.5)

#### Table 1 ####
## T-tests comparing average counts of each taxonomic group ##

#Diptera
Dipt_count_zebs <- as.integer(treatment_sites$corr_diptera)
Dipt_count_nozebs <- control_site$Diptera
str(Dipt_count_nozebs)
str(Dipt_count_zebs)

mean(Dipt_count_zebs)
mean(Dipt_count_nozebs)

shapiro.test(Dipt_count_zebs) #p<0.05
shapiro.test(Dipt_count_nozebs) #p<0.05
#data is non-normal, proceed with wilcox.test

#Amphipoda

Amph_count_zebs <- as.integer(treatment_sites$corr_amphipoda)
Amph_count_nozebs <- control_site$Amphipoda
str(Amph_count_nozebs)
str(Amph_count_zebs)

shapiro.test(Amph_count_nozebs) #p<0.05
shapiro.test(Amph_count_zebs) #p<0.05
#data s non-normal, wilcox

mean(Amph_count_zebs)
mean(Amph_count_nozebs)

#Trichoptera
Tri_count_zebs <- as.integer(treatment_sites$corr_trichoptera)
Tri_count_nozebs <- control_site$Trichoptera
str(Tri_count_nozebs)
str(Tri_count_zebs)

shapiro.test(Tri_count_nozebs) #p<0.05
shapiro.test(Tri_count_zebs) #p<0.05
#non-normal, wilcox

mean(Tri_count_zebs)
mean(Tri_count_nozebs)

#Ephemeroptera
Ephe_count_zebs <- as.integer(treatment_sites$corr_ephemeroptera)
Ephe_count_nozebs <- control_site$Ephemeroptera
str(Ephe_count_nozebs)
str(Ephe_count_zebs)

shapiro.test(Ephe_count_nozebs) #p<0.05
shapiro.test(Ephe_count_zebs) #p<0.05
#non normal, wilcox

mean(Ephe_count_zebs)
mean(Ephe_count_nozebs)

length(Ephe_count_zebs)
length(Ephe_count_nozebs)

ephe_new_nozebs <- Filter(function(x) x != 0, Ephe_count_nozebs)
(ephe_new_nozebs)
ephe_new_zebs <- Filter(function(x) x != 0, Ephe_count_zebs)
(ephe_new_zebs)

(ephe_wilcox <- wilcox.test(ephe_new_nozebs, ephe_new_zebs))#p>0.05, means are not significantly different


#Coleptera
Cole_count_zebs <- as.integer(treatment_sites$corr_coleptera)
Cole_count_nozebs <- control_site$Coleptera
str(Cole_count_nozebs)
str(Cole_count_zebs)

shapiro.test(Cole_count_nozebs) #p<0.05
shapiro.test(Cole_count_zebs) #p<0.05
#non normal, wilcox

mean(Cole_count_zebs)
mean(Cole_count_nozebs)

#Odonata
Odo_count_zebs <- as.integer(treatment_sites$corr_odonata)
Odo_count_nozebs <- control_site$Odonata
str(Odo_count_nozebs)
str(Odo_count_zebs)

shapiro.test(Odo_count_nozebs) #p<0.05
shapiro.test(Odo_count_zebs) #p<0.05
#non normal, wilcox

mean(Odo_count_zebs)
mean(Odo_count_nozebs)

#Plecoptera
Plec_count_zebs <- as.integer(treatment_sites$corr_plecoptera)
Plec_count_nozebs <- control_site$Plecoptera
str(Plec_count_nozebs)
str(Plec_count_zebs)

shapiro.test(Plec_count_nozebs) #p<0.05
shapiro.test(Plec_count_zebs) #p<0.05
#non normal, wilcox

mean(Plec_count_zebs)
mean(Plec_count_nozebs)

#Gastropoda
Gast_count_zebs <- as.integer(treatment_sites$corr_gastropoda)
Gast_count_nozebs <- control_site$Gastropoda
str(Gast_count_nozebs)
str(Gast_count_zebs)

shapiro.test(Gast_count_nozebs) #p<0.05
shapiro.test(Gast_count_zebs) #p<0.05
#non normal, wilcox

mean(Gast_count_zebs)
mean(Gast_count_nozebs)

#Platyhelminthes
Plat_count_zebs <- as.integer(treatment_sites$corr_platyhelminthes)
Plat_count_nozebs <- control_site$Platyhelminthes
str(Plat_count_nozebs)
str(Plat_count_zebs)

shapiro.test(Plat_count_nozebs) #p<0.05
shapiro.test(Plat_count_zebs) #p<0.05
#non normal, wilcox

mean(Plat_count_zebs)
mean(Plat_count_nozebs)

#Unionidea
Uni_count_zebs <- as.integer(treatment_sites$corr_unionidae)
Uni_count_nozebs <- control_site$Sphaeriidae
str(Uni_count_nozebs)
str(Uni_count_zebs)

shapiro.test(Uni_count_nozebs) #p<0.05
shapiro.test(Uni_count_zebs) #p<0.05
#non normal, wilcox

mean(Uni_count_zebs)
mean(Uni_count_nozebs)

#Tests
(dipt_wilcox <- wilcox.test(Dipt_count_zebs, Dipt_count_nozebs))#p<0.001, means are significantly different***
(amp_wilcox <- wilcox.test(Amph_count_zebs, Amph_count_nozebs))#p<0.05, means are significantly different***
(tri_wilcox <- wilcox.test(Tri_count_zebs, Tri_count_nozebs))#p<0.05, means are significantly different***
(odo_wilcox <- wilcox.test(Odo_count_zebs, Odo_count_nozebs))#p<0.05, means are significantly different***
(ephe_wilcox <- wilcox.test(Ephe_count_zebs, Ephe_count_nozebs))#p>0.05, means are not significantly different
(cole_wilcox <- wilcox.test(Cole_count_zebs, Cole_count_nozebs))#p>0.05, means are significantly different
(plec_wilcox <- wilcox.test(Plec_count_zebs, Plec_count_nozebs))#p<0.05, means are significantly different***
(gast_wilcox <- wilcox.test(Gast_count_zebs, Gast_count_nozebs))#p<0.05, means are significantly different***
(plat_wilcox <- wilcox.test(Plat_count_zebs, Plat_count_nozebs))#p<0.05, means are significantly different***
(uni_wilcox <- wilcox.test(Uni_count_zebs, Uni_count_nozebs))#p<0.05, means are significantly different***

mean_dipt_nozebs <- mean(Dipt_count_nozebs)
mean_dipt_zebs <- mean(Dipt_count_zebs)

mean(Dipt_count_zebs)
mean(Dipt_count_nozebs)
mean(Amph_count_nozebs)
mean(Amph_count_zebs)
mean(Tri_count_nozebs)
mean(Tri_count_zebs)
mean(Ephe_count_zebs)
mean(Ephe_count_nozebs)
mean(Odo_count_nozebs)
mean(Odo_count_zebs)
mean(Cole_count_zebs)
mean(Cole_count_nozebs)
mean(Plec_count_nozebs)
mean(Plec_count_zebs)
mean(Gast_count_nozebs)
mean(Gast_count_zebs)
mean(Plat_count_nozebs)
mean(Plat_count_zebs)
mean(Uni_count_nozebs)
mean(Uni_count_zebs)

(Cole_count_nozebs)
(Cole_count_zebs)

boxplot(Cole_count_nozebs, Cole_count_zebs, Ephe_count_zebs, Ephe_count_nozebs)
hist(Cole_count_zebs, col="blue")
hist(Cole_count_nozebs, add=T, col="red")
hist(Cole_count_nozebs)

stripchart(Cole_count_nozebs, method="jitter", vert=T)
stripchart(Cole_count_zebs, method="jitter", vert=T)

boxplot(Dipt_count_zebs, Dipt_count_nozebs, Amph_count_zebs, Amph_count_nozebs,Tri_count_zebs, Tri_count_nozebs,Ephe_count_zebs, Ephe_count_nozebs,Cole_count_zebs, Cole_count_nozebs,Plec_count_zebs, Plec_count_nozebs,Gast_count_zebs, Gast_count_nozebs,Plat_count_zebs, Plat_count_nozebs, names=c("Inf. Diptera", "Uninf. Diptera", "Inf. Amphipoda", "Uninf. Amphipoda", "Inf. Trichoptera", "Uninf. Trichoptera", "Inf. Ephemeroptera", "Uninf. Ephemeroptera", "Uninf. Coleptera", "Uninf. Coleptera", "Inf. Plecoptera", "Uninf. Plecoptera", "Inf. Gastropoda", "Uninf. Gastropoda", "Inf. Platyhelminthes", "Uninf. Platyhelminthes"),cex.axis=0.45)

boxplot(Dipt_count_zebs, Dipt_count_nozebs, Amph_count_zebs, Amph_count_nozebs,Tri_count_zebs, Tri_count_nozebs,Ephe_count_zebs, Ephe_count_nozebs,Cole_count_zebs, Cole_count_nozebs,Plec_count_zebs, Plec_count_nozebs,Gast_count_zebs, Gast_count_nozebs,Plat_count_zebs, Plat_count_nozebs, names=c("Inf. Diptera", "Uninf. Diptera", "Inf. Amphipoda", "Uninf. Amphipoda", "Inf. Trichoptera", "Uninf. Trichoptera", "Inf. Ephemeroptera", "Uninf. Ephemeroptera", "Uninf. Coleptera", "Uninf. Coleptera", "Inf. Plecoptera", "Uninf. Plecoptera", "Inf. Gastropoda", "Uninf. Gastropoda", "Inf. Platyhelminthes", "Uninf. Platyhelminthes"),cex.axis=0.45, ylim=c(0,6000))

boxplot(Dipt_count_zebs, Dipt_count_nozebs, Amph_count_zebs, Amph_count_nozebs,Tri_count_zebs, Tri_count_nozebs,Ephe_count_zebs, Ephe_count_nozebs, names=c("Inf. Diptera", "Uninf. Diptera", "Inf. Amphipoda", "Uninf. Amphipoda", "Inf. Trichoptera", "Uninf. Trichoptera", "Inf. Ephemeroptera", "Uninf. Ephemeroptera"),cex.axis=0.45, ylim=c(0,6000))

boxplot(Cole_count_zebs, Cole_count_nozebs,Plec_count_zebs, Plec_count_nozebs,Gast_count_zebs, Gast_count_nozebs,Plat_count_zebs, Plat_count_nozebs, names=c("Uninf. Coleptera", "Uninf. Coleptera", "Inf. Plecoptera", "Uninf. Plecoptera", "Inf. Gastropoda", "Uninf. Gastropoda", "Inf. Platyhelminthes", "Uninf. Platyhelminthes"),cex.axis=0.45, ylim=c(0,6000))

#AFDM

str(AFDM)

control_AFDM <- AFDM[AFDM$Site==1,]

treatment_AFDM <- AFDM[AFDM$Site!=1,]

shapiro.test(control_AFDM$percent_carbon) #p>0.05
shapiro.test(treatment_AFDM$percent_carbon) #p<0.05
#non normal, use wilcox.test

wilcox.test(control_AFDM$percent_carbon, treatment_AFDM$percent_carbon) #p>0.05, populations are not the same

