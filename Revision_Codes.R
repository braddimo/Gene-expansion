library(phytools)
library(cluster)
library(ape)
library(seriation)
library(RColorBrewer)
library(ggfortify)
library(ggridges)
library(ggpubr)
library(phangorn)
library(dplyr)
library(tximportData)
library(tximport)
library(DESeq2)
library(readr)
library(evemodel)

####################################################################
#####################  Prepare data for Cafe #######################
####################################################################
#Make species tree produced from STAG ultrametric for use in Cafe
tree <- read.tree("Species_tree_dates.tre")
ultrametric_tree <- force.ultrametric(tree, method=c("extend"))
text<-write.tree(ultrametric_tree,digits=7)
Ultra_tree <- write.csv(text, file = "Ultrametric_tree.tre")
#Note before Cafe you need to unlabel the tree 
#Read in the OrthoFinder output file and remove all the OGs which don't have atleast 1 transcript in all species
#Compare distributions of transcript size classes to make sure everything is comparable
Abundances <- as.data.frame(read.csv("Size_distribution.csv"))
p <- ggplot(Abundances, aes(x= Size_class, y = Number, group = Species))+
  #  geom_errorbar(aes(ymin=Coral-Coral_std, ymax=Coral+Coral_std), width=.1,position=position_dodge(0.1)) +
  geom_line(aes(color=Species))+ geom_point()
p +scale_colour_viridis(
  alpha = 1,
  begin = 0,
  end = 1,
  direction = 1,
  discrete = TRUE,
  option = "D"
)+theme_classic()

#Distributions look good so Cafe is ready to run
#Cafe is written in Java, so it needs to be run in a compatible environment
#If using a linux system I would reccoment installing Cafe with Conda so that its ready to go.
#Once Cafe is compiled you are ready to run
#I included all of the input files as well as the script I used in the folder "Cafe"
#Now you can calculate gene family turnover with the scripts called 'Cafe.sh'


####################################################################
################  Go enrichments on Cafe results ###################
####################################################################

#Navigate to the Cafe_results folder
#Use this framework to take the rapidly evolving Orthogroups/Gene families and link them to transcript IDs and GO annotations
#Test for the effects of assembly quality on Cafe results
Model <- read.csv("Model_quality_check.csv")
lm <- lm(Rapid_Ogs~Transcript_number+CDAs+Complete.single.copy, data = Model)
summary(lm)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)          150.922600 137.980001   1.094    0.471
#Transcript_number     -0.001564   0.008677  -0.180    0.886
#CDAs                   0.815951   0.250097   3.263    0.189
#Complete.single.copy  -1.322024   0.447600  -2.954    0.208

#Residual standard error: 15.7 on 1 degrees of freedom
#Multiple R-squared:  0.919,	Adjusted R-squared:  0.6759 
#F-statistic:  3.78 on 3 and 1 DF,  p-value: 0.3575

#Looks like the assembly quality overall is not signficiantly influencing the Cafe results, nor are any of the individual parameters

plot_summs(lm, scale = TRUE)

lm1 <- lm(Rapid_contigs~Transcript_number+CDAs+Complete.single.copy, data = Model)
summary(lm1)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)          -3564.4452  2439.9982  -1.461    0.382
#Transcript_number        0.3779     0.1534   2.463    0.246
#CDAs                     4.5685     4.4226   1.033    0.490
#Complete.single.copy    -5.5632     7.9152  -0.703    0.610

#Residual standard error: 277.7 on 1 degrees of freedom
#Multiple R-squared:  0.9084,	Adjusted R-squared:  0.6335 
#F-statistic: 3.305 on 3 and 1 DF,  p-value: 0.3794

plot_summs(lm1, scale = TRUE)

#Navigate into each species folder and run through this pipeline to identify which contigs are in rapidly evolving HTFS to do GO enrichments
#Note that some of file names might change slightly between species
Mcav_OGs <- read.csv("Orthogroups_Mcav.csv")
Mcav_expanding_OGs <- read.csv("Mcav_expanding_OGs.csv")
expanding <- inner_join(Mcav_OGs, Mcav_expanding_OGs, by ="OG")
write.csv(expanding,file = "Mcav_expanding_contigs.csv")
Mcav_contracting_OGS <- read.csv("Mcav_contracting_OGs.csv")
Mcav_contractions <- inner_join(Mcav_OGs, Mcav_contracting_OGS, by ="OG")
write.csv(Mcav_contractions, file = "Mcav_contracting_contigs.csv")
Mcav_annotations <- read.csv("Mcav_annotations.csv")
Mcav_expanding_contigs <- read.csv("Mcav_expanding_contig_IDs.csv")
Mcav_contracting_contigs <- read.csv("Mcav_contracting_contig_IDs.csv")
Mcav_expansions <- Mcav_annotations$Mcav_contig%in%Mcav_expanding_contigs$Mcav_coral_only_denovo
Mcav_contractions <- Mccav_annotations$Mcav_contig%in%Mcavv_contracting_contigs$Mcav_coral_only_denovo
Mcav_results <- cbind(Mcav_annotations,Mcav_expansions)
write.csv(Mcav_results, file = "Mcav_Annotated_gene_evolution.csv")
#For GO enrichment change all the TRUEs to 1s and all the Falses to 0
#Repeat this workflow for each species
#Time to Run GOMWU script can be found at: https://github.com/z0on/GO_MWU
#For each species take the expansions and contractions and run them seperately, make sure for each to use the entire assembly background
#Parameters are smallest = 50, largest = 0.1 and cluster =0.25

#plot interesting GO terms
Immune_bubble <- read.csv("Bubble_plot.csv")
ggplot(Immune_bubble, aes(x=Species, y=Term, size = log_pvalue, color=Color)) +
  theme_minimal()+
  geom_point(alpha=0.7)+
  scale_color_manual(values=c('#425AA8','#6DC282', '#BCB189'))+
  scale_size(range = c(0.1, 1), name="pvalue")


####################################################################
################  Gene expression quantification ###################
####################################################################

#Navigate to the Expression directory
#Expression was quantified with the Salmon in bash using an index size of 31
#This this framework for each species by navigating in thier directory
Mcav_samples <- read.table("Mcav_samples.txt", header = TRUE)
Mcav_samples
Mcav_files <- file.path(Mcav_samples$Samples, "quant.sf")
#For the species with only 4 samples change this line to ("sample", 1:4)
names(Mcav_files) <- paste0("sample", 1:5)
all(file.exists(Mcav_files))
#This needs to say "TRUE" to continue
Mcav_tx2gene <- read.csv("Mcav_txt2_gene.csv")
head(Mcav_tx2gene)
Mcav_txi <- tximport(Mcav_files, type = 'salmon',tx2gene=Mcav_tx2gene)
write.csv(Mcav_txi, file = "Mcav_counts.csv")
#For downstream analysis only keep the counts and also note that the sample order in the output file is the same order as in the sample read in file
#Repeat this step for each species

####################################################################
##################### Expression divergence ########################
####################################################################

#First we need to combine all of the counts and normalize in DESEQ2 in order to run the evolutionary Anova model (EVE)
#For these we are interested in counts not TPM
#Navigate to the folder labeled EVE where all of these files should be
Mcav <- read.csv("Mcav_counts.csv")
Cnat <- read.csv("Cnat_counts.csv")
Orb <- read.csv("Ofav_counts.csv")
Past <- read.csv("Past_counts.csv")
Ssid <- read.csv("Ssid_counts.csv")
Mcav_Cnat <- inner_join(Mcav,Cnat, by ="OG")
MCO <- inner_join(Mcav_Cnat, Orb, by ="OG")
MCOP <- inner_join(MCO,Past,by ="OG")
MCOPS <- inner_join(MCOP,Ssid, by ="OG")
write.csv(MCOPS, file = "Species_expression.csv")

#Read files in to DESEQ2
countData <- as.matrix(read.csv("Species_expression.csv", row.names="OG"))
colnames(countData)
colData <- read.csv("Design.csv",header = TRUE)
colData
#Normalize by species ID
#Because Salmon gives abundances in terms of decimals so you need to round each transcript abundance
dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~Species)
#Remove low abundnace transctipts
dds <- dds[ rowMeans(counts(dds)) > 10, ]
#Noramlize
rld <- rlog(dds)
#Extract normalized counts
rrld <- assay(rld)
write.csv(rrld, file ="Normalized_counts.csv" )
#Read counts in
#I also provided the normalized coutns in case you didn't want to run everything
exprMat <- as.matrix(read.csv("Normalized_counts.csv", row.names = "OG"))
#Same species tree we already used in Cafe
speciesTree <- read.tree("Species_tree_dates.txt")
#Before running EVE make sure that your samples have the species name with and underscore then the sample number and that the species name matches what is in the tree
#Example for Mcav the sample names are 'Mcav_1', 'Mcav_2' and so on
plot(speciesTree)
speciesTree$tip.label
colnames(exprMat)
colSpecies <- sub("_.*$","",colnames(exprMat))
colSpecies
res <- betaSharedTest(tree = speciesTree, gene.data = exprMat, colSpecies = colSpecies)
#Eve takes a really long time to run ~1.5 hours on this size of data so I provide the output
results <- res$indivBetaRes$par
pval = pchisq(res$LRT,df = 1,lower.tail = F)
pvalue <- pval
Eve_results <- cbind(rownames(exprMat),results,pvalue)
write.csv(Eve_results, file = "Expression_shifts.csv")
#Note that to convert Beta values to Expression divergence -log10 normlize the beta values, that way the expression divergence is signed and normal.
#It is also worth looking at the beta values checking which genes have maxed of beta values (99) and considering removing these

Rapid_OGs <- read.csv("Rapid_OGs.csv")
SCO <- read.csv("Single_copy.csv")
Rapid_div <- Expression_divergence$OG%in%Rapid_OGs$OG
SCO_div <- Expression_divergence$OG%in%SCO$OG
Divergence_Class <- cbind(Rapid_div,SCO_div,Expression_divergence)
write.csv(Divergence_Class, file = "Labeled_EVE_results.csv")

#Now combine these results with the Cafe results
# Perform Fishers exact test on how many rapdily evolving OG also have significant expression shifts


Revised_violin <- read.csv("EVE_results_Violing_plot.csv")
p <- ggplot(Revised_violin, aes(x=Class, y=Expression_divergence, color=Class)) + 
  geom_violin(trim=FALSE, size =2) +
  geom_jitter(shape=16, size=2,position=position_jitter(0.4)) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme_minimal()+
  scale_color_manual(values=c("#33099F", "#BCB189"))
p

Rapid_OG_expression <- filter(Divergence_Class, Rapid_div == 'TRUE')
mean(Rapid_OG_expression$Expression_divergence) 
SCO_expression <- filter(Divergence_Class, SCO_div == 'TRUE')
mean(SCO_expression$Expression_divergence)
t.test(Rapid_OG_expression$Expression_divergence,SCO_expression$Expression_divergence)
random_genes <- Divergence_Class[sample(1:nrow(Divergence_Class),217,
                                  replace=FALSE),]
mean(random_genes$Expression_divergence)
t.test(Rapid_OG_expression$Expression_divergence,random_genes$Expression_divergence)

