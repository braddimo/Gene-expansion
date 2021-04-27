library(phytools)
library(cluster)
library(ape)
library(seriation)
library(RColorBrewer)
library(ggfortify)
library(ggridges)
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

#Compare distributions of transcript size classes to make sure everything is comparable
Abundances <- as.data.frame(read.csv("Size_distribution.csv"))
p <- ggplot(Abundances, aes(x= Size_class, y = Number, group = Species))+
  #  geom_errorbar(aes(ymin=Coral-Coral_std, ymax=Coral+Coral_std), width=.1,position=position_dodge(0.1)) +
  geom_line(aes(color=Species))+ geom_point()
#scale_color_brewer(palette="Paired")+theme_minimal()
p + scale_color_manual(values=c("#040D4F", "#276D96", "#5AC4ED","#21663D","#D7CF64")) + theme_classic()

#Distributions look good so Cafe is ready to run
#Cafe is written in Java, so it needs to be run in a compatible environment
#If using a linux system I would reccoment installing Cafe with Conda so that its ready to go.
#Once Cafe is compiled you are ready to run
#I included all of the input files as well as the script I used in the folder "Cafe"
#To estimate the error model run this:
#python  caferror.py -i ./Cafe_run.sh -dreports/caferror_files  -v 0 -f 1python
#Estimating the error model takes quite awhile so I included the final result in this folder in a file called 'Error_model_scores.txt'
#The best score was 1.220703125e-05 so specify this rate when running Cafe
#Now you can calculate gene family turnover with the scripts called 'Cafe_run_final.sh'
#After running Cafe one species (M.cav) has the highest number of rapidly evolving gene families so re-run cafe to compare if this species should be given its own evolutionary rate
#The online tutorial explains this process super well and I followed it exactly so I would refer to the cafe manual to see how to do this
#The results of the turnover rate comparisons are found in this folder under the title: 
#The log-liklihood of tghe single rate of turnover outpreformed the double rate so stick with a sinlge rate.


####################################################################
################  Go enrichments on Cafe results ###################
####################################################################

#Navigate to the Cafe_results folder
#Use this framework to take the rapidly evolving Orthogroups/Gene families and link them to transcript IDs and GO annotations
#Test for the effects of assembly quality on Cafe results
Model <- read.csv("Rapid_OGs_linear_models.csv")
lm <- lm(Rapid_Ogs~Contig_number_filtered+Busco_missing+Busco_duplicated, data = Model)
summary(lm)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)            -2.765e+02  5.423e+01  -5.099   0.1233  
#Contig_number_filtered  4.211e-02  6.356e-03   6.625   0.0954 .
#Busco_missing          -8.085e-02  1.076e-01  -0.751   0.5899  
#Busco_duplicated       -2.311e+00  2.889e-01  -8.001   0.0792 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 8.522 on 1 degrees of freedom
#Multiple R-squared:  0.9886,	Adjusted R-squared:  0.9544 
#F-statistic:  28.9 on 3 and 1 DF,  p-value: 0.1357

#Looks like the assembly quality overall is not signficiantly influencing the Cafe results, nor are any of the individual parameters

Mcav_OGs <- read.csv("Orthogroups_Mcav.csv")
Mcav_expanding_OGs <- read.csv("Mcav_expanding_orthogroups.csv")
expanding <- inner_join(Mcav_OGs, Mcav_expanding_OGs, by ="Orthogroup")
write.csv(expanding,file = "Mcav_expanding_contigs.csv")
Mcav_contracting_OGS <- read.csv("Mcav_contracting_orthogroups.csv")
Mcav_contractions <- inner_join(Mcav_OGs, Mcav_contracting_OGS, by ="Orthogroup")
write.csv(Mcav_contractions, file = "Mcav_contracting_contigs.csv")
Mcav_annotations <- read.csv("Mcav_annotation.emapper.annotations.csv")
Mcav_expanding_contigs <- read.csv("Mcav_expanding_contigs.csv")
Mcav_contracting_contigs <- read.csv("Mcav_contracting_contigs.csv")
Mcav_expansions <- Mcav_annotations$Mcav_contig%in%Mcav_expanding_contigs$Mcav_contig
Mcav_contractions <- Mcav_annotations$Mcav_contig%in%Mcav_contracting_contigs$Mcav_contig
Mcav_results <- cbind(Mcav_annotations,Mcav_expansions,Mcav_contractions)
write.csv(Mcav_results, file = "Mcav_Annotated_gene_evolution.csv")
#For GO enrichment change all the TRUEs to 1s and all the Falses to 0
#Repeat this workflow for each species
#Time to Run GOMWU script can be found at: https://github.com/z0on/GO_MWU
#For each species take the expansions and contractions and run them seperately, make sure for each to use the entire proteome background
#Parameters are smallest = 50, largest = 0.1 and cluster =0.25

#plot interesting GO terms
Immune_bubble <- read.csv("GO_bubble.csv")
ggplot(Immune_bubble, aes(x=Species, y=Term, size = log_pvalue, color=Color)) +
  theme_minimal()+
  geom_point(alpha=0.7)+
  scale_color_manual(values=c('#425AA8','#6DC282', '#BCB189'))+
  scale_size(range = c(0.1, 10), name="pvalue")


####################################################################
################  Gene expression quantification ###################
####################################################################

#Navigate to the Expression directory
#Expression was quantified with the Salmon in bash using an index size of 31
#This this framework for each species by navigating in thier directory
Mcav_samples <- read.table("Mcav_samples.txt", header = TRUE)
Mcav_samples
Mcav_files <- file.path(Mcav_samples$Mcav_samples, "quant.sf")
#For the species with only 4 samples change this line to ("sample", 1:4)
names(Mcav_files) <- paste0("sample", 1:5)
all(file.exists(Mcav_files))
#This needs to say "TRUE" to continue
Mcav_tx2gene <- read.csv("Mcav_txt2_gene.csv")
head(Mcav_tx2gene)
Mcav_txi <- tximport(Mcav_files, type = 'salmon',tx2gene=Mcav_tx2gene)
write.csv(Mcav_txi, file = "Mcav_counts.csv")
#For downstream analysis only keep the counts and also note that the sample order in the output file is the same order as in the sample read in file

####################################################################
##################### Expression divergence ########################
####################################################################

#First we need to combine all of the counts and normalize in DESEQ2 in order to run the evolutionary Anova model
#For these we are interested in counts not TPM
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
speciesTree <- read.tree("Species_tree.txt")
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
pvalue <- pval
Eve_results <- cbind(rownames(exprMat),results,pvalue)
write.csv(Eve_results, file = "Expression_shifts.csv")
#Not that to convert Beta values to Expression divergence -log10 normlize the beta values, that way the expression divergence is signed and normal.

#Now combine these results with the Cafe results
# Perform Fishers exact test on how many rapdily evolving OG also have significant expression shifts
Expression_shifts <- matrix(c(15,313,41,1660),
                            nrow = 2,
                            dimnames = list(c("Hits","Total"),
                                            c("Hits","Total")))
fisher.test(Expression_shifts, alternative = 'two.sided')

Eve_scatter <- read.csv("Eve_scatterplot.csv")
Eve_family_size <- inner_join(OG_counts,Eve_scatter, by ="OG")
colnames(df_Eve_scatter)
Expression_divergence_lm <- lm(Eve_family_size$Expression_divergence~Eve_family_size$Total)
summary(Expression_divergence_lm)
#Multiple R-squared:  0.002771,	Adjusted R-squared:  0.002588 
#F-statistic: 15.14 on 1 and 5449 DF,  p-value: 0.000101

#The relationship between gene family size class and expression divergence is significant, but the amount of varaince explained is really low

ggscatter(df_Eve_scatter,x= 'Expression_divergence', y='log_pvalue',
          color = 'Class',
          palette = c("aliceblue", "deepskyblue3", "lightgreen"),
          size = (2.0),
)


#Compare the rates of evolution between rapidly evolving immune gene families, rapidly evolving gene families and single copy gene families
Beta_class <- read.csv("Immunity_rapid_Eve.csv")
aov <- aov(log_beta~Ortho_size_class, data = Beta_class)
summary(aov)
TukeyHSD(aov)

p <- ggplot(Beta_class, aes(x=Ortho_size_class, y=Inverted_log_beta, color=Ortho_size_class)) + 
  geom_violin(trim=FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.4)) +
  stat_summary(fun.y=mean, geom="point", size=2, color="cyan") +
  theme_minimal()+
  scale_color_manual(values=c("#33099F", "#120355","#BCB189"))
p