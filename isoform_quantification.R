setwd('/Users/andre/Documents/GenomeTech/PTSD_experiment')

# Data tables are more efficient and faster than data frames
library(data.table)
library(reshape2)
library(cowplot)
library(ggplot2)
library(DESeq2)

pdf("expected_observed_sequin_abundance.pdf")

# import count matrix as a data table
counts <- fread('all_counts_matrix.tsv')
# import seauin reference file
sequins <- fread('rnasequin_isoforms_2.4.tsv')

# The first thing I did was run DESeq2 (differential expression analysis) the count matrix
# The input to the DESeq2 needs to be formatted in a particular way
# Each sample must be a column and each row a different gene/transcript

# I changed the data so that the ids would be row names
row_names <- counts[,ids]
# and then deleted the ids column
# I also transformed the data into a matrix which is the input format required by DESeq2
cts <- as.matrix(counts[,-1])

# It is also necessary to create a table with the experimental design
# where the row names correspond to the names of each sample
row.names(cts) <- row_names
# and the columns correspond to the different variables in the experimental design
# in this case we only have on variable Condition, with two possible values (control or PTSD)
coldata <- data.frame(condition=as.factor(rep(c("Control","PTSD"),each=3)))
rownames(coldata) <- colnames(counts_matrix)

# I then just give those inputs (count matrix and design table) and transform them into
# the required format read by the Deseq function
# You also need to specify the variable that you are using to compare the samples, 
# which is the "condition" in this case
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# With all that setup you can just run the differential analysis 
dds <- DESeq(dds)

# and use this function to colect the results
res <- results(dds)

# this function here can be use to shrink the effect of non-differential transcripts and
# highlight the differentially expressed ones when you make an MA plot
resLFC <- lfcShrink(dds, coef="condition_PTSD_vs_Control", type="apeglm")

# this command orders the transcripts by p-value 
resOrdered <- res[order(res$pvalue),]
# and you can find how many were significant after applying multiple hypotheses correction
sum(res$padj < 0.05, na.rm=TRUE)

# I subseted only the trancripts that were significant
signif <- res[res$padj < 0.05,]
# I then extract their names
signif_names <- data.table(NAME=rownames(signif))
# I wanted to find how many of those were sequins, so I extracted names that started with R
signif_sequins <- signif_names[startsWith(NAME,"R")]
# There is a difference in the sequin name between the count data and the sequin reference file
# Example: "R1_103_2" "R1_103_2_R1_103"
# So I just adjusted the names from the count data so that they would match
signif_sequins [,ID:=paste(tstrsplit(NAME,"_")[1:3],collapse ="_"),by=NAME]
# then merged the datasets
merge(signif_sequins,sequins,by.x="ID",by.y="NAME")

# I then wanted to separate the sequins between control (same abundance) and differentially expressed.
control_sequin <- sequins[MIX_A==MIX_B,NAME]
diff_sequin <- sequins[!(MIX_A==MIX_B),NAME]

# This functions plot log fold change for trascripts, highlighting the differentially expressed ones
plotMA(res, ylim=c(-2,2))
# you can also plot the shrinkage data
plotMA(resLFC, ylim=c(-2,2))

# You can also the normalized counts for individual transcripts
# I chose to visualize the counts for the 7 sequins that were differentially expressed
plotCounts(dds, gene="R1_103_2_R1_103", intgroup="condition")
plotCounts(dds, gene="R1_13_2_R1_13", intgroup="condition")
plotCounts(dds, gene="R2_150_1_R2_150", intgroup="condition")
plotCounts(dds, gene="R2_150_2_R2_150", intgroup="condition")
plotCounts(dds, gene="R2_53_1_R2_53", intgroup="condition")
plotCounts(dds, gene="R2_53_2_R2_53", intgroup="condition")
plotCounts(dds, gene="R2_7_2_R2_7", intgroup="condition")

# Here I just wanted to see the relationship between the two replicates
# to see if there was a strong correlation
# both between the controls 
ggplot(counts,aes(log2(Control_1_ConditionA_batch1),log2(Control_2_ConditionA_batch1)))+
  geom_point(alpha=0.6)+
  theme_cowplot()+
  geom_abline(slope=1,color="red")+
  theme(aspect.ratio = 1)

# and PTSD condition
ggplot(counts,aes(log2(PTSD_1_ConditionB_batch1),log2(PTSD_2_ConditionB_batch1)))+
  geom_point(alpha=0.6)+
  theme_cowplot()+
  geom_abline(slope=1,color="red")+
  theme(aspect.ratio = 1)

# I just wanted to see the overall count distribution for each sample,
# to visualize any major differences in sequencing depth
# just doing that same correction with the sequin names
counts[startsWith(ids,"R"),NAME:=paste(tstrsplit(ids,"_")[1:3],collapse ="_"),by=ids]
# reformating the data to make the plot, in a long format, so that all samples and their respective trancripts are
# just in two columns instead of 6.
counts <- data.table(melt(counts,id.vars=c("NAME","ids"),variable.name="SAMPLE",value.name = "COUNT"))
# classifying samples as control and PTSD
counts[,GROUP:="control"]
counts[,SAMPLE:=as.character(SAMPLE)]
counts[startsWith(SAMPLE,"PTSD"),GROUP:="PTSD"]

# Making boxplots
ggplot(counts,aes(SAMPLE,log2(COUNT),color=GROUP))+
  geom_boxplot()+
  theme_cowplot()+
  geom_hline(yintercept = log2(10),linetype="dashed")+
  theme(aspect.ratio = 1,axis.text.x = element_text(angle=45,hjust=1))

# I then wanted to specifically look at the sequins
counts_sequins <- counts[startsWith(ids,"R")]

# in evaluate the relationship between expected and observed abundance
sequins <- data.table(melt(sequins,id.vars = c("NAME","LENGTH"),variable.name = "MIX",value.name = "ABUNDANCE"))

# for the control samples samples that have MIX A
control <- merge(counts_sequins[GROUP=="control"],sequins[MIX=="MIX_A"],by="NAME")
# and for the PTSD sampled that had MIX B
ptsd <- merge(counts_sequins[GROUP=="PTSD"],sequins[MIX=="MIX_B"],by="NAME")

table(unique(control$NAME) %in% diff_sequin)
# ploting each group and also showing regression line.
ggplot(control,aes(ABUNDANCE,COUNT))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(MIX~SAMPLE)+
  theme_cowplot()+
  theme(aspect.ratio = 1)


ggplot(ptsd,aes(ABUNDANCE,COUNT))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(MIX~SAMPLE)+
  theme_cowplot()+
  theme(aspect.ratio = 1)

dev.off()


