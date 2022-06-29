library(data.table)
library(VennDiagram)
library(UpSetR)
library(readxl)

overlap <- fread('male_v_female_mod_regions_overlap.tab')

gene_name <- fread('../../../reference_files/gencode.vM25.transcripts.dict.cleaned.tab',header=FALSE)
name <- gene_name[,.(V1,gene_name=V3)]

overlap <- merge(overlap,name)

siter_per_gene <- overlap[,.N,by=gene_name]

ggplot(siter_per_gene,aes(N))+
  geom_histogram()+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 1)+
  labs(x="Number of modification sites",y="Frequency")

setEnrichrSite("Enrichr")

websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021","KEGG_2019_Mouse","MGI_Mammalian_Phenotype_Level_4_2021")
if (websiteLive) {
  enriched <- enrichr(siter_per_gene$gene_name, dbs)
}

bp <- data.table(enriched[["GO_Biological_Process_2021"]])

write.table(bp[-1*log10(Adjusted.P.value)>2],'intersection_biological_process_go_terms.tab',col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


bp_plot <- bp[1:20,.(Term,Adjusted.P.value)]
bp_plot[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
bp_plot <- bp_plot[order(Adjusted.P.value,decreasing = TRUE)]
bp_plot$Term <- factor(bp_plot$Term,levels=bp_plot$Term)

ggplot(bp_plot[-1*log10(Adjusted.P.value)>2],aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightblue")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1/2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=c(0,2,4))+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Biological Process")

mf <- data.table(enriched[["GO_Molecular_Function_2021"]])
mf_plot <- mf[1:20,.(Term,Adjusted.P.value)]
mf_plot[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
mf_plot <- unique(mf_plot[order(Adjusted.P.value,decreasing = TRUE)])
mf_plot$Term <- factor(mf_plot$Term,levels=mf_plot$Term)

ggplot(mf_plot,aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightcoral")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2)+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Molecular Function")


cc <- data.table(enriched[["GO_Cellular_Component_2021"]])
cc_plot <- cc[1:20,.(Term,Adjusted.P.value)]
cc_plot[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
cc_plot <- unique(cc_plot[order(Adjusted.P.value,decreasing = TRUE)])
cc_plot$Term <- factor(cc_plot$Term,levels=cc_plot$Term)

ggplot(cc_plot,aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightgreen")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2)+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Cellular Compartment")


phenotype <- data.table(enriched[["MGI_Mammalian_Phenotype_Level_4_2021"]])
phenotype_plot <- phenotype[1:20,.(Term,Adjusted.P.value)]
phenotype_plot[,Term:=tstrsplit(Term," MP",fixed=TRUE)[1]]
phenotype_plot <- unique(phenotype_plot[order(Adjusted.P.value,decreasing = TRUE)])
phenotype_plot$Term <- factor(phenotype_plot$Term,levels=phenotype_plot$Term)

ggplot(phenotype_plot,aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="orange")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2)+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="MGI_Mammalian_Phenotype_Level_4_2021")

overlap[,ID:=paste("OV",seq(1,.N),sep="_")]

male <- fread('male_mod_regions.bed')
female <- fread('female_mod_regions.bed')

female_nov <- female[!(V4 %in% overlap$V4)]
female_nov <- female_nov[,.(ID=paste("F_NOV",seq(1,.N),sep="_"),Female=1,Male=0)]

male_nov <- male[!(V4 %in% overlap$V8)]
male_nov <- male_nov[,.(ID=paste("M_NOV",seq(1,.N),sep="_"),Female=0,Male=1)]

all_overlap <- rbind(overlap[,.(ID,Female=1,Male=1)],male_nov,female_nov)

upset(all_overlap,sets=c("Male","Female"))

deg <- data.table(read_excel("/Users/andre/Documents/GenomeTech/PTSD_experiment/deg_genes.xlsx"))

male <- merge(male,name)
female <- merge(female,name)

deg_male_mods <- male[gene_name %in% deg[sex=="male" & padj<0.05,external_gene_name]]

write.table(deg_male_mods,'deg_male_mods.tab',col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


deg_female_mods <- female[gene_name %in% deg[sex=="female" & padj<0.05,external_gene_name]]

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021","KEGG_2019_Mouse","MGI_Mammalian_Phenotype_Level_4_2021")
if (websiteLive) {
  enriched <- enrichr(deg_male_mods$gene_name, dbs)
}

bp <- data.table(enriched[["GO_Biological_Process_2021"]])
bp_plot <- bp[1:20,.(Term,Adjusted.P.value)]
bp_plot[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
bp_plot <- bp_plot[order(Adjusted.P.value,decreasing = TRUE)]
bp_plot$Term <- factor(bp_plot$Term,levels=bp_plot$Term)

ggplot(bp_plot,aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightblue")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2)+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Biological Process")
