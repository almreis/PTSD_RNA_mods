library(data.table)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggpubr)
library(ggseqlogo)
library(enrichR)
library(Gviz)
library(RColorBrewer)
library(bedr)

# RNA mods analysis for males
# The same script applies for females 

a <- fread('male_mod_regions.counts.clean.tab')
names(a) <- c('kmer','count')

b_list <- lapply(1:10,function(x) fread(sprintf('male_control%d_regions.counts.clean.tab',x)))
lapply(1:10,function(x) b_list[[x]][,V3:=sprintf("control%d",x)])
b <- rbindlist(b_list)
names(b) <- c('kmer','count','label')

a[,label:="mods"]

dt <- rbind(a,b)

dt[str_detect(kmer,"[AGT][AG]AC[ACT]"),Motif:="DRACH"]

drach_kmers <- unique(dt[Motif=="DRACH",kmer])

drach <- dt[Motif=="DRACH",.(SUM=sum(count)),by=list(label,Motif)]
total <- dt[,.(TOTAL=sum(count)),by=list(label)]

drach <- merge(drach,total)
drach[,FRACTION:=SUM/TOTAL]
drach[,TYPE:="control"]
drach[label=="mods",TYPE:="mods"]

set.seed(4567)
random <- dt[is.na(Motif) & label=="mods"][round(runif(18, min=1, max=dim(dt[is.na(Motif) & label=="mods"])[1]))]
random <- dt[kmer %in% random$kmer]
random_kmers <- random$kmer
random[,Motif:="Random"]
random <- random[,.(SUM=sum(count)),by=list(label,Motif)]
random <- mergrna_mods_male_analysis.Re(random,total)
random[,FRACTION:=SUM/TOTAL]
random[,TYPE:="control"]
random[label=="mods",TYPE:="mods"]

all <- rbind(drach,random)
all.summary <- all[,.(FRACTION=mean(FRACTION)),by=Motif]

ggplot(all,aes(Motif,FRACTION))+
  geom_jitter(aes(color=TYPE),width = 0.1,alpha=0.6,size=3)+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_crossbar(data=all.summary, aes(ymin = FRACTION, ymax=FRACTION),
                size=0.2,col="black", width = .5,alpha=0.6)+
  labs(x="Motif",y="Relative frequency",color=NULL)+
  scale_y_continuous(breaks=c(0.018,0.020,0.022))

ggplot() + geom_logo( drach_kmers, method='p') + theme_logo() + theme(aspect.ratio = 1/3)
ggplot() + geom_logo( random_kmers, method='p' ) + theme_logo() + theme(aspect.ratio = 1/3)

########################

a <- fread('male_mod_regions.counts.clean.tab')
names(a) <- c('kmer','mods')

b_list <- lapply(1:10,function(x) fread(sprintf('male_control%d_regions.counts.clean.tab',x)))
lapply(1:10,function(x) b_list[[x]][,V3:=sprintf("control%d",x)])
b <- rbindlist(b_list)
names(b) <- c('kmer','count','label')
b <- b[,.(control=mean(count)),by=kmer]

dt <- merge(a,b,by="kmer")

dt[str_detect(kmer,"[AGT][AG]AC[ACT]"),Motif:="DRACH"]

ggplot(dt,aes(mods,control,color=Motif))+
  geom_point(alpha=0.6)+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_abline(slope=1,color="red")+
  labs(x="putative mod regions",y="control regions",title="5-mer count")

ggplot(dt,aes(log2(mods),log2(control),color=Motif))+
  geom_point(alpha=0.6)+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_abline(slope=1,color="red")+
  labs(x="putative mod regions (log2)",y="control regions (log2)",title="5-mer count")
  
drach <- data.table(melt(dt[Motif=="DRACH"],id.vars = c("kmer","Motif")))

set.seed(22)
random <- dt[is.na(Motif)][round(runif(18, min=1, max=dim(dt[is.na(Motif)])[1]))]
random[,Motif:="random"]
random <- data.table(melt(random,id.vars = c("kmer","Motif")))

motifs <- rbind(drach,random)

relative_pos <- fread('male_mod_regions.center_mass.features.tab')
relative_pos[,Relative_Pos:=(V2-V6)/(V7-V6)]
relative_pos$V8 <- factor(relative_pos$V8,levels=c("5UTR","CDS","3UTR"))
relative_pos[,Type:=tstrsplit(V4,"_",fixed=TRUE)[1]]

relative_pos[V8=="CDS",Relative_Pos:=Relative_Pos+1]
relative_pos[V8=="3UTR",Rerna_mods_male_analysis.Rlative_Pos:=Relative_Pos+2]


ggplot(relative_pos,aes(V8,fill=V8))+
  geom_bar()+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x=NULL,y="Frequency")+
  guides(fill=FALSE)

feature_per_size <- merge(relative_pos[,.N,by=V8],relative_pos[,.(SIZE=sum(V7-V6)),by=V8])
feature_per_size[,PER_NT:=N/(SIZE/1000)]

ggplot(feature_per_size,aes(V8,PER_NT,fill=V8))+
  geom_bar(stat="identity")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x=NULL,y="Mods/kb")+
  guides(fill=FALSE)

ggplot(relative_pos,aes(Relative_Pos))+
  geom_density(aes(y=..count../max(..count..)),color="darkblue")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1,axis.text.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,2))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  annotate("text", x = 0.5, y = 1.1, label = "5'UTR")+
  annotate("text", x = 1.5, y = 1.1, label = "CDS")+
  annotate("text", x = 2.5, y = 1.1, label = "3'UTR")+
  geom_vline(xintercept =c(1,2),linetype="dashed")+
  labs(y="Relative frequency",x=NULL)

relative_pos <- fread('male_mod_regions.center_mass.features.tab')
relative_pos[,Relative_Pos:=(V2-V6)/(V7-V6)]
relative_pos$V8 <- factor(relative_pos$V8,levels=c("5UTR","CDS","3UTR"))
relative_pos[,Type:=tstrsplit(V4,"_",fixed=TRUE)[1]]
relative_pos[,SIZE:=V7-V6]
step <- 0.025
relative_pos[,Pos.Label:=cut(Relative_Pos,seq(0,1,step),labels=seq(step,1,step))]
relative_pos[,BIN_SIZE:=SIZE*step]
relative_pos_count <- relative_pos[,.N,by=list(V8,Pos.Label,BIN_SIZE)]
relative_pos_count[,DENSITY:=N/(BIN_SIZE/1000)]
relative_pos_count[,Pos.Label:=as.numeric(as.character(Pos.Label))]
relative_pos_count[V8=="CDS",Pos.Label:=Pos.Label+1]
relative_pos_count[V8=="3UTR",Pos.Label:=Pos.Label+2]

max <- quantile(relative_pos_count$DENSITY,0.98)

ggplot(relative_pos_count,aes(Pos.Label,DENSITY))+
  geom_smooth(se=FALSE)+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1,axis.text.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  annotate("text", x = 0.5, y = max, label = "5'UTR")+
  annotate("text", x = 1.5, y = max, label = "CDS")+
  annotate("text", x = 2.5, y = max, label = "3'UTR")+
  geom_vline(xintercept =c(1,2),linetype="dashed")+
  labs(y="Average density (Mods/kb)",x=NULL)
  

#############

relative_pos <- fread('intersection_male_mod_regions_RMBase.features.tab')
relative_pos[,Relative_Pos:=(V2-V6)/(V7-V6)]
relative_pos$V8 <- factor(relative_pos$V8,levels=c("5UTR","CDS","3UTR"))
relative_pos[,Type:=tstrsplit(V4,"_",fixed=TRUE)[1]]

relative_pos[V8=="CDS",Relative_Pos:=Relative_Pos+1]
relative_pos[V8=="3UTR",Relative_Pos:=Relative_Pos+2]

ggplot(relative_pos[Type=="m6A"],aes(Relative_Pos))+
  geom_density(aes(y=..count../max(..count..)),color="darkblue")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1,axis.text.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,2))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  annotate("text", x = 0.5, y = 1.1, label = "5'UTR")+
  annotate("text", x = 1.5, y = 1.1, label = "CDS")+
  annotate("text", x = 2.5, y = 1.1, label = "3'UTR")+
  geom_vline(xintercept =c(1,2),linetype="dashed")+
  labs(y="Relative frequency",x=NULL)
  
m6a_seqs <- fread('intersection_male_mod_regions_RMBase.seqs.tab',header=FALSE)

ggplot() + geom_logo( m6a_seqs$V2) + theme_logo() + theme(aspect.ratio = 1/3)

gene_list <- fread('male_mod_regions.center_mass.features.tab')
gene_list <- unique(gene_list[,.(V1)])

gene_name <- fread('../../../reference_files/gencode.vM25.transcripts.dict.cleaned.tab',header=FALSE)

gene_list <- merge(gene_list,gene_name)

setEnrichrSite("Enrichr")

websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
if (websiteLive) {
  enriched <- enrichr(gene_list$V3, dbs)
}

bp <- data.table(enriched[["GO_Biological_Process_2021"]])

write.table(bp[1:20],'male_biological_process_go_terms.tab',col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")

bp_plot <- bp[1:20,.(Term,Adjusted.P.value)]
bp_plot[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
bp_plot <- bp_plot[order(Adjusted.P.value,decreasing = TRUE)]
bp_plot$Term <- factor(bp_plot$Term,levels=bp_plot$Term)

bp_plot <- bp_plot[order(Adjusted.P.value,decreasing = FALSE)]

ggplot(bp_plot,aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightblue")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1/2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=c(0,2,4))+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Biological Process")

ggplot(bp_plot[1:10],aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightblue")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1/2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=c(0,2,4))+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Biological Process")

##############

mf <- data.table(enriched[["GO_Molecular_Function_2021"]])
mf_plot <- mf[1:20,.(Term,Adjusted.P.value)]
mf_plot[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
mf_plot <- unique(mf_plot[order(Adjusted.P.value,decreasing = TRUE)])
mf_plot$Term <- factor(mf_plot$Term,levels=mf_plot$Term)

ggplot(mf_plot,aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightcoral")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  coord_flip()+rna_mods_male_analysis.R
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Molecular Function")

####################

cc <- data.table(enriched[["GO_Cellular_Component_2021"]])
cc_plot <- cc[1:20,.(Term,Adjusted.P.value)]
cc_plot[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
cc_plot <- unique(cc_plot[order(Adjusted.P.value,decreasing = TRUE)])
cc_plot$Term <- factor(cc_plot$Term,levels=cc_plot$Term)

ggplot(cc_plot,aes(Term,-1*log10(Adjusted.P.value)))+
  geom_bar(stat="identity",fill="lightgreen")+
  geom_hline(yintercept = -1*log10(0.01),linetype="dashed")+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  coord_flip()+
  labs(x=NULL,y="Adjusted p-value (Log10)",title="GO Cellular Compartment")


######################

options(ucscChromosomeNames=FALSE)

ref_dt <- fread('../../../reference_files/gencode.vM25.transcripts.renamed.fa.fai',header=FALSE)
features_dt <- fread('../../../reference_files/gencode.vM25.transcripts.renamed.all_features.protein_coding.bed')
mods <- fread('male_closest_mods_nanocompore_xpore.clean.tab')
nano <- fread('male_nanocompore_filtered.bed')
xp <- fread('male_xpore_filtered.bed')

regions_dt <- fread('male_mod_regions.bed')
rmbase_dt <- fread('RMBase_mm10_all_mods.filtered.tx_coords.bed')
rmbase_intersect_dt <- fread('intersection_male_mod_regions_RMBase.bed')

exons_dt <- fread('male_gencode.vM25.annotation.exons.tx_coords.bed')

transcript <- "ENSMUST00000088666.3"

plot_transcript <- function(transcript,title,mods,rmbase_intersect,rmbase_dt,ref_dt,
                            features_dt,exons_dt,nano,xp,intersection_dt,regions_dt) {
  intersection_dt <- rbind(unique(mods[V4>0 & V4<=10,.(V1,V2,V3=V2+1,V4="nano")]),
                           unique(mods[V4>0 & V4<=10,.(V1,V2=V3,V3=V3+1,V4="xpore")]))
  
  all_colors <- brewer.pal(9, "Set3")
  
  color_features <- all_colors[1:3]
  color_xpore <- all_colors[4]
  color_nano <-  all_colors[5]
  color_overlap <-  all_colors[6]
  color_collapsed <-  all_colors[7]
  color_rmbase <-  all_colors[8:9]
  color_rmbase[1] <- "#e75480"
  
  color_exons <- c("#FED8B1","#FFB862")
  
  alpha_value <- 1
  
  a <- rmbase_intersect[V1==transcript]
  b <- rmbase_dt[V1==transcript]
  b[V2 %in% a$V2,V4:=color_rmbase[1]]
  b[is.na(V4),V4:=color_rmbase[2]]
  
  ref<- with(ref_dt[V1==transcript], GRanges(V1, IRanges(0, V2), "+"))
  ref_track <- GenomeAxisTrack(ref, lwd=4, fontsize=20)
  
  features <-  with(features_dt[V1==transcript], GRanges(V1, IRanges(V2, V3), "+",id=V4))
  features_track <- AnnotationTrack(features, name = "Features", width=15, min.height=1,shape="box",col=NULL,fill=color_features,alpha=alpha_value)
  
  exons <- with(exons_dt[V1==transcript], GRanges(V1, IRanges(V2, V3), "+",id=V4))
  color_exons_rep <- rep(color_exons,length.out=dim(exons_dt[V1==transcript])[1])
  exons_track <- AnnotationTrack(exons, name = "Exons", width=15, min.height=1,shape="box",col=NULL,fill=color_exons_rep,alpha=alpha_value)
  
  nanocompore <- with(unique(nano[V1==transcript,.(V1,V2,V3)]),GRanges(V1, IRanges(V2, V3), "+"))
  nanocompore_track <- AnnotationTrack(nanocompore, name = "Nano", width=15, min.height=1,shape="box",col=NULL,fill=color_nano,alpha=alpha_value)
  
  xpore <- with(unique(xp[V1==transcript,.(V1,V2,V3)]),GRanges(V1, IRanges(V2, V3), "+"))
  xpore_track <- AnnotationTrack(xpore, name = "Xpore", width=15, min.height=1,shape="box",col=NULL,fill=color_xpore,alpha=alpha_value)
  
  intersection <- with(intersection_dt[V1==transcript],GRanges(V1, IRanges(V2, V3), "+",id=V4))
  colors <- intersection$id
  colors[colors=="nano"] <-color_nano
  colors[colors=="xpore"] <- color_xpore
  intersection_track <- AnnotationTrack(intersection, name = "Inter", width=15, min.height=1,shape="box",col=NULL,fill=colors,alpha=alpha_value)
  
  regions <- with(regions_dt[V1==transcript],GRanges(V1, IRanges(V2, V3), "+"))
  regions_track <- AnnotationTrack(regions, name = "Mods", width=15, min.height=1,shape="box",col=NULL,fill=color_collapsed,alpha=alpha_value)
  
  rmbase <- with(b[V1==transcript],GRanges(V1, IRanges(V2, V3), "+"))
  rmbase_track <- AnnotationTrack(rmbase, name = "RMBase", width=15, min.height=1,shape="box",col=NULL,fill=b$V4,alpha=alpha_value)
  
  plotTracks(c(ref_track, exons_track,features_track,xpore_track,nanocompore_track,intersection_track,regions_track,rmbase_track),
             sizes=c(2,1,1,1,1,1,1,1),main=title,stacking="dense")
}

pdf('male_rna_mods_transcript_examples.pdf')

plot_transcript(transcript="ENSMUST00000199856.1",title="Agap3-205",mods,rmbase_intersect_dt,rmbase_dt,ref_dt,features_dt,exons_dt,nano,xp,intersection_dt,regions_dt)

plot_transcript(transcript="ENSMUST00000031843.6",title="Npy",mods,rmbase_intersect_dt,rmbase_dt,ref_dt,features_dt,exons_dt,nano,xp,intersection_dt,regions_dt)

plot_transcript(transcript="ENSMUST00000070375.7",title="Penk",mods,rmbase_intersect_dt,rmbase_dt,ref_dt,features_dt,exons_dt,nano,xp,intersection_dt,regions_dt)

plot_transcript(transcript="ENSMUST00000079278.4",title="Nrsn2",mods,rmbase_intersect_dt,rmbase_dt,ref_dt,features_dt,exons_dt,nano,xp,intersection_dt,regions_dt)

dev.off()

summary_plot <- data.table(V1=c("Xpore","Nanocompore","Overlap","Collapsed mods","RMBase"),V2=rep(0,5))
summary_plot[V1=="Nanocompore",V2:=dim(unique(nano[,.(V1,V2,V3)]))[1]]
summary_plot[V1=="Xpore",V2:=dim(unique(xp[,.(V1,V2,V3)]))[1]]
summary_plot[V1=="Overlap",V2:=dim(mods[V4>0 & V4<=10])[1]]
summary_plot[V1=="Collapsed mods",V2:=dim(regions_dt)[1]]
summary_plot[V1=="RMBase",V2:=dim(rmbase_intersect_dt)[1]]

summary_plot$V1 <-factor(summary_plot$V1,levels=c("RMBase","Collapsed mods","Overlap","Nanocompore","Xpore"))

ggplot(summary_plot,aes(V1,V2,fill=V1))+
  geom_bar(stat="identity",width = 0.7)+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("RMBase"=color_rmbase[1],"Collapsed mods"=color_collapsed,"Overlap"="#BE99A3","Nanocompore"=color_nano,"Xpore"=color_xpore))+
  labs(y="Frequency",x=NULL)+
  geom_text(aes(y=V2+1500,label=prettyNum(V2,big.mark=",",scientific=FALSE)))+
  coord_flip()+
  guides(fill="none")

#####################################

regions_dt <- fread('male_mod_regions.bed')
regions_dt[,CENTER:=round((V3+V2)/2)]
center_mass <- regions_dt[,.(V1,V2=CENTER,V3=CENTER+1,V4)]

exons_dt <- fread('male_gencode.vM25.annotation.exons.tx_coords.bed')
exons_dt[,V5:=V4[1],by=V1]
exons_dt[,V6:=tail(V4,1),by=V1]
exons_dt[,V7:="middle"]
exons_dt[V4==V5,V7:="first"]
exons_dt[V4==V6,V7:="last"]
exons_coords <- exons_dt[,.(V1,V2,V3,V4=paste(V4,V7,sep="_"))]

write.table(center_mass,'male_mod_regions.center_mass.bed',quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")
write.table(exons_coords,'male_gencode.vM25.annotation.exons.tx_coords.v2.bed',quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")

ov <- fread('male_mod_regions.center_mass.exon_intersect.tab')
ov[,EXON:=tstrsplit(V8,"_",fixed=TRUE)[2]]
ov[,TYPE:=tstrsplit(V8,"_",fixed=TRUE)[3]]
ov[TYPE!="first",ACCEPTOR:=V2-V6]
ov[TYPE!="last",DONOR:=V2-V7]

ov$TYPE <- factor(ov$TYPE,levels=c("first","middle","last"))

ggplot(ov,aes(TYPE,fill=TYPE))+
  geom_bar(position = position_dodge())+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 3,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill="none")+
  labs(x="Exon type",y="Frequency")

size <- unique(ov[,.(V1,V8,TYPE,SIZE=V7-V6)])

ggplot(size,aes(log2(SIZE),color=TYPE))+
  geom_density()+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x="Exon length (Log2)",y="Density",color="Exon type")

ov_plot <- rbind(ov[!is.na(DONOR),.(Distance=DONOR,LABEL="donor site")],
      ov[!is.na(ACCEPTOR),.(Distance=ACCEPTOR,LABEL="acceptor site")])

ggplot(ov_plot,aes(log2(abs(Distance)),color=LABEL))+
  geom_density()+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x="Distance to splice site (log2)",y="Density",color=NULL)

ggplot(ov_plot,aes(Distance,fill=LABEL))+
  geom_histogram(alpha=0.6,color="darkgray")+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(-1000,1000))+
  labs(x="Distance to splice site",y="Density",fill=NULL)

ggplot(ov_plot,aes(Distance,color=LABEL))+
  geom_density(alpha=0.6)+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(-1000,1000))+
  labs(x="Distance to splice site",y="Density",color=NULL)

scatter_ov <- rbind(ov[,.(V8,TYPE,SIZE=V7-V6,DISTANCE=abs(ACCEPTOR),LABEL="acceptor site")],
      ov[,.(V8,TYPE,SIZE=V7-V6,DISTANCE=abs(DONOR),LABEL="donor site")])

ggplot(scatter_ov,aes(log2(DISTANCE),log2(SIZE),color=TYPE))+
  geom_point(alpha=0.6)+
  facet_wrap(~LABEL)+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_abline(slope=1,linetype="dashed")+
  labs(x="Distance to splice site (log2)",y="Exon length (Log2)",color="Exon type")

ggplot(ov[TYPE=="middle"],aes(log2(abs(ACCEPTOR)),log2(abs(DONOR))))+
  geom_point(alpha=0.6)+
  theme_bw(base_size=15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_abline(slope=1,linetype="dashed",color="red")+
  labs(x="Distance to acceptor site (log2)",y="Distance to donor site (log2)")

