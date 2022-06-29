library(stringr)
library(data.table)
library(reshape2)

male <- fread('~/Documents/PTSD/rna_modifications/comparison_xpore_nanocompore/male/male_xpore_filtered.bed')
names(male) <- c("id","position","V3","V4")
male_xpore <- fread('~/Documents/PTSD/rna_modifications/xpore/male/diffmod/diffmod.filtered.tab')
male_xpore_m6a <- fread('~/Documents/PTSD/rna_modifications/comparison_xpore_nanocompore/male/male_xpore_mods_overlap_rmbase.bed')
male_xpore_m6a <- male_xpore_m6a[,.(id=V1,position=V2)]
male_xpore_m6a[,LABEL:=paste(id,position,sep=":")]

male <- merge(male,male_xpore,by=c("id","position"))

male_rate <- male[,.(id,position,`mod_rate_control-rep1`,`mod_rate_control-rep2`,`mod_rate_control-rep3`,
        `mod_rate_stress-rep1`,`mod_rate_stress-rep2`,`mod_rate_stress-rep3`)]

male_rate_long <- data.table(melt(male_rate,id.vars = c("id","position")))

male_rate_long[str_detect(variable,"control"),condition:="control"]
male_rate_long[str_detect(variable,"stress"),condition:="stress"]

male_rate_long[,sex:="male"]
male_rate_long[,LABEL:=paste(id,position,sep=":")]

male_rate_long[,Type:="other"]
male_rate_long[LABEL %in% male_xpore_m6a$LABEL,Type:="m6A"]

#####################

female <- fread('~/Documents/PTSD/rna_modifications/comparison_xpore_nanocompore/female/female_xpore_filtered.bed')
names(female) <- c("id","position","V3","V4")
female_xpore <- fread('~/Documents/PTSD/rna_modifications/xpore/female/diffmod/diffmod.filtered.tab')

female_xpore_m6a <- fread('~/Documents/PTSD/rna_modifications/comparison_xpore_nanocompore/female/female_xpore_mods_overlap_rmbase.bed')
female_xpore_m6a <- female_xpore_m6a[,.(id=V1,position=V2)]
female_xpore_m6a[,LABEL:=paste(id,position,sep=":")]

female <- merge(female,female_xpore,by=c("id","position"))

female_rate <- female[,.(id,position,`mod_rate_control-rep1`,`mod_rate_control-rep2`,`mod_rate_control-rep3`,
                     `mod_rate_stress-rep1`,`mod_rate_stress-rep2`,`mod_rate_stress-rep3`)]

female_rate_long <- data.table(melt(female_rate,id.vars = c("id","position")))

female_rate_long[str_detect(variable,"control"),condition:="control"]
female_rate_long[str_detect(variable,"stress"),condition:="stress"]

female_rate_long[,sex:="female"]
female_rate_long[,LABEL:=paste(id,position,sep=":")]

female_rate_long[,Type:="other"]
female_rate_long[LABEL %in% female_xpore_m6a$LABEL,Type:="m6A"]

all_rate <- rbind(female_rate_long,male_rate_long)

ggplot(all_rate,aes(condition,value,fill=sex))+
  #facet_wrap(~sex)+
  geom_boxplot()+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  #guides(fill="none")+
  labs(x=NULL,y="Modification rate")

ggplot(all_rate,aes(condition,value,fill=sex))+
  facet_wrap(~Type)+
  geom_boxplot()+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  #guides(fill="none")+
  labs(x=NULL,y="Modification rate")

all_rate_wide <- data.table(dcast(all_rate,id+position+sex+Type~variable,value.var = "value"))
all_rate_wide[,control:=mean(c(`mod_rate_control-rep1`,`mod_rate_control-rep2`,`mod_rate_control-rep3`),na.rm=TRUE),by=list(id,position)]
all_rate_wide[,stress:=mean(c(`mod_rate_stress-rep1`,`mod_rate_stress-rep2`,`mod_rate_stress-rep3`),na.rm=TRUE),by=list(id,position)]
all_rate_wide[,ratio:=stress/control]

ggplot(all_rate_wide,aes(ratio,color=Type))+
  facet_wrap(~sex)+
  geom_density()+
  theme_bw(base_size = 15)+
  theme(aspect.ratio = 1,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(0,3))
  #guides(fill="none")+
  labs(x=NULL,y="Modification rate")

ggplot(all_rate[sex=="male"],aes(condition,value,fill=condition))+
  #facet_wrap(~sex)+
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill="none")+
  labs(x=NULL,y="Modification rate")

ggplot(all_rate[sex=="female"],aes(condition,value,fill=condition))+
  #facet_wrap(~sex)+
  geom_boxplot()+
  theme_bw(base_size = 15)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme(aspect.ratio = 2,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill="none")+
  labs(x=NULL,y="Modification rate")

all_rate[,.(mean(value,na.rm = TRUE)),by=list(condition,sex)]
