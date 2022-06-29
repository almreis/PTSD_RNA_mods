library(data.table)
library(readxl)
library(reshape2)

male <- data.table(read_excel("/Users/andre/Documents/GenomeTech/PTSD_experiment/gsea_analysis_biological_process.xlsx",sheet = "Male"))
female <- data.table(read_excel("/Users/andre/Documents/GenomeTech/PTSD_experiment/gsea_analysis_biological_process.xlsx",sheet = "Female"))
all <- rbind(male[,.(Term,Adjusted.P.value,Sex="Male")],female[,.(Term,Adjusted.P.value,Sex="Female")])
all <- data.table(dcast(all,Term~Sex,value.var = "Adjusted.P.value"))
all[,Term:=tstrsplit(Term,"(",fixed=TRUE)[1]]
all[,Mean:=ifelse(Male,Male+Female)/2]
all <- all[order(Mean,Male,Female)]
top10 <- all[1:10]
top10 <- data.table(melt(top10[,.(Term,Male,Female)],id.vars = "Term",variable.name = "Sex",value.name = "Adjuste.P.Value"))
top10$Term <- factor(top10$Term,levels=rev(all[1:10,Term]))


ggplot(top10,aes(Sex,Term,fill=-log10(Adjuste.P.Value)))+
  geom_tile()+
  theme_bw()+
  coord_equal()+
  labs(x=NULL,y=NULL,fill="P-value")
