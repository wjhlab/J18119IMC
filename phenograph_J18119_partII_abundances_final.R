####Part 2 Density analysis####
####REQUIRED LIBRARIES####
require(scales);require(readxl);require(plyr);require(dplyr);
require(Rphenograph);require(Hmisc); require(ComplexHeatmap); require(pals); require(matrixStats);
require(reshape2); require(ggplot2); require(readxl); require(colorspace); require(ggpubr); require(limma)

#load directory and files
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

metaDataFile = paste0(work,"/Config/J18119_metadata_v2.xlsx")
panelDataFile = paste0(work,"/Config/cleanpanel.xlsx")
dataDirectory = paste0(work,"/Data")

##Set up levels ======
Timepointlevels <- c("PreTx","OnTx")
Tissuelevels <- c("Liver","Lung","Abdomen","Primary")
Tissuelevel <- c("Liver","Lung")
OverallReslevels <- c("PD","SD","PR")
Responselevels <- c("NR","R")
Benefitlevels <- c("NB","B")
Racelevels <- c("C","AA")
samplevels <- c("p30_on_01", "p30_on_02", "p30_pre_01", "p30_pre_02",
                "p24_on_01", "p24_on_02", "p24_pre_01", "p24_pre_02",
                "p28_pre_01", "p28_pre_02", 
                "p56_on_01", "p56_pre_01", "p56_pre_02", "p56_pre_03",
                "p38_on_01", "p38_on_02", "p38_on_03", "p38_pre_01",
                "p20_on_01", "p20_pre_01", "p20_pre_02", "p20_pre_03",
                "p01_on_01", "p01_on_02", "p01_pre_01", "p01_pre_02", "p01_pre_03",
                "p52_pre_01", 
                "p06_pre_01", "p06_pre_02",
                "p26_on_01", "p26_on_02", "p26_on_03",
                "p54_on_01", "p54_on_02", "p54_on_03", "p54_pre_01", "p54_pre_02", "p54_pre_03",
                "p50_on_01", "p50_on_02", "p50_on_03", "p50_pre_01", "p50_pre_02",
                "p45_on_01", "p45_on_02", "p45_on_03", "p45_pre_01", "p45_pre_02", "p45_pre_03",
                "p40_on_01", "p40_on_02", "p40_pre_01", "p40_pre_02",
                "p53_on_01", "p53_on_02", "p53_on_03", "p53_pre_01", "p53_pre_02", "p53_pre_03",
                "p04_on_01", "p04_on_02", "p04_pre_01", "p04_pre_02", "p04_pre_03",
                "p03_pre_01", "p03_pre_02", "p03_pre_03", "p43_on_01",
                "p43_on_02", "p43_on_03", "p43_pre_01", "p43_pre_02", "p43_pre_03",
                "p44_pre_01", "p44_pre_02", 
                "p18_pre_01", "p18_pre_02", 
                "p11_pre_01", "p11_pre_02", 
                "p23_on_01", "p23_on_02", "p23_pre_02", "p23_pre_03", "p23_pre_01",
                "p09_on_01", "p09_on_02", "p09_on_03", "p09_pre_02", "p09_pre_03", "p09_pre_01",
                "p55_pre_01", "p55_pre_02", "p55_pre_03", "p55_pre_04")

caselevels <- c("p30", "p24", "p28", "p56", "p38", "p20", "p01", "p52", "p06", "p26", "p54",
                "p50", "p45", "p40", "p53", "p04", "p03", "p43", "p44", "p18", "p11", "p23",
                "p09", "p55")

##Read-in metadata and clean =======
ifelse(grepl(metaDataFile,pattern='.xlsx'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
md$File_name <- factor(md$File_name)
md$File_order <- factor(md$File_order)
md$Sample_ID <- factor(md$Sample_ID)
md$AvgSample_ID <- factor(md$AvgSample_ID)
md$pubCase <- factor(md$pubCase)
md$Timepoint <- factor(md$Timepoint, levels=Timepointlevels)
md$Tissue <- factor(md$Tissue, levels=Tissuelevels )
md$OverallRes <- factor(md$OverallRes, levels=OverallReslevels)
md$Response <- factor(md$Response, levels=Responselevels)
md$Benefit <- factor(md$Benefit, levels=Benefitlevels)
md$Race <- factor(md$Race, levels=Racelevels)

#LOAD if output is previously saved
output<-readRDS("backup_output_0.RDS")
data_full <- data.frame(output[[1]])
data <- data.matrix(output[[2]])
data01 <- output[[3]]
csv_full <- output[[4]]

##colors
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(data_full$clustercombined2)))
hex <- hue_pal()(9)
colorassignedbroad <- c(rep(hex[1],1), #lym
                        rep(hex[5],1), #other
                        rep(hex[2],3), #Epith
                        rep(hex[3],5), #myeloids
                        rep(hex[4],1), #stroma
                        rep(hex[1],3), #lym
                        rep(hex[10],1)) #hep
clusterlevels=c("Hep","Epi_I", "Epi_II", "Epi_III", "Epi_KI67", "Epi_pSTAT3",
                "Tc","Tc1","TcM","Th","Th1","ThM","Treg","DPT","NKT","B",
                "Mac_I","Mac_II","DC","Neutrophil","Stroma","UA")
clusternames<-clusterlevels
names(colorassigned)<-clusternames

####DATAFRAME SETUP####

#how many of each cell types are there in the dataset?
counts_table <- table(data_full$clustercombined2, data_full$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(props,"props.csv")

ggdf <- melt(data.frame(cluster = rownames(props), props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_ids")
ggdf$sample_id <- factor(ggdf$sample_ids, levels=samplevels)
ggdf$Timepointlevels <- factor(md$Timepoint[match(ggdf$sample_ids,md$Sample_ID)], levels=Timepointlevels)
ggdf$Tissuelevels <- factor(md$Tissue[match(ggdf$sample_ids,md$Sample_ID)], levels=Tissuelevels)
ggdf$Tissuelevel <- factor(md$Tissue[match(ggdf$sample_ids,md$Sample_ID)], levels=Tissuelevel)
ggdf$OverallReslevels <- factor(md$OverallRes[match(ggdf$sample_ids,md$Sample_ID)], levels=OverallReslevels)
ggdf$Responselevels <- factor(md$Response[match(ggdf$sample_ids,md$Sample_ID)], levels=Responselevels)
ggdf$Benefitlevels <- factor(md$Benefit[match(ggdf$sample_ids,md$Sample_ID)], levels=Benefitlevels)
ggdf$Racelevels <- factor(md$Race[match(ggdf$sample_ids,md$Sample_ID)], levels=Racelevels)
ggdf$caselevels <- factor(md$pubCase[match(ggdf$sample_ids,md$Sample_ID)], levels=caselevels)

#### Plots ####
### % CELLS by Timepoint ####
ggp2<-ggplot(ggdf,aes(x=Timepointlevels,y=proportion,fill=Timepointlevels))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2)+
  scale_shape_manual(values=c(1:10,1:18,1:8))+
  facet_wrap(~cluster,ncol=6,scales="free")+
  ylab("% of Cells")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11))+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,size = 0.1,
               position = position_dodge(0.9),
               colour = "black")+
  stat_compare_means(aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
                     method = "t.test",
                     label.x = 1.4,
                     size = 3,         # Increase size of the p-value labels
                     color = "black",  # Change color of the p-value labels
                     hide.ns = TRUE    # Hide non-significant results ("ns")
  )
pdf('WH_Abundance_box_Timepoint.pdf',width=14,height=14)
ggp2
dev.off()
 
### % CELLS by response only for lung and liver and at baseline ####
ggdf_pre <- ggdf[ggdf$Timepointlevels=="PreTx",]

ggp2<-ggplot(ggdf_pre,aes(x=Responselevels,y=proportion,fill=Responselevels))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2, aes(shape=Timepointlevels))+
  scale_shape_manual(values=c(1,19))+
  facet_wrap(~cluster,ncol=6,scales="free")+
  ylab("% of Cells")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11))+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.4,size = 0.1,
               position = position_dodge(0.9),
               colour = "black")+
  stat_compare_means(aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
                     method = "t.test",
                     label.x = 1.4,
                     size = 3,         # Increase size of the p-value labels
                     color = "black",  # Change color of the p-value labels
                     hide.ns = TRUE    # Hide non-significant results ("ns")
  )
pdf('WH_Abundance_box_ResponsePre.pdf',width=14,height=14)
ggp2
dev.off()

### % CELLS at baseline lung vs liver ####
ggdf_tissue_pre <- ggdf_pre[ggdf_pre$Tissuelevels %in% c("Lung","Liver"),]
ggdf_tissue_pre$Tissuelevels <- factor(ggdf_tissue_pre$Tissuelevels, levels=c("Liver","Lung"))

ggp2<-ggplot(ggdf_tissue_pre,aes(x=Tissuelevels,y=proportion,fill=Tissuelevels))+
  geom_boxplot(outlier.shape=NA, lwd=0.25, size=0.25)+
  geom_jitter(width=0.1, size=1, aes(shape=Timepointlevels))+
  scale_shape_manual(values=c(1,19))+
  facet_wrap(~cluster,ncol=6,scales="free")+
  ylab("% of Cells")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
        strip.text.x = element_text(size = 11))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  stat_compare_means(aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
                     method = "t.test",
                     label.x = 1.4,
                     label.y.npc = 1,
                     size = 3,         # Increase size of the p-value labels
                     color = "black",  # Change color of the p-value labels
                     hide.ns = TRUE    # Hide non-significant results ("ns")  
  )
pdf('WH_Abundance_box_TissuelevelsPre.pdf',width=8,height=10)
ggp2
dev.off()

### % Cells by Case Levels - only pre among liver ####
ggdf_pre <- ggdf[ggdf$Timepointlevels == "PreTx",]
ggdf_pre <- ggdf_pre[ggdf_pre$Tissuelevels == "Liver",]

ggp2<-ggplot(ggdf_pre,aes(x=caselevels,y=proportion,fill=caselevels))+
  geom_boxplot(outlier.shape=NA, lwd=0.5)+
  geom_jitter(width=0.2, aes(shape=Timepointlevels))+
  scale_shape_manual(values=c(1,19))+
  facet_wrap(~cluster,ncol=3,scales="free")+
  ylab("% of Cells")+
  theme(#axis.text.x = element_blank(),
    axis.text.y = element_text(size=10, color="black"),
    axis.title.x = element_blank(),
    axis.line.x = element_line(size=0.25, color="black"),
    axis.line.y = element_line(size=0.25, color="black"),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size=12, color="black"),
    strip.background = element_rect(fill=NA),
    strip.text = element_text(size=7, color="black"),
    panel.background = element_rect(fill="white"),
    legend.title = element_blank(),
    legend.key.size = unit(2,'lines'),
    legend.text = element_text(size=8),
    legend.key = element_rect(fill="white"),
    strip.text.x = element_text(size = 11))
pdf('WH_Abundance_box_LiverCasePre.pdf',width=14,height=14)
ggp2
dev.off()






