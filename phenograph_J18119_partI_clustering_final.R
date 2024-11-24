####Part I Clustering####
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

## Read-in metadata and clean =======
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

## Make sure all files in metadata present in datadirectory
if(!all(md$File_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])){
  print(paste('ERR: not all filenames in metadata present in data folder - missing',
              md$file_name[!which(data$File_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),
                                                                                     pattern = '.csv')])],'Subsetting...'))
  md <- md[-c(!which(md$File_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.csv')])),]
}

## Read csv into csv_raw =========
csv_raw <- lapply(paste0(dataDirectory,"/",md$File_name),read.csv)
csv_raw_full <- plyr::ldply(csv_raw, rbind)
#csv_raw_full$ImageId <- md$sample_id[match(csv_raw_full$ImageId,md$ImageId)]

#clean csv_raw_full dataframe to csv_full containing analysis markers only
cleanpanel <- read_xlsx(panelDataFile)
colnames(csv_raw_full)[match(cleanpanel$original_names, colnames(csv_raw_full))] <- cleanpanel$clean_names
csv_raw_full <- csv_raw_full[,cleanpanel$clean_names]
panel <- cleanpanel$clean_names[cleanpanel$analysis > 0]
csv_full <- csv_raw_full[,colnames(csv_raw_full) %in% panel]

#sort panels into different categories
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]
#Tcellmarkers <- cleanpanel$clean_names[cleanpanel$tcell == 1]
#Macmarkers <- cleanpanel$clean_names[cleanpanel$mac == 1]
#NKmarkers <- cleanpanel$clean_names[cleanpanel$nk == 1]
#Stromamarkers <- cleanpanel$clean_names[cleanpanel$stroma == 1]

#LOAD if output is previously saved
output<-readRDS("backup_output_0.RDS")
data_full <- data.frame(output[[1]])
data <- data.matrix(output[[2]])
data01 <- output[[3]]
csv_full <- output[[4]]

####SKIP IF LOADING SAVED OUTPUT======

#Cluster heatmap for unannotated clusters
data_full <- csv_full
data <- data.matrix(data_full[,union(subtype_markers,functional_markers)])
data <- asinh(data[, union(subtype_markers,functional_markers)] / 0.8)

#phenograph clustering of data
rng <- colQuantiles(data, probs = c(0.01, 0.99))
data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data01[data01 < 0] <- 0; data01[data01 > 1] <- 1;data01 <-data01[,union(subtype_markers,functional_markers)]

set.seed(1234)
phenographout<-Rphenograph(data01)
data_full$cluster<-factor(membership(phenographout[[2]]))

####PICK UP HERE IF LOADING SAVED OUTPUT====
cluster_mean <- data.frame(data01, cluster = data_full$cluster, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_mat<-as.matrix(cluster_mean[,union(subtype_markers,functional_markers)])

rownames(cluster_mean_mat)<-1:nrow(cluster_mean_mat)

cluster_scaled<-t(scale(t(cluster_mean_mat)))

rownames(cluster_scaled)<-1:nrow(cluster_scaled)

#rename sample_ids to metadata image/sample ids (image ids are automated from histocat and are randomly generated)
cellspersample<-rle(data_full$sample_id)
sample_ids <- as.vector(rep(md$Sample_ID, cellspersample$lengths))
idconversion <- data.frame(new=rle(sample_ids)$values,old=cellspersample$values)
write.csv(idconversion,"idconversion.csv")
data_full$sample_id <- sample_ids

####save as RDS file before any annotations====
backup_output0 <- list(data_full, data, data01, csv_full)
saveRDS(backup_output0, "backup_output0.RDS")

####Annotations begin here====

## Annotation for the original clusters
annotation_row <- data.frame(Cluster = factor(cluster_mean$cluster))
rownames(annotation_row) <- rownames(cluster_mean)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean$cluster)))
names(colorassigned)<- sort(unique(cluster_mean$cluster))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean$cluster),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean$cluster,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data_full$cluster)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

pdf("WH_clusterheatmap_unannotated.pdf",width=10,height=12)
Heatmap(cluster_scaled,
        column_title="Phenograph Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled)*unit(4, "mm"), 
        height = nrow(cluster_scaled)*unit(4, "mm"))
dev.off() 

####DIAGNOSTICS####
## Spot check - number of cells per sample
cell_table <- table(data_full$sample_id)
ggdf <- data.frame(sample_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$case <- factor(md$pubCase[match(ggdf$sample_id,md$Sample_ID)], levels=caselevels)
ggdf$tissue <- factor(md$Tissue[match(ggdf$sample_id,md$Sample_ID)], levels=Tissuelevels)


ggp<-ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) + 
  geom_bar(stat = 'identity') + 
  #geom_text(aes(label = cell_counts), angle = 45, hjust = 0.5, vjust = -0.5, size = 2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('WH_diagnostics_cellcounts.pdf',width=24, height=24);ggp; 
dev.off()

## Multi-dimensional scaling plot to show similarities between samples
## Get the mean marker expression per sample
expr_mean_sample_tbl <- data.frame(sample_id = data_full$sample_id, data) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))
expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$case <- factor(md$pubCase[match(ggdf$sample_id,md$Sample_ID)], levels=caselevels)
ggdf$tissue <- factor(md$Tissue[match(ggdf$sample_id,md$Sample_ID)], levels=Tissuelevels)
ggdf$timepoint <- factor(md$Timepoint[match(ggdf$sample_id,md$Sample_ID)], levels=Timepointlevels)
ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = tissue, shape = timepoint)) +
  #ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = tumor, shape = timegroup)) +
  geom_point(size = 2.5) +
  #geom_text(aes(label = patient_id)) +
  theme_bw()+
  theme(#plot.background = element_rect(fill="white"),
        #panel.background = element_rect(fill="white"),
        #panel.grid = element_blank(),
        #axis.line = element_line(color="white"),
        #axis.text = element_text(color="white"),
        #axis.title = element_text(color="white"),
        legend.background = element_rect(fill="white"),
        #legend.text = element_text(color="white"),
        legend.key = element_rect(fill="white"))
pdf('WH_diagnostics_MDS_tissue.pdf',width=6, height=6);ggp; dev.off()

#cluster heatmap for merged annotations ===========
clusterMergeFile = paste0(work,"/Config/J18119_unverified_merged.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("UA",
                "Epith",
                "Hep",
                "Stroma",
                "Mac_I",
                "Th",
                "Mac_II",
                "DC",
                "Neutrophil",
                "Tc")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(data_full$cluster, cluster_merging$original_cluster)
data_full$cluster1m <- cluster_merging$new_cluster[mm1]

cluster_mean_merged <- data.frame(data01, cluster = data_full$cluster1m, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))

cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])

cluster_scaled_merged<-t(scale(t(cluster_mean_merged_mat)))

rownames(cluster_scaled_merged)<-1:nrow(cluster_scaled_merged)

## Annotation for the merged clusters

if(!is.null(clusterMergeFile)){
  ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
  cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
  annotation_row$Merged <- cluster_merging$new_cluster
  color_clusters2 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Merged))
  names(color_clusters2) <- levels(cluster_merging$new_cluster)
  annotation_colors$Merged <- color_clusters2
}

## Colors for the heatmap
legend_breaks = seq(from = 0, to = 1, by = 0.2)
####
clusternames<-clusterlevels
####
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))

names(colorassigned)<-clusternames
####
rownames(cluster_scaled_merged)<-cluster_mean_merged$cluster

color_list = list(clusters=colorassigned)

color_list_byoriginal = colorassigned[match(unique(cluster_merging$new_cluster),names(colorassigned))]

cp<-rowAnnotation(#clusters=clusternames,
  col=color_list,
  gp = gpar(col = "white", lwd = .5),
  counts= anno_barplot(
    as.vector(table(data_full$cluster1m)),
    gp = gpar(fill=colorassigned),
    border = F,
    bar_width = 0.75, 
    width = unit(2,"cm")))

pdf("WH_clusterheatmap_unverified_merged2.pdf",width=10,height=4)
Heatmap(cluster_scaled_merged[clusterlevels,],
        column_title="Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled_merged)*unit(4, "mm"), 
        height = nrow(cluster_scaled_merged)*unit(4, "mm"))
dev.off()

####Subcluster T cells====
data01_T<-as.data.frame(data01)
data01_T$cluster1m <- data_full$cluster1m
data01_T_subset <- data01_T[data01_T$cluster1m %in% c("Th","Tc"),]
set.seed(1234)
phenographout<-Rphenograph(as.matrix(data01_T_subset[,union(subtype_markers,functional_markers)]))
data01_T_subset$clusterT<-factor(membership(phenographout[[2]]))

cluster_mean_T <- data.frame(data01_T_subset[,union(subtype_markers,functional_markers)], clusterT = data01_T_subset$clusterT, check.names = FALSE) %>%
  group_by(clusterT) %>% summarize_all(list(mean))
cluster_mean_T_mat<-as.matrix(cluster_mean_T[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_T_mat)<-1:nrow(cluster_mean_T_mat)
cluster_scaled_T<-t(scale(t(cluster_mean_T_mat)))
rownames(cluster_scaled_T)<-1:nrow(cluster_scaled_T)

annotation_row <- data.frame(Cluster = factor(cluster_mean_T$clusterT))
rownames(annotation_row) <- rownames(cluster_mean_T)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean_T$clusterT)))
names(colorassigned)<- sort(unique(cluster_mean_T$clusterT))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean_T$clusterT),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean_T$clusterT,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data01_T_subset$clusterT)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

tcellmarkers<-c("CD4","CD3","CD45RO","CD68","CD8","CXCR3","TOX","LAG3","PD1","CD137","FOXP3","CD45RA","CD20","KI67","GZMB","CD57")
pdf("WH_clusterheatmap_unannotated_T.pdf",width=10,height=12)
Heatmap(cluster_scaled_T[,tcellmarkers],
        column_title="Phenograph Clusters T",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled_T)), to = round(max(cluster_scaled_T)))),
        width = ncol(cluster_scaled_T[,tcellmarkers])*unit(4, "mm"), 
        height = nrow(cluster_scaled_T)*unit(4, "mm"))
dev.off() 

clusterMergeFileT = paste0(work,"/Config/J18119_merged_T.xlsx") #create dummy merger numbers prior to annotation
cluster_mergingT <- read_excel(clusterMergeFileT)
clusterlevelsT=c("Tc",
                 "Tc1",
                 "TcM",
                 "Th",
                 "Th1",
                 "ThM",
                 "Treg",
                 "DPT",
                 "NKT",
                 "B",
                 "Mac_I")

mm2 <- match(data01_T_subset$clusterT, cluster_mergingT$original_cluster)
data01_T_subset$cluster1mT <- cluster_mergingT$new_cluster[mm2]

##heatmap just for the T cells after annotating with original expressions
tcellmarkers2<-c("CD4","CD3","CD45RO","CD8","CXCR3","TOX","LAG3","PD1","CD137","FOXP3","CD45RA","KI67","GZMB")
cluster_mean_T <- data.frame(data01_T_subset[,union(subtype_markers,functional_markers)], clusterT = data01_T_subset$cluster1mT, check.names = FALSE) %>%
  group_by(clusterT) %>% summarize_all(list(mean))
cluster_mean_T_mat<-as.matrix(cluster_mean_T[cluster_mean_T$clusterT%nin%c("B","Mac_I","NKT","DPT"),union(subtype_markers,functional_markers)])
rownames(cluster_mean_T_mat)<-cluster_mean_T$clusterT[cluster_mean_T$clusterT%nin%c("B","Mac_I","NKT","DPT")]
cluster_scaled_T<-t(scale(t(cluster_mean_T_mat)))
pdf("WH_clusterheatmap_annotated_T.pdf",width=10,height=12)
pheatmap(cluster_scaled_T[,tcellmarkers2],
         scale="column",
         color = rev(brewer.rdbu(100)),
         cellwidth=10,
         cellheight=10,
         border_color = NA,
         treeheight_row = 20,
         treeheight_col = 20)
dev.off() 

##add in T cell annotations
data_full$clustercombined <- data_full$cluster1m
data_full[rownames(data01_T_subset),]$clustercombined <- data01_T_subset$cluster1mT

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(union(clusterlevels, clusterlevelsT)))
clusternames<-union(clusterlevels, clusterlevelsT)
names(colorassigned)<-clusternames
cluster_mean_merged <- data.frame(data01, cluster = data_full$clustercombined, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])
cluster_scaled_merged<-t(scale(t(cluster_mean_merged_mat)))
rownames(cluster_scaled_merged)<-cluster_mean_merged$cluster

color_list = list(clusters=colorassigned)

cp<-rowAnnotation(#clusters=clusternames,
  col=color_list,
  gp = gpar(col = "white", lwd = .5),
  counts= anno_barplot(
    as.vector(table(data_full$clustercombined)),
    gp = gpar(fill=colorassigned),
    border = F,
    bar_width = 0.75, 
    width = unit(2,"cm")))

pdf("WH_clusterheatmap_combined.pdf",width=10,height=6)
Heatmap(cluster_scaled_merged,
        column_title="Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = F,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled_merged)*unit(4, "mm"), 
        height = nrow(cluster_scaled_merged)*unit(4, "mm"))
dev.off()

####Subcluster Epithelial cells====
data01_epi<-as.data.frame(data01)
data01_epi$cluster1m <- data_full$cluster1m
data01_epi_subset <- data01_epi[data01_epi$cluster1m %in% c("Epith"),]
set.seed(1234)
phenographout<-Rphenograph(as.matrix(data01_epi_subset[,union(subtype_markers,functional_markers)]))
data01_epi_subset$clusterepi<-factor(membership(phenographout[[2]]))

cluster_mean_epi <- data.frame(data01_epi_subset[,union(subtype_markers,functional_markers)], clusterepi = data01_epi_subset$clusterepi, check.names = FALSE) %>%
  group_by(clusterepi) %>% summarize_all(list(mean))
cluster_mean_epi_mat<-as.matrix(cluster_mean_epi[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_epi_mat)<-1:nrow(cluster_mean_epi_mat)
cluster_scaled_epi<-t(scale(t(cluster_mean_epi_mat)))
rownames(cluster_scaled_epi)<-1:nrow(cluster_scaled_epi)

annotation_row <- data.frame(Cluster = factor(cluster_mean_epi$clusterepi))
rownames(annotation_row) <- rownames(cluster_mean_epi)
color_clusters1 <- kovesi.rainbow_bgyrm_35_85_c69(nlevels(annotation_row$Cluster))
names(color_clusters1) <- levels(annotation_row$Cluster)
annotation_colors <- list(Cluster = color_clusters1)

legend_breaks = seq(from = 0, to = 1, by = 0.2)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_mean_epi$clusterepi)))
names(colorassigned)<- sort(unique(cluster_mean_epi$clusterepi))
color_list = list(clusters=colorassigned)
color_list_byoriginal = colorassigned[match((cluster_mean_epi$clusterepi),names(colorassigned))]

rAbar<-rowAnnotation(clusters=cluster_mean_epi$clusterepi,
                     col=color_list,
                     gp = gpar(col = "white", lwd = .5),
                     counts= anno_barplot(
                       as.vector(table(data01_epi_subset$clusterepi)),
                       gp = gpar(fill=colorassigned),
                       border = F,
                       bar_width = 0.75, 
                       width = unit(2,"cm")))

epicellmarkers<-c("pSTAT3","CDX2","PanCK","KI67","HLADR","CD86","PDL1","HepPar","FGL1","ADAM10","ADAM17","VIM")
pdf("WH_clusterheatmap_unannotated_epi.pdf",width=10,height=12)
Heatmap(cluster_scaled_epi[,epicellmarkers],
        column_title="Phenograph Clusters Epi",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = NA,
        rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = rAbar,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled_epi)), to = round(max(cluster_scaled_epi)))),
        width = ncol(cluster_scaled_epi[,epicellmarkers])*unit(4, "mm"), 
        height = nrow(cluster_scaled_epi)*unit(4, "mm"))
dev.off() 

clusterMergeFileEpi = paste0(work,"/Config/J18119_merged_Epi2.xlsx") #create dummy merger numbers prior to annotation
cluster_mergingEpi <- read_excel(clusterMergeFileEpi)
clusterlevelsEpi=c("Epi_I",
                   "Epi_II",
                   "Epi_III",
                   "Epi_KI67",
                   "Epi_pSTAT3",
                   "UA")

mm3 <- match(data01_epi_subset$clusterepi, cluster_mergingEpi$original_cluster)
data01_epi_subset$cluster1mEpi <- cluster_mergingEpi$new_cluster[mm3]

##heatmap just for the T cells after annotating with original expressions
epicellmarkers2<-c("pSTAT3","CDX2","PanCK","KI67","HLADR","PDL1","FGL1","ADAM10","ADAM17")
cluster_mean_epi <- data.frame(data01_epi_subset[,union(subtype_markers,functional_markers)], clusterepi = data01_epi_subset$cluster1mEpi, check.names = FALSE) %>%
  group_by(clusterepi) %>% summarize_all(list(mean))
cluster_mean_epi_mat<-as.matrix(cluster_mean_epi[cluster_mean_epi$clusterepi%nin%c("UA"),union(subtype_markers,functional_markers)])
rownames(cluster_mean_epi_mat)<-cluster_mean_epi$clusterepi[cluster_mean_epi$clusterepi%nin%c("UA")]
cluster_scaled_epi<-t(scale(t(cluster_mean_epi_mat)))
pdf("WH_clusterheatmap_annotated_epi.pdf",width=10,height=12)
pheatmap(cluster_scaled_epi[,epicellmarkers2],
         scale="column",
         color = rev(brewer.rdbu(100)),
         cellwidth=10,
         cellheight=10,
         border_color = NA,
         treeheight_row = 20,
         treeheight_col = 20)
dev.off() 


cluster_mean_epi <- data.frame(data01_epi_subset[,union(subtype_markers,functional_markers)], clusterepi = data01_epi_subset$clusterepi, check.names = FALSE) %>%
  group_by(clusterepi) %>% summarize_all(list(mean))
cluster_mean_epi_mat<-as.matrix(cluster_mean_epi[,union(subtype_markers,functional_markers)])
rownames(cluster_mean_epi_mat)<-1:nrow(cluster_mean_epi_mat)
cluster_scaled_epi<-t(scale(t(cluster_mean_epi_mat)))
rownames(cluster_scaled_epi)<-1:nrow(cluster_scaled_epi)


##add in epi annotations
data_full$clustercombined2 <- data_full$clustercombined
data_full[rownames(data01_epi_subset),]$clustercombined2 <- data01_epi_subset$cluster1mEpi

clusterlevelscomb<-union(clusterlevels, clusterlevelsT)
clusterlevelscomb<-union(clusterlevelscomb, clusterlevelsEpi)
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusterlevelscomb))
clusternames<-clusterlevelscomb
names(colorassigned)<-clusternames
cluster_mean_merged <- data.frame(data01, cluster = data_full$clustercombined2, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_merged_mat<-as.matrix(cluster_mean_merged[,union(subtype_markers,functional_markers)])
cluster_scaled_merged<-t(scale(t(cluster_mean_merged_mat)))
rownames(cluster_scaled_merged)<-cluster_mean_merged$cluster

color_list = list(clusters=colorassigned)

cp<-rowAnnotation(#clusters=clusternames,
  col=color_list,
  gp = gpar(col = "white", lwd = .5),
  props= anno_barplot(
    as.vector(table(data_full$clustercombined2)/sum(table(data_full$clustercombined2))*100),
    gp = gpar(fill=colorassigned),
    border = F,
    bar_width = 0.75, 
    width = unit(2,"cm")))


canonmarkers<-c("ADAM10","KI67","CDX2","PanCK","HepPar","VIM","SMA","CD20","CD68","CD3","CD4","CD8","FOXP3","DCLAMP","CXCR3","CD15","CD163","CD57")

pdf("WH_clusterheatmap_combined2+.pdf",width=10,height=6)
Heatmap(cluster_scaled_merged[,canonmarkers],
        column_title="Phenograph Merged Clusters",
        name = "scaled",
        col=rev(brewer.rdbu(100)),
        cluster_columns = T,
        cluster_rows = T,
        border = F,
        #rect_gp = gpar(col = "white", lwd = .5),
        right_annotation = cp,
        show_row_names = T,
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=10),
        heatmap_legend_param = list(at=seq(from = round(min(cluster_scaled)), to = round(max(cluster_scaled)))),
        width = ncol(cluster_scaled_merged[,canonmarkers])*unit(4, "mm"), 
        height = nrow(cluster_scaled_merged[,canonmarkers])*unit(4, "mm"))
dev.off()

####resave====

backup_output <- list(data_full, data, data01, csv_full)
saveRDS(backup_output, "backup_output_0.RDS")

#data_full has all annotations
#data is expression profile only
#data01 is 0 -- 1 percentile normazliation

#in data_full
#column clustercombined2 has both detailed T and Epi clustering
#column clustercombined has detailed T clustering only
#column cluster1m just has broad clustering only

