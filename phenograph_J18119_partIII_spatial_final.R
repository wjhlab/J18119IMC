####Part 3 Distance Relationships####
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
cbbPalette <- c("#E69F00", "#999999", "#56B4E9", "#FF99FF","#F0E442", "#0072B2", "#999978", "#009E73", "#CC0000")
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

#how many of each cell types are there in the dataset?

counts_table <- table(data_full$clustercombined2, data_full$sample_id)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdft <- melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
              id.vars = "cluster", value.name = "counts", 
              variable.name = "sample_id")
ggdft$sample_id <- factor(ggdft$sample_id, levels=samplevels)
ggdft$cluster <- factor(ggdft$cluster, levels=clusterlevels)
ggdft$Timepoint <- factor(md$Timepoint[match(ggdft$sample_id,md$Sample_ID)], levels=Timepointlevels)
ggdft$Response <- factor(md$Response[match(ggdft$sample_id,md$Sample_ID)], levels=Responselevels)
ggdftPreonly <- ggdft[ggdft$Timepoint=="PreTx",]

totalcounts<-ggdft %>% group_by(cluster, Timepoint, Response) %>% summarize_at(vars(counts),funs(sum))
write.csv(totalcounts,"Totalcounts.csv")

totalcountspreonly <- ggdftPreonly %>% group_by(cluster, Timepoint, Response) %>% summarize_at(vars(counts),funs(sum))
write.csv(totalcountspreonly,"Totalcountspreonly.csv")

totalcounts_Response_Pre <- ggdftPreonly %>% group_by(cluster, Response) %>% summarize_at(vars(counts),funs(sum))
write.csv(totalcounts_Response_Pre,"Totalcounts_Response_Pre.csv")

#percentage of each respective total

totalcounts_Response_NR_Pre <- totalcounts_Response_Pre[totalcounts_Response_Pre$Response=="NR",]
totalcounts_Response_R_Pre <- totalcounts_Response_Pre[totalcounts_Response_Pre$Response=="R",]

totalResponse_NR_Pre<-sum(totalcounts_Response_NR_Pre$counts)
totalResponse_R_Pre<-sum(totalcounts_Response_R_Pre$counts)

pct_Response_NR_Pre <- totalcounts_Response_NR_Pre$counts/totalResponse_NR_Pre*100
names(pct_Response_NR_Pre)<- totalcounts_Response_NR_Pre$cluster
pct_Response_R_Pre <- totalcounts_Response_R_Pre$counts/totalResponse_R_Pre*100
names(pct_Response_R_Pre)<- totalcounts_Response_R_Pre$cluster

####EXTRACTING MEAN EXPRESSIONS OF KEY MARKERS AND COLORLEGENDS BASED ON EXPRESSION####
#list out the cell types and create a legends data frame
allcelltypes<-clusterlevels
legendctype<-as.data.frame(cbind(paste0("ctype",1:length(allcelltypes)),allcelltypes))
legendctype$maintype<-1
legendctype$maintype[str_detect(legendctype$allcelltypes,"B")]<-"B"
legendctype$maintype[str_detect(legendctype$allcelltypes,"DC")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"DPT")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Epi_I")]<-"Epi"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Epi_II")]<-"Epi"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Epi_III")]<-"Epi"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Epi_KI67")]<-"Epi"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Epi_pSTAT3")]<-"Epi"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Hep")]<-"Hep"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mac_I")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mac_II")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"NKT")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Neutrophil")]<-"Neutrophil"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Stroma")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tc")]<-"Tc_Th"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tc1")]<-"Tc_Th"
legendctype$maintype[str_detect(legendctype$allcelltypes,"TcM")]<-"Tc_Th"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Th")]<-"Tc_Th"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Th1")]<-"Tc_Th"
legendctype$maintype[str_detect(legendctype$allcelltypes,"ThM")]<-"Tc_Th"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Treg")]<-"Treg"
legendctype$maintype[str_detect(legendctype$allcelltypes,"UA")]<-"Other"

####SETUP DATAFRAME####

#identify which cell types in either of the data subsets are very rare - have less than 0.1% cells out of each subset
#also exclude Hep and UA clusters = ctype1, ctype22

exclude_Response_NR_Pre<-which(pct_Response_NR_Pre<0.1)
exclude_Response_NR_Pre<-legendctype$V1[match(names(exclude_Response_NR_Pre), legendctype$allcelltypes)]
exclude_Response_NR_Pre<-c(exclude_Response_NR_Pre, "ctype1")
exclude_Response_NR_Pre<-c(exclude_Response_NR_Pre, "ctype22")
exclude_Response_R_Pre<-which(pct_Response_R_Pre<0.1)
exclude_Response_R_Pre<-legendctype$V1[match(names(exclude_Response_R_Pre), legendctype$allcelltypes)]
exclude_Response_R_Pre<-union(exclude_Response_R_Pre, "ctype1")
exclude_Response_R_Pre<-union(exclude_Response_R_Pre, "ctype22")

#get out expression levels, X, and Y coords

#sort panels into different categories
cleanpanel <- read_xlsx(panelDataFile)
subtype_markers <- cleanpanel$clean_names[cleanpanel$subtype == 1]
functional_markers <- cleanpanel$clean_names[cleanpanel$functional == 1]
otherparameters <- cleanpanel$clean_names[cleanpanel$other ==1]
cluster_by <- cleanpanel$clean_names[cleanpanel$cluster_by == 1]

expr0<-data.frame(data_full[,c(union(subtype_markers,functional_markers),"CellId","X_coord","Y_coord")],
                  cluster=data_full$clustercombined2,
                  sample_id=data_full$sample_id)

expr0$Timepoint <- factor(md$Timepoint[match(expr0$sample_id,md$Sample_ID)], levels=Timepointlevels)
expr0$Response <- factor(md$Response[match(expr0$sample_id,md$Sample_ID)], levels=Responselevels)

expr0<-expr0[expr0$cluster!="UA",]

expr0_Timepoint_pre <- expr0[expr0$Timepoint=="PreTx",]
expr0_Response_NR_Pre <- expr0_Timepoint_pre[expr0_Timepoint_pre$Response=="NR",]
expr0_Response_R_Pre <- expr0_Timepoint_pre[expr0_Timepoint_pre$Response=="R",]

####CREATE DISTANCE MATRICES FOR PROGRESSION CRITERIA####

Responseid_NR_Pre<-unique(expr0_Response_NR_Pre$sample_id)
Responseid_R_Pre<-unique(expr0_Response_R_Pre$sample_id)

##Cells for Non-Responders before Treatment
expr_Response_NR_Pre<-c()

for(k in 1:length(Responseid_NR_Pre)){
  expr_k<-expr0_Response_NR_Pre[expr0_Response_NR_Pre$sample_id==Responseid_NR_Pre[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-legendctype$allcelltypes
  dummy<-as.data.frame(dummy)
  expr_k<-data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  expr_Response_NR_Pre <- rbind(expr_Response_NR_Pre, expr_k)
}
colnames(expr_Response_NR_Pre)<-c(colnames(expr0_Response_NR_Pre),legendctype$V1)
expr_Response_NR_Pre_m<-as.matrix(expr_Response_NR_Pre[,colnames(expr_Response_NR_Pre)[str_detect(colnames(expr_Response_NR_Pre),"ctype")]])
expr_Response_NR_Pre$ctype_no<-paste0("ctype",match(expr_Response_NR_Pre$cluster,legendctype$allcelltypes))
rownames(expr_Response_NR_Pre_m)<-expr_Response_NR_Pre$ctype_no
expr_Response_NR_Pre_m[is.infinite(expr_Response_NR_Pre_m)]<-NA

expr_Response_NR_Pre_combined <- cbind(expr0_Response_NR_Pre,expr_Response_NR_Pre)

expr_Response_NR_Pre<-expr_Response_NR_Pre_m
expr_Response_NR_Pre[expr_Response_NR_Pre == 0] <- NA

saveRDS(expr_Response_NR_Pre,'backup_dist_expr_Response_NR_Pre.rds')

##Cells for Responders before Treatment
expr_Response_R_Pre<-c()

for(k in 1:length(Responseid_R_Pre)){
  expr_k<-expr0_Response_R_Pre[expr0_Response_R_Pre$sample_id==Responseid_R_Pre[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-legendctype$allcelltypes
  dummy<-as.data.frame(dummy)
  expr_k<-data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  expr_Response_R_Pre <- rbind(expr_Response_R_Pre, expr_k)
}
colnames(expr_Response_R_Pre)<-c(colnames(expr0_Response_R_Pre),legendctype$V1)
expr_Response_R_Pre_m<-as.matrix(expr_Response_R_Pre[,colnames(expr_Response_R_Pre)[str_detect(colnames(expr_Response_R_Pre),"ctype")]])
expr_Response_R_Pre$ctype_no<-paste0("ctype",match(expr_Response_R_Pre$cluster,legendctype$allcelltypes))
rownames(expr_Response_R_Pre_m)<-expr_Response_R_Pre$ctype_no
expr_Response_R_Pre_m[is.infinite(expr_Response_R_Pre_m)]<-NA

expr_Response_R_Pre_combined <- cbind(expr0_Response_R_Pre,expr_Response_R_Pre)

expr_Response_R_Pre<-expr_Response_R_Pre_m
expr_Response_R_Pre[expr_Response_R_Pre == 0] <- NA

saveRDS(expr_Response_R_Pre,'backup_dist_expr_Response_R_Pre.rds')

####SKIP TO HERE IF DIST MATRICES ARE ALREADY SAVED####
#load previously saved distance matrix

expr_Response_NR_Pre<- readRDS('backup_dist_expr_Response_NR_Pre.rds')
expr_Response_R_Pre<- readRDS('backup_dist_expr_Response_R_Pre.rds')

#create expression + distance data frames

expr_Response_NR_Pre_combined <- cbind(expr0_Response_NR_Pre,expr_Response_NR_Pre)
expr_Response_R_Pre_combined <- cbind(expr0_Response_R_Pre,expr_Response_R_Pre)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_Response_NR_Prerm0<-expr_Response_NR_Pre[,colnames(expr_Response_NR_Pre) %nin% colnames(expr_Response_NR_Pre)[colSums(expr_Response_NR_Pre, na.rm = T)==0]]
expr_Response_R_Prerm0<-expr_Response_R_Pre[,colnames(expr_Response_R_Pre) %nin% colnames(expr_Response_R_Pre)[colSums(expr_Response_R_Pre, na.rm = T)==0]]

#remove rows/columns where there are very rare Cell types

expr_Response_NR_Prerm<-expr_Response_NR_Prerm0[rownames(expr_Response_NR_Prerm0) %nin% exclude_Response_NR_Pre, colnames(expr_Response_NR_Prerm0) %nin% exclude_Response_NR_Pre]
expr_Response_R_Prerm<-expr_Response_R_Prerm0[rownames(expr_Response_R_Prerm0) %nin% exclude_Response_R_Pre, colnames(expr_Response_R_Prerm0) %nin% exclude_Response_R_Pre]

####NETWORK VISUALIZATION####  

#for Cells of Non-Responders before Treatment
mat_Response_NR_Pre=aggregate(x=expr_Response_NR_Prerm, by=list(rownames(expr_Response_NR_Prerm)), FUN=mean, na.rm=T)
groupnames<-mat_Response_NR_Pre$Group.1
mat_Response_NR_Pre<-as.matrix(mat_Response_NR_Pre[,2:ncol((mat_Response_NR_Pre))])
rownames(mat_Response_NR_Pre)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_Response_NR_Pre)<-str_sub(colnames(mat_Response_NR_Pre),6,) #simplify colnames
dist_Response_NR_Pre<-mat_Response_NR_Pre[colnames(mat_Response_NR_Pre),]

#for Cells of Responders before Treatment
mat_Response_R_Pre=aggregate(x=expr_Response_R_Prerm, by=list(rownames(expr_Response_R_Prerm)), FUN=mean, na.rm=T)
groupnames<-mat_Response_R_Pre$Group.1
mat_Response_R_Pre<-as.matrix(mat_Response_R_Pre[,2:ncol((mat_Response_R_Pre))])
rownames(mat_Response_R_Pre)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_Response_R_Pre)<-str_sub(colnames(mat_Response_R_Pre),6,) #simplify colnames
dist_Response_R_Pre<-mat_Response_R_Pre[colnames(mat_Response_R_Pre),]

## Manual Color Editing

#Color Editing
custom_colors <- c("ctype1" = "#56B4E9",
                   "ctype2" = "#999999",
                   "ctype3" = "#999999",
                   "ctype4" = "#999999",
                   "ctype5" = "#999999",
                   "ctype6" = "#999999",
                   "ctype7" = "#009E73",
                   "ctype8" = "#009E73",
                   "ctype9" = "#009E73",
                   "ctype10" = "#009E73",
                   "ctype11" = "#009E73",
                   "ctype12" = "#009E73",
                   "ctype13" = "#CC0000",
                   "ctype14" = "#FF99FF",
                   "ctype15" = "#FF99FF",
                   "ctype16" = "#E69F00",
                   "ctype17" = "#F0E442",
                   "ctype18" = "#F0E442",
                   "ctype19" = "#F0E442",
                   "ctype20" = "#0072B2",
                   "ctype21" = "#999978",
                   "ctype22" = "plum")

legendctypeex<-legendctype[legendctype$maintype!="Other",] #excluding Other subtypes

#Make sure cols matches your cell type colors
#Use following command to verify hex codes in cols match cell type order:
#levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"]

#for Cells of Non-Responders Before Treatment
pdf("Distance_Relationships_Response_NR_Pre.pdf",width=8,height=6)

xx<-dist_Response_NR_Pre
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlist_Response_NR_Pre<-custom_colors[as.numeric(colnames(mat_Response_NR_Pre))] #set color
names(colorlist_Response_NR_Pre)<-as.character(colnames(mat_Response_NR_Pre)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_Response_NR_Pre[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Cells of Non-Responders Before Treatment 1/yy (from each cluster node)",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_Response_NR_Pre,1/dist_Response_R_Pre, na.rm=T), #IMPORTANT! 
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_Response_NR_Pre[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% c("Hep","Other")], 
       col = cbbPalette , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)
dev.off()

#for Cells of Responders Before Treatment
pdf("Distance_Relationships_Response_R_Pre.pdf",width=8,height=6)

xx<-dist_Response_R_Pre
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlist_Response_R_Pre<-custom_colors[as.numeric(colnames(mat_Response_R_Pre))] #set color
names(colorlist_Response_R_Pre)<-as.character(colnames(mat_Response_R_Pre)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_Response_R_Pre[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Cells of Responders Before Treatment 1/yy (from each cluster node)",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_Response_R_Pre,1/dist_Response_NR_Pre, na.rm=T), #IMPORTANT! 
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_Response_R_Pre[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% c("Hep","Other")], 
       col = cbbPalette , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)
dev.off()

####VIOLIN PLOTS OF DISTANCES####

###TcM-Treg
df_Response_NR_Pre<-melt(expr_Response_NR_Pre)
df_Response_NR_Pre$Response <- "NR"
df_Response_R_Pre<-melt(expr_Response_R_Pre)
df_Response_R_Pre$Response <- "R"
df_combined_Pre<-rbind(df_Response_NR_Pre,df_Response_R_Pre)
colnames(df_combined_Pre)<-c("var1", "var2","Distance (µm)", "Response")
df_combined_Pre<-df_combined_Pre[!is.na(df_combined_Pre$`Distance (µm)`),]
df_combined_c1ind_Pre<-df_combined_Pre[df_combined_Pre$var1=="ctype9",]

pdf('Distance_violinplots_TcM_Pre.pdf',width=3,height=3)
ggplot(df_combined_c1ind_Pre[df_combined_c1ind_Pre$var2=="ctype13",], aes(x=Response, y=`Distance (µm)`, fill=Response))+
  ylim(0,500)+ggtitle('TcM_Treg')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))+theme(
    plot.title = element_text(hjust = 0.5))+stat_summary(fun = "mean",
                                                         geom = "crossbar", 
                                                         width = 0.4,
                                                         size = 0.1,
                                                         position = position_dodge(0.9),
                                                         colour = "black") +
  stat_compare_means(
    aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
    method = "t.test",
    label.x = 1.4,
    size = 2,         # Increase size of the p-value labels
    color = "black",  # Change color of the p-value labels
    hide.ns = TRUE)    # Hide non-significant results ("ns")

dev.off()

###ThM-Treg
df_Response_NR_Pre<-melt(expr_Response_NR_Pre)
df_Response_NR_Pre$Response <- "NR"
df_Response_R_Pre<-melt(expr_Response_R_Pre)
df_Response_R_Pre$Response <- "R"
df_combined_Pre<-rbind(df_Response_NR_Pre,df_Response_R_Pre)
colnames(df_combined_Pre)<-c("var1", "var2","Distance (µm)", "Response")
df_combined_Pre<-df_combined_Pre[!is.na(df_combined_Pre$`Distance (µm)`),]
df_combined_c1ind_Pre<-df_combined_Pre[df_combined_Pre$var1=="ctype12",]

pdf('Distance_violinplots_ThM_Pre.pdf',width=3,height=3)
ggplot(df_combined_c1ind_Pre[df_combined_c1ind_Pre$var2=="ctype13",], aes(x=Response, y=`Distance (µm)`, fill=Response))+
  ylim(0,500)+ggtitle('ThM_Treg')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))+theme(
    plot.title = element_text(hjust = 0.5))+stat_summary(fun = "mean",
                                                         geom = "crossbar", 
                                                         width = 0.4,
                                                         size = 0.1,
                                                         position = position_dodge(0.9),
                                                         colour = "black") +
  stat_compare_means(
    aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
    method = "t.test",
    label.x = 1.4,
    size = 2,         # Increase size of the p-value labels
    color = "black",  # Change color of the p-value labels
    hide.ns = TRUE)   # Hide non-significant results ("ns")

dev.off()

###Th-Treg
df_Response_NR_Pre<-melt(expr_Response_NR_Pre)
df_Response_NR_Pre$Response <- "NR"
df_Response_R_Pre<-melt(expr_Response_R_Pre)
df_Response_R_Pre$Response <- "R"
df_combined_Pre<-rbind(df_Response_NR_Pre,df_Response_R_Pre)
colnames(df_combined_Pre)<-c("var1", "var2","Distance (µm)", "Response")
df_combined_Pre<-df_combined_Pre[!is.na(df_combined_Pre$`Distance (µm)`),]
df_combined_c1ind_Pre<-df_combined_Pre[df_combined_Pre$var1=="ctype10",]

pdf('Distance_violinplots_Th_Pre.pdf',width=3,height=3)
ggplot(df_combined_c1ind_Pre[df_combined_c1ind_Pre$var2=="ctype13",], aes(x=Response, y=`Distance (µm)`, fill=Response))+
  ylim(0,500)+ggtitle('Th_Treg')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))+theme(
    plot.title = element_text(hjust = 0.5))+stat_summary(fun = "mean",
                                                         geom = "crossbar", 
                                                         width = 0.4,
                                                         size = 0.1,
                                                         position = position_dodge(0.9),
                                                         colour = "black") +
  stat_compare_means(
    aes(label = ifelse(..p.signif.. == "", "", paste0(..p.format.., " ", ..p.signif..))),
    method = "t.test",
    label.x = 1.4,
    size = 2,         # Increase size of the p-value labels
    color = "black",  # Change color of the p-value labels
    hide.ns = TRUE)    # Hide non-significant results ("ns")

dev.off()
