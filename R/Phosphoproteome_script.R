#### Loading libraries ####
library("tidyverse")        
library("ggpubr")
library("PhosR")
library("limma")
library("data.table")
library("SummarizedExperiment")
library("prodlim")
library("clusterProfiler")
library("naniar")
library("GGally")
library("Biobase")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("ggnewscale")
library("ggridges")
library("enrichplot")
library("ggrepel")
library("directPA")
library("calibrate")
library("reactome.db")
library("annotate")
library("KSEAapp")
library("gplots")
library("pheatmap")
library("ggseqlogo")

#START
setwd("C:/Users/mct330/Desktop/Work/TINX_STUDY/TINX_R_analyses/")

#### Importing and preprocessing data ####

#loading stable phospho sites
Human_SPSs <- read.csv("Human_SPSs.txt", sep = "\t", header = F, row.names = NULL, dec = ".", stringsAsFactors = F)

#loading data ready for normalisation.
ppe <- read.csv("Combined_ins_ex.txt", sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)

#### converting to numeric values
ppe$Localization.prob <- as.numeric(gsub(",", ".", gsub("\\.", "", ppe$Localization.prob)))

# delete reverse matches and potential contaminants
del <- which(ppe[,"Reverse"]=="+" | ppe[,"Potential.contaminant"]=="+")
ppe <- ppe[-del,]

#Only keeping class 1 sites
class_1_sites <- ppe$Localization.prob >= 0.75
ppe <- ppe[class_1_sites,]

#Rename columns after experiment name"
names(ppe)[1:10] <- c("Basal_1", "Insulin_1","Basal_2", "Insulin_2","Basal_3", "Insulin_3","Basal_4", "Insulin_4","Basal_5", "Insulin_5")
names(ppe)[11:20] <- c("Rest_1", "Exercise_1","Rest_2", "Exercise_2","Rest_3", "Exercise_3","Rest_4", "Exercise_4","Rest_5", "Exercise_5")

Gene_index <- ppe[,c("Gene.names","Amino.acid","Position")]
Gene_index <- paste(Gene_index$Gene.names,";",Gene_index$Amino.acid,Gene_index$Position, ";")
Gene_index <- gsub("[[:space:]]", "", Gene_index) 
ppe$Gene_index <- Gene_index

#### Log2 transform, removing missing values, medianScaling, ratio calculation ####
ppe[, 1:20][ppe[, 1:20] == 0] <- NA
ppe[,1:20] = as.matrix(log2(ppe[,1:20]))

total2_grp = c(rep(1,10),rep(2,10))

com_tot <- as.matrix(ppe[,1:20])
rownames(com_tot) <- ppe$Unique.identifier

filt_com_tot <- selectGrps(com_tot, total2_grp, 1, n=1)
UID_names <- rownames(filt_com_tot)

rownames(com_tot) <- Gene_index
filt_com_tot <- selectGrps(com_tot, total2_grp, 1, n=1)

filt_com_tot <- medianScaling(filt_com_tot)

new_ratio <- cbind(filt_com_tot[,c(2,4,6,8,10)] - filt_com_tot[,c(1,3,5,7,9)], filt_com_tot[,c(12,14,16,18,20)] - filt_com_tot[,c(11,13,15,17,19)])

com_grp <- c("ins","ins","ins","ins","ins","Ex","Ex","Ex","Ex","Ex")

#### Batch correction/normalisation to phosphosites ####
set.seed(123)
design = model.matrix(~ com_grp - 1)
ctl = which(rownames(new_ratio) %in% Human_SPSs[,1])
ppe_final_ins = RUVphospho(new_ratio, M=design, k = 2, ctl = ctl)

new_proc_ins <- data.frame(cbind(filt_com_tot[,c(1,3,5,7,9)],filt_com_tot[,c(1,3,5,7,9)]+ppe_final_ins[,1:5]))
new_proc_ex <- data.frame(cbind(filt_com_tot[,c(11,13,15,17,19)],filt_com_tot[,c(11,13,15,17,19)]+ppe_final_ins[,6:10]))
colnames(new_proc_ins) <- c("Basal_1","Basal_2","Basal_3","Basal_4","Basal_5","Insulin_1","Insulin_2","Insulin_3","Insulin_4","Insulin_5")
colnames(new_proc_ex) <- c("Rest","Rest","Rest","Rest","Rest","Exercise","Exercise","Exercise","Exercise","Exercise")

rownames(new_proc_ins) <- UID_names
rownames(new_proc_ex) <- UID_names

#new_proc_ins <- new_proc_ins[,1:10]
Data_setup_ins <- matrix(c("Basal","Basal","Basal","Basal","Basal","Insulin","Insulin","Insulin","Insulin","Insulin","1","2","3","4","5","1","2","3","4","5"), ncol = 10, byrow = TRUE) %>% t() %>% 
  'colnames<-' (c("stimulus","Paired")) %>% data.frame()

Data_setup_ex <- matrix(c("Rest","Rest","Rest","Rest","Rest","Exercise","Exercise","Exercise","Exercise","Exercise","1","2","3","4","5","1","2","3","4","5"), ncol = 10, byrow = TRUE) %>% t() %>% 
  'colnames<-' (c("stimulus","Paired")) %>% data.frame()

new_proc_ins <- selectGrps(new_proc_ins, Data_setup_ins$stimulus, 1, n=2)
new_proc_ex <- selectGrps(new_proc_ex, Data_setup_ex$stimulus, 1, n=2)


#### PCA Insulin ####
t_exprs <- as.data.frame(t(new_proc_ins))
t_exprs <- cbind(t_exprs, Data_setup_ins$stimulus)
colnames(t_exprs)[9473] <- "Group"
pca_res <- prcomp(t_exprs[,1:9472], scale. = TRUE)

new_data_frame <- data.frame(pca_res$x)
new_data_frame$Group <- factor(Data_setup_ins$stimulus)
new_data_frame$Individual <- factor(Data_setup_ins$Paired)

theme<-theme(axis.title = element_text(size = 20), axis.line = element_line(colour = "black")
             ,panel.background = element_blank(),
             panel.border=element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),
             axis.ticks=element_line(colour="black"), legend.title = element_text(size=14),
             legend.text = element_text(size = 12))


percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(new_data_frame), "(", paste( as.character(percentage), "%", ")", sep="") )

ggplot(new_data_frame, aes(x=PC1,y=PC2, colour = Group, label=Individual)) + geom_point(aes(),size=3) + 
  theme + labs(color='Group') + 
  scale_color_manual(values=c("grey", "#b54a02")) + geom_text_repel(size=5) +
  geom_line(aes(group = Individual), colour="#b54a02", size=0.1) +
  xlab(percentage[1]) + ylab(percentage[2])

#### PCA Exercise ####
t_exprs <- as.data.frame(t(new_proc_ex))
t_exprs <- cbind(t_exprs, Data_setup_ex$stimulus)
colnames(t_exprs)[9318] <- "Group"
pca_res <- prcomp(t_exprs[,1:9317], scale. = TRUE)

new_data_frame <- data.frame(pca_res$x)
new_data_frame$Group <- factor(Data_setup_ex$stimulus)
new_data_frame$Individual <- factor(Data_setup_ex$Paired)

theme<-theme(axis.title = element_text(size = 20), axis.line = element_line(colour = "black")
             ,panel.background = element_blank(),
             panel.border=element_blank(), panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),
             axis.ticks=element_line(colour="black"), legend.title = element_text(size=14),
             legend.text = element_text(size = 12))

percentage <- round(factoextra::get_eig(pca_res)$variance.percent, 2)
percentage <- paste(colnames(new_data_frame), "(", paste( as.character(percentage), "%", ")", sep="") )

ggplot(new_data_frame, aes(x=PC1,y=PC2, colour = Group, label=Individual)) + geom_point(aes(),size=3) + 
  theme + labs(color='Group') + 
  scale_color_manual(values=c("#4444ab", "grey")) + geom_text_repel(size=5) +
  geom_line(aes(group = Individual), colour="#4444ab", size=0.1) +
  xlab(percentage[1]) + ylab(percentage[2])


#### insulin data ####
Data_setup_ins <- matrix(c("Basal","Basal","Basal","Basal","Basal","Insulin","Insulin","Insulin","Insulin","Insulin","1","2","3","4","5","1","2","3","4","5"), ncol = 10, byrow = TRUE) %>% t() %>% 
  'colnames<-' (c("stimulus","Paired")) %>% data.frame()

Paired <- factor(Data_setup_ins$Paired)
Treat<- factor(Data_setup_ins$stimulus, levels = c("Basal","Insulin"))
design_ins <- model.matrix(~ Treat -1) #Creating a design matrix
cm <- makeContrasts(TreatInsulin - TreatBasal, levels=design_ins)
corfit <- duplicateCorrelation(new_proc_ins, design_ins, block=Paired)
Insulin_table <- lmFit(new_proc_ins,design_ins,block=Paired,correlation=corfit$consensus) %>% contrasts.fit(., cm) %>% eBayes() %>% topTable(., 
                                                                                                                    coef="TreatInsulin - TreatBasal", 
                                                                                                                    number = Inf,
                                                                                                                    sort.by = "none")
Insulin_table$UID <- rownames(Insulin_table)
length(which(Insulin_table$adj.P.Val < 0.05)) #Number of differentially regulated sites
Insulin_table_final <- cbind(Insulin_table,ppe[which(ppe$Unique.identifier %in% rownames(Insulin_table)),c("Leading.proteins","Protein.names","Gene.names","Sequence.window","Amino.acid","Position","Gene_index","Unique.identifier","Multiplicity")])
total_Insulin_table_final <- cbind(Insulin_table_final,new_proc_ins)

#volcano plot
theme_vol <-theme(axis.title = element_text(size = 16), axis.line = element_line(colour = "black", size = 2)
                  ,panel.background = element_blank(),
                  panel.border= element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),strip.background=element_blank(),
                  axis.text.x=element_text(colour="black", size = 16),
                  axis.text.y=element_text(colour="black", size = 16),
                  axis.ticks=element_line(colour="black"), legend.position="none")

Insulin_table_final$diffexpressed <- "NO"
Insulin_table_final$diffexpressed[Insulin_table_final$logFC > 0 & Insulin_table_final$adj.P.Val < 0.05] <- "UP"
Insulin_table_final$diffexpressed[Insulin_table_final$logFC < 0 & Insulin_table_final$adj.P.Val < 0.05] <- "DOWN"
Insulin_table_final$delabel <- NA
Insulin_table_final$delabel[Insulin_table_final$diffexpressed != "NO"] <- Insulin_table_final$Gene_index[Insulin_table_final$diffexpressed != "NO"]
plot_Insulin_table_final <- ggplot() + geom_point(data=Insulin_table_final[Insulin_table_final$adj.P.Val <= 0.05,], aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel),alpha = 0.5, size = 1.5, stroke=NA)+ xlim(-2, 4) + scale_color_manual(values=c("grey43","#b54a02")) +
  geom_point(data=Insulin_table_final[Insulin_table_final$adj.P.Val > 0.05,], aes(x=logFC, y=-log10(P.Value), label = delabel), col="lightgrey",alpha=0.2,size=1, stroke=NA) + geom_text_repel() + ggtitle("Insulin vs Basal") + theme_vol
plot_Insulin_table_final + theme(aspect.ratio=1.25/1)

#Sites up/and downregulated
table(Insulin_table_final$logFC > 0 & Insulin_table_final$adj.P.Val < 0.05)["TRUE"]
table(Insulin_table_final$logFC < 0 & Insulin_table_final$adj.P.Val < 0.05)["TRUE"]

#Exporting for supplementary table
write.table(total_Insulin_table_final,"Insulin_table.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")

#### exercise data ####
Data_setup_ex <- matrix(c("Rest","Rest","Rest","Rest","Rest","Exercise","Exercise","Exercise","Exercise","Exercise","1","2","3","4","5","1","2","3","4","5"), ncol = 10, byrow = TRUE) %>% t() %>% 
  'colnames<-' (c("stimulus","Paired")) %>% data.frame()
Paired <- factor(Data_setup_ex$Paired)
Treat<- factor(Data_setup_ex$stimulus, levels = c("Rest","Exercise"))
design_ex <- model.matrix(~ Treat -1)
cm <- makeContrasts(TreatExercise - TreatRest, levels=design_ex)
corfit <- duplicateCorrelation(new_proc_ex, design_ex, block=Paired)
Exercise_table <- lmFit(new_proc_ex,design_ex,block=Paired,correlation=corfit$consensus) %>% contrasts.fit(., cm) %>% eBayes() %>% topTable(., 
                                                                                                                                           coef="TreatExercise - TreatRest", 
                                                                                                                                           number = Inf,
                                                                                                                                           sort.by = "none")
Exercise_table$UID <- rownames(Exercise_table)
length(which(Exercise_table$adj.P.Val < 0.05)) #Number of differentially regulated sites
Exercise_table_final <- cbind(Exercise_table,ppe[which(ppe$Unique.identifier %in% rownames(Exercise_table)),c("Leading.proteins","Protein.names","Gene.names","Sequence.window","Amino.acid","Position","Gene_index","Unique.identifier","Multiplicity")])
total_exercise_table_final <- cbind(Exercise_table_final,new_proc_ex)
colnames(total_exercise_table_final)[17:26] <- c("Rest_1","Rest_2","Rest_3","Rest_4","Rest_5","Exercise_1","Exercise_2","Exercise_3","Exercise_4","Exercise_5")

#volcano plot
Exercise_table_final$diffexpressed <- "NO"
Exercise_table_final$diffexpressed[Exercise_table_final$logFC > 0 & Exercise_table_final$adj.P.Val <= 0.05] <- "UP"
Exercise_table_final$diffexpressed[Exercise_table_final$logFC < 0 & Exercise_table_final$adj.P.Val <= 0.05] <- "DOWN"
Exercise_table_final$delabel <- NA
Exercise_table_final$delabel[Exercise_table_final$diffexpressed != "NO"] <- Exercise_table_final$Gene_index[Exercise_table_final$diffexpressed != "NO"]
plot_Exercise_table_final <- ggplot(data=Exercise_table_final, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + geom_point(alpha = 0.4, size = 3)+ xlim(-4, 4) + scale_color_manual(values=c("grey43", "lightgrey","#4444ab")) + ggtitle("Exercise vs Rest") + theme_vol + geom_text_repel()
plot_Exercise_table_final <- ggplot() + geom_point(data=Exercise_table_final[Exercise_table_final$adj.P.Val <= 0.05,], aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel),alpha = 0.5, size = 1.5, stroke=NA)+ xlim(-4, 4) + scale_color_manual(values=c("grey43","#4444ab")) +
  geom_point(data=Exercise_table_final[Exercise_table_final$adj.P.Val > 0.05,], aes(x=logFC, y=-log10(P.Value), label = delabel),col="lightgrey",alpha = 0.2, size = 1, stroke=NA)+ xlim(-4, 4) + ggtitle("Exercise vs Rest") + theme_vol + geom_text_repel()
plot_Exercise_table_final + theme(aspect.ratio = 1.25/1)

table(Exercise_table_final$logFC > 0 & Exercise_table_final$adj.P.Val < 0.05)["TRUE"]
table(Exercise_table_final$logFC < 0 & Exercise_table_final$adj.P.Val < 0.05)["TRUE"]

#Exporting for supplementary table
write.table(total_exercise_table_final,"Exercise_table.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")


#### combining two data sets ####
comb_test <- inner_join(Insulin_table_final,Exercise_table_final, by="UID", suffix = c(".Insulin",".Exercise"))
comb_test_sign <- comb_test[which(comb_test$adj.P.Val.Insulin < 0.05 & comb_test$adj.P.Val.Exercise < 0.05),]
comb_test_sign_UP_UP <- comb_test[which(comb_test$adj.P.Val.Insulin < 0.05 & comb_test$logFC.Insulin > 0 & comb_test$adj.P.Val.Exercise < 0.05 & comb_test$logFC.Exercise > 0),]
comb_test_sign_down_down <- comb_test[which(comb_test$adj.P.Val.Insulin < 0.05 & comb_test$logFC.Insulin < 0 & comb_test$adj.P.Val.Exercise < 0.05 & comb_test$logFC.Exercise < 0),]
comb_test_sign_UPins_downEx <- comb_test[which(comb_test$adj.P.Val.Insulin < 0.05 & comb_test$logFC.Insulin > 0 & comb_test$adj.P.Val.Exercise < 0.05 & comb_test$logFC.Exercise < 0),]
comb_test_sign_UPEx_downIns <- comb_test[which(comb_test$adj.P.Val.Insulin < 0.05 & comb_test$logFC.Insulin < 0 & comb_test$adj.P.Val.Exercise < 0.05 & comb_test$logFC.Exercise > 0),]

#Exporting for supplementary table
write.table(comb_test,"Combined_table.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")
write.table(comb_test_sign,"Combined_table_significant.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")

#### making heatmap ####
ONLY_significant <- comb_test[comb_test$adj.P.Val.Insulin < 0.05 & comb_test$adj.P.Val.Exercise < 0.05,]
ONLY_significant <- cbind(ONLY_significant,new_proc_ins[which(rownames(new_proc_ins) %in% ONLY_significant$UID),],new_proc_ex[which(rownames(new_proc_ex) %in% ONLY_significant$UID),])
colnames(ONLY_significant)[46:55] <- c("Rest_1","Rest_2","Rest_3","Rest_4","Rest_5","Exercise_1","Exercise_2","Exercise_3","Exercise_4","Exercise_5")
ONLY_significant <- cbind(ONLY_significant[,34:35],ONLY_significant[,c(41:45)]-ONLY_significant[,c(36:40)],ONLY_significant[,c(51:55)]-ONLY_significant[,c(46:50)])
ONLY_significant <- ONLY_significant[,c(28:49)]
colnames(ONLY_significant)[1:2] <- c("diffexpressed","delabel")

#plotting heatmap of sites regulated by both insulin/exercise
long_only_sign <- gather(ONLY_significant, ID, Intensity, -delabel, -diffexpressed)
long_only_sign$delabel <- word(long_only_sign$delabel,1,sep=";") 
heat_mat <- as.matrix(ONLY_significant[,3:12])
col<- colorRampPalette(c("darkgreen", "white", "purple"))(256)
heatmap.2(heat_mat, scale = "none", col = col, 
          trace = "none", density.info = "none")

#### Making delta ins/ex plots Figure 2 ####
#IDs to extract: UID9484 (SH3BP4 S246), UID802 (CLASP2 S370), UID24492 (SORBS1 S949/S953), UID13916 (AKT1S1 T246)
Ins_site <- total_Insulin_table_final[which(total_Insulin_table_final$UID == "UID6399"),c(17:26)]
Ex_site <- total_exercise_table_final[which(total_exercise_table_final$UID == "UID6399"),c(17:26)]

Ins_and_Ex_site <- cbind(Ins_site,Ex_site)
Long_Ins_and_Ex_site <- gather(Ins_and_Ex_site, ID,Intensity)

Long_Ins_and_Ex_site$subject <- word(Long_Ins_and_Ex_site$ID,2,sep="_") 
Long_Ins_and_Ex_site$condition <- word(Long_Ins_and_Ex_site$ID,1,sep="_") 
Long_Ins_and_Ex_site$type <- c(rep("Pre",5),rep("Post",5),rep("Pre",5),rep("Post",5))  

Long_Ins_and_Ex_site$study_day <- "Insulin"
Long_Ins_and_Ex_site$study_day[11:20] <- "Exercise"
Long_Ins_and_Ex_site$study_day <- factor(Long_Ins_and_Ex_site$study_day, levels=c("Insulin","Exercise"))

delta_values <- as.data.frame(cbind(Long_Ins_and_Ex_site[6:10,]$Intensity-Long_Ins_and_Ex_site[1:5,]$Intensity,
                                    Long_Ins_and_Ex_site[16:20,]$Intensity-Long_Ins_and_Ex_site[11:15,]$Intensity))
colnames(delta_values) <- c("Insulin","Exercise")

df_long <- gather(delta_values, key = "Type", value = "Value", Insulin, Exercise)

df_long$Type <- factor(df_long$Type, levels=c("Insulin","Exercise"))

# Plot
ggplot(df_long, aes(x = Type, y = Value, color=Type)) + 
  geom_point(size=3) + theme_bw() + scale_color_manual(values=c("darkred","darkblue")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 2))
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1.4/1) +
  ylab("delta (Post-Pre)") + xlab("") + ggtitle("CLASP1 S559") +
  geom_hline(yintercept = 0, linetype=2)

#### Preparing data for enrichment analysis ####
foreground_ins <- Insulin_table_final[which(Insulin_table_final$adj.P.Val < 0.05 & Insulin_table_final$logFC > 0),]$Sequence.window
background_ins <- Insulin_table_final[which(Insulin_table_final$adj.P.Val > 0.05),]$Sequence.window
foreground_ex <- Exercise_table_final[which(Exercise_table_final$adj.P.Val < 0.05 & Exercise_table_final$logFC > 0),]$Sequence.window
background_ex <- Exercise_table_final[which(Exercise_table_final$adj.P.Val > 0.05),]$Sequence.window

write.table(foreground_ins,"foreground_ins.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(background_ins,"background_ins.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(foreground_ex,"foreground_ex.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(background_ex,"background_ex.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

sig_test_ins <- subset(comb_test, adj.P.Val.Insulin < 0.05 & adj.P.Val.Exercise < 0.05 &
                         logFC.Insulin > 0 & logFC.Exercise > 0)
sig_test_ins <- sig_test_ins$Sequence.window.Insulin
back_sign <- comb_test$Sequence.window.Insulin

write.table(sig_test_ins,"foreground_ins_ex_up.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(back_sign,"background_ins_ex_up.txt", row.names = FALSE, col.names = FALSE, sep = "\t")


#### Insulin-activated kinases ####
KinLib_insulin <- read.delim("Final_Insulin_enrichment-analysis-result-table.txt",
                             sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)
up_KinLib_insulin <- KinLib_insulin[which(KinLib_insulin$enrichment_value_log2 > 0),]

converted_gene_names <- read.delim("Converted_gene_names.txt",sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)
colnames(converted_gene_names)[1] <- "kinase"
muscle_match_list <- inner_join(converted_gene_names,up_KinLib_insulin, by = "kinase")

protein_lib <- read.delim("Expression_atlas.tsv",
                          sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)

Found_muscle_match <- muscle_match_list[which(muscle_match_list$Approved.symbol %in% protein_lib$Gene),]

Found_muscle_match$BH_adj.p.val <- p.adjust(Found_muscle_match$p_value, method="BH")

ggplot(Found_muscle_match, aes(x=enrichment_value, y=p_value_log10_abs, label=kinase, color=enrichment_value)) + geom_point(size=2, alpha=0.4, stroke=NA) +
  geom_text_repel() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  theme + xlab("Enrichment score") + ylab("-log10 p.value") +
  geom_hline(yintercept=-log10(0.0025),linetype=2, size =0.25, alpha=0.4,color="black") +
  ggtitle("Insulin-activated kinases") + scale_colour_gradient(low = "#132B43",
                                                               high = "#b54a02")

write.table(Found_muscle_match,"Insulin_activated_kinases.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")

#### Exercise-activated kinases ####
KinLib_exercise <- read.delim("Final_Exercise_enrichment-analysis-result-table.txt",
                              sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)
up_KinLib_exercise <- KinLib_exercise[which(KinLib_exercise$enrichment_value_log2 > 0),]

muscle_match_list_exercise <- inner_join(converted_gene_names,up_KinLib_exercise, by = "kinase")

Found_muscle_match <- muscle_match_list_exercise[which(muscle_match_list_exercise$Approved.symbol %in% protein_lib$Gene),]

Found_muscle_match$BH_adjust_p.val <- p.adjust(Found_muscle_match$p_value, method="BH")

ggplot(Found_muscle_match, aes(x=enrichment_value, y=p_value_log10_abs, label=kinase, color=enrichment_value)) + geom_point(size=2, alpha=0.4, stroke=NA) +
  geom_text_repel() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  theme + xlab("Enrichment score") + ylab("-log10 p.value") +
  geom_hline(yintercept=-log10(0.021),linetype=2, size =0.25, alpha=0.4,color="black") +
  ggtitle("Exercise-activated kinases") + scale_colour_gradient(low = "#132B43",
                                                                high = "#4444ab")
write.table(Found_muscle_match,"Exercise_activated_kinases.txt", row.names = FALSE, col.names = TRUE, sep = "\t",dec=",")

#### Upregulated by insulin and exercise ####
KinLib_UPUP <- read.delim("Final_Exe_Ins_UP_enrichment-analysis-result-table.txt",
                          sep = "\t", header = T, row.names = NULL, dec = ".", stringsAsFactors = F)
up_KinLib_UPUP <- KinLib_UPUP[which(KinLib_UPUP$enrichment_value_log2 > 0),]

muscle_match_list <- inner_join(converted_gene_names,up_KinLib_UPUP, by = "kinase")

Found_muscle_match <- muscle_match_list[which(muscle_match_list$Approved.symbol %in% protein_lib$Gene),]

Found_muscle_match$BH_adj.p.val <- p.adjust(Found_muscle_match$p_value, method="BH")

Significant_IE_up <- Found_muscle_match[Found_muscle_match$BH_adj.p.val <= 0.05,]

Significant_IE_up$kinase <- as.factor(Significant_IE_up$kinase)

ggplot(Significant_IE_up,aes(x=p_value_log10_abs,y=Significant_IE_up$kinase)) + 
  geom_segment(aes(xend = 0, yend = Significant_IE_up$kinase), color = "grey50") +
  geom_point(aes(size=enrichment_value),color="#853c91",stroke=NA) + theme + xlab("-log10(p.value)") + ylab("") +
  scale_size(range = c(5, 10)) +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) +
  labs(size = "Enrichment score")

#### Preparing data for logo motif ####
new_seqs <- word(comb_test_sign_UP_UP$Sequence.window.Exercise,1,sep=";")
new_seqs <- str_sub(new_seqs, 9, -9)
new_seqs<- new_seqs[new_seqs != ""]

ggseqlogo( new_seqs, method="bits",seq_type='aa')

#### -3 R position ####
sum(str_count(str_sub(new_seqs, 5, 5)))
sum(str_count(str_sub(new_seqs, 5, 5), "R"))

new_seqs_background <- word(comb_test$Sequence.window.Insulin,1,sep=";")
new_seqs_background <- str_sub(new_seqs_background, 9, -9)
new_seqs_background <- new_seqs_background[new_seqs_background != ""]

matrix_data <- matrix(c(sum(str_count(str_sub(new_seqs, 5, 5), "R")), sum(str_count(str_sub(new_seqs, 5, 5))),
                        sum(str_count(str_sub(new_seqs_background, 5, 5), "R")), sum(str_count(str_sub(new_seqs_background, 5, 5)))), nrow=2)
result <- fisher.test(matrix_data)
print(result$p.value)

#### -5 R position ####
sum(str_count(str_sub(new_seqs, 3, 3)))
sum(str_count(str_sub(new_seqs, 3, 3), "R"))

new_seqs_background <- word(comb_test$Sequence.window.Insulin,1,sep=";")
new_seqs_background <- str_sub(new_seqs_background, 9, -9)
new_seqs_background <- new_seqs_background[new_seqs_background != ""]

sum(str_count(str_sub(new_seqs_background, 3, 3)))
sum(str_count(str_sub(new_seqs_background, 3, 3), "R"))

matrix_data <- matrix(c(sum(str_count(str_sub(new_seqs, 3, 3), "R")), sum(str_count(str_sub(new_seqs, 3, 3))),
                        sum(str_count(str_sub(new_seqs_background, 3, 3), "R")), sum(str_count(str_sub(new_seqs_background, 3, 3)))), nrow=2)
result <- fisher.test(matrix_data)
print(result$p.value)


#### -4 Q position ####
sum(str_count(str_sub(new_seqs, 3, 3)))
sum(str_count(str_sub(new_seqs, 3, 3), "R"))

new_seqs_background <- word(comb_test$Sequence.window.Insulin,1,sep=";")
new_seqs_background <- str_sub(new_seqs_background, 9, -9)
new_seqs_background <- new_seqs_background[new_seqs_background != ""]

sum(str_count(str_sub(new_seqs_background, 3, 3)))
sum(str_count(str_sub(new_seqs_background, 3, 3), "R"))

matrix_data <- matrix(c(sum(str_count(str_sub(new_seqs, 4, 4), "Q")), sum(str_count(str_sub(new_seqs, 4, 4))),
                        sum(str_count(str_sub(new_seqs_background, 4, 4), "Q")), sum(str_count(str_sub(new_seqs_background, 4, 4)))), nrow=2)
result <- fisher.test(matrix_data)
print(result$p.value)