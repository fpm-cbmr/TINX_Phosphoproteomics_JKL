#### Libraries ####
library(tidyverse)         
library(limma)             
library(ggrepel)           
library(clusterProfiler)   
library(org.Hs.eg.db)      
library(PhosR)             
library(biomaRt)           
library(GOSemSim)          

#### Figure 3G - In-situ contractions +/- insulin Reps1 phos & GU ####
In_situ_data <- read.delim("In_situ_input_R.txt", dec=",")
In_situ_data <- In_situ_data[!c(In_situ_data$WB_log2 < 22),]

Rest_Ins_subset <- In_situ_data[In_situ_data$Contraction == "Rest" & In_situ_data$Insulin == "Insulin",]
shapiro.test(Rest_Ins_subset$WB_scaled)
cor.test(Rest_Ins_subset$GU,Rest_Ins_subset$WB_scaled, method="pearson")
(0.7918247)^2

ggplot(Rest_Ins_subset, aes(x=WB_scaled, y=GU)) + geom_point(aes(color = Contraction, shape=Insulin),size=5) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.25, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1, color="black"),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16),legend.position = "none", aspect.ratio=1/1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) + 
  scale_color_manual(values=c("darkred")) + scale_shape_manual(values = c(16, 21)) +
  ggtitle("Spearman correlation") +
  xlab("p-Reps1 S709 [AU]") + ylab("μmol/g/hr") +
  annotate(geom="text",x=1,y=40, label="R2 = 0.63", color = "black",size=5) +
  annotate(geom = "text",x=1,y=30,label="p = 0.006",size=5) + xlim(0.5,1.2)

#### Figure S3F - In-situ contractions +/- insulin delta correlations ####
In_situ_delta_data <- read.delim("In_situ_contraction_delta_correlations.txt",dec=",")

In_situ_delta_data_SALINE <- In_situ_delta_data[In_situ_delta_data$Condition == "Saline",]
cor.test(In_situ_delta_data_SALINE$delta_GU,In_situ_delta_data_SALINE$scaled_phos,method="pearson")
-(0.2535771)^2

ggplot(In_situ_delta_data_SALINE, aes(x=scaled_phos, y=delta_GU)) + geom_point(size=3,color="darkblue") +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.25, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1,color="black"),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  ggtitle("Glucose uptake") +
  xlab(" delta p-Reps1 S709 [AU]") + ylab("delta μmol/g/hr") +
  annotate(geom="text",x=1,y=550, label="R2 = -0.064", color = "black",size=5) +
  annotate(geom = "text",x=1,y=520,label="p = 0.58",size=5)

#### Figure 3E - Barplots GU data in-sito contractions +/- insulin ####
GU_data_in_situ <- read.delim("In_situ_GU_data.txt", dec=",")
GU_data_in_situ$ID <- factor(GU_data_in_situ$ID)

# Adjusting factor levels
GU_data_in_situ$Insulin <- factor(GU_data_in_situ$Insulin, levels = c("Saline", "Insulin"))
GU_data_in_situ$Contraction <- factor(GU_data_in_situ$Contraction, levels = c("Rest", "Contraction"))

mean_df <- GU_data_in_situ %>%
  group_by(Insulin, Contraction) %>%
  summarise(MeanIntensity = mean(GU, na.rm = TRUE), .groups = "drop")


# Split data into 'Rest' and 'Contraction'
rest_data <- GU_data_in_situ %>% filter(Contraction == "Rest")
contraction_data <- GU_data_in_situ %>% filter(Contraction == "Contraction")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(rest_data, contraction_data, by = c("ID", "Insulin"), suffix = c("_rest", "_contraction"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("Rest" = "darkred", "Contraction" = "darkblue")  # Colors for points


# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Contraction, y = MeanIntensity, fill = Contraction, color=Contraction), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = GU_data_in_situ, aes(x = Contraction, y = GU, color = Contraction), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("Rest"), xend = factor("Contraction"), y = GU_rest, yend = GU_contraction, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Insulin) +
  
  # Customizations
  labs(x = "", y = "Glucose Uptake (GU)", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 3/1)

print(p)

#### Figure 3F - Reps1 in-situ data ####
REPS1_data_in_situ <- read.delim("Reps1_In_situ_data.txt", dec=",")
REPS1_data_in_situ$ID <- factor(REPS1_data_in_situ$ID)

# Adjusting factor levels
REPS1_data_in_situ$Insulin <- factor(REPS1_data_in_situ$Insulin, levels = c("Saline", "Insulin"))
REPS1_data_in_situ$Contraction <- factor(REPS1_data_in_situ$Contraction, levels = c("Rest", "Contraction"))


# Split data into 'Rest' and 'Contraction'
rest_data <- REPS1_data_in_situ %>% filter(Contraction == "Rest")
contraction_data <- REPS1_data_in_situ %>% filter(Contraction == "Contraction")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(rest_data, contraction_data, by = c("ID", "Insulin"), suffix = c("_rest", "_contraction"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("Rest" = "darkred", "Contraction" = "darkblue")  # Colors for points


mean_df <- REPS1_data_in_situ %>%
  group_by(Insulin, Contraction) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Contraction, y = MeanIntensity, fill = Contraction, color=Contraction), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = REPS1_data_in_situ, aes(x = Contraction, y = Intensity, color = Contraction), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("Rest"), xend = factor("Contraction"), y = Intensity_rest, yend = Intensity_contraction, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Insulin) +
  
  # Customizations
  labs(x = "", y = "Glucose Uptake (GU)", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 3/1)

print(p)

#### Figure 3B - Human WB validation of p-REPS1 S709 Basal/Insulin ####
Human_ins_phos_REPS1 <- read.delim("pREPS1_S709_WB_Insulin.txt",dec = ",")

Human_ins_phos_REPS1$ID <- as.factor(Human_ins_phos_REPS1$ID)
Human_ins_phos_REPS1$Clamp <- factor(Human_ins_phos_REPS1$Clamp, levels=c("Pre","Post"))

mean_df <- Human_ins_phos_REPS1 %>%
  group_by(Clamp) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

ggplot(data=Human_ins_phos_REPS1, aes(x=Clamp,y=Intensity)) +
  geom_bar(data = mean_df, aes(y = MeanIntensity, fill = Clamp, color=Clamp),width=0.60, 
           stat = "identity", position = position_dodge(width = 0.8)) + 
  geom_point(aes(color=Clamp), size=5) +
  geom_line(aes(group=ID),alpha = 0.6, size=0.3) +
  labs(x = "", y = "Fold change from Basal", title = "p-REPS1 S709") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = c("black","darkred")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 1.5/1) +
  scale_x_discrete(labels=c('Basal', 'Insulin'))


#### Figure 3C - Human WB validation of p-REPS1 S709 Rest/Exercise ####
Human_Ex_phos_REPS1 <- read.delim("pREPS1_S709_WB_Exercise.txt",dec = ",")

Human_Ex_phos_REPS1$ID <- as.factor(Human_Ex_phos_REPS1$ID)
Human_Ex_phos_REPS1$Exercise <- factor(Human_Ex_phos_REPS1$Exercise, levels=c("Pre","Post"))

mean_df <- Human_Ex_phos_REPS1 %>%
  group_by(Exercise) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

ggplot(data=Human_Ex_phos_REPS1, aes(x=Exercise,y=Intensity)) +
  geom_bar(data = mean_df, aes(y = MeanIntensity, fill = Exercise, color=Exercise),width=0.60, 
           stat = "identity", position = position_dodge(width = 0.8)) + 
  geom_point(aes(color=Exercise), size=5) +
  geom_line(aes(group=ID),alpha = 0.6, size=0.3) +
  labs(x = "", y = "Fold change from Rest", title = "p-REPS1 S709") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = c("black","darkblue")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 1.5/1) +
  scale_x_discrete(labels=c('Rest', 'Exercise'))


#### Figure 3D - Clamp p-REPS1 S709 correlation with steady-state clamp ####
Clamp_REPS1_correlation <- read.delim("Clamp_REPS1_correlation.txt",dec=",")
shapiro.test(Clamp_REPS1_correlation$AUC) ## test for normality
shapiro.test(Clamp_REPS1_correlation$pReps1) ## test for normality

cor.test(Clamp_REPS1_correlation$pReps1, Clamp_REPS1_correlation$AUC, method="pearson")
(0.9289736)^2

ggplot(Clamp_REPS1_correlation, aes(x=pReps1, y=AUC)) + geom_point(aes(color = "darkred"),size=3) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.25, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1, color="black"),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16),legend.position = "none", aspect.ratio=1/1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) + 
  scale_shape_manual(values = c(16, 21)) +
  scale_color_manual(values=c("darkred")) +
  ggtitle("Glucose uptake") +
  xlab("Δ p-REPS1 S709 (Post-Pre clamp)") + ylab("leg glucose uptake \n (μmol/min/kg)") +
  annotate(geom="text",x=1.5,y=30, label="R2 = 0.86", color = "black",size=5) +
  annotate(geom = "text",x=1.5,y=20,label="p < 0.001",size=5)


#### Figure 3J - Reps1 KD GU data ####
KD_Reps1_data <- read.delim("GU_data_Reps1_KD.txt",dec=",")
KD_Reps1_data$KD <- as.factor(KD_Reps1_data$KD)
KD_Reps1_data$Insulin <- factor(KD_Reps1_data$Insulin, levels=c("Basal","Insulin"))

mean_df <- KD_Reps1_data %>%
  group_by(Insulin, KD) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = KD, y = MeanIntensity, fill = Insulin, color = Insulin),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = KD, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Insulin),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = KD_Reps1_data, aes(x = KD, y = Intensity, color = Insulin), 
             position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from siCTR Basal", title = "p-Reps1 S709") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("black", "darkred")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1) +
  scale_x_discrete(labels = c('siCTR', 'SiReps1'))

#### Figure 3I - Reps1 KD WB confirmation ####
Reps1_KD_wB <- read.delim("Reps1_KD_protein.txt",dec=",")
mean_df <- Reps1_KD_wB %>%
  group_by(KD) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = KD, y = MeanIntensity, fill=KD, color = KD),
           width=0.6, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = KD, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = KD),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = Reps1_KD_wB, aes(x = KD, y = Intensity, color = KD), 
             position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8), size = 5) +
  labs(x = "", y = "Fold change from siCTR Basal", title = "Reps1 protein") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("black", "darkred")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none"
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1.4/1) +
  scale_x_discrete(labels = c('siCTR', 'SiReps1')) + ylim(0,1.2)

#### Figure 5B left - LFD_HFD_data  BW ####
LFD_HFD_data <- read.delim("LFD_HFD_BW_data.txt",dec=",")
LFD_HFD_data$Diet <- factor(LFD_HFD_data$Diet, levels=c("LFD","HFD"))
mean_df <- LFD_HFD_data %>%
  group_by(Diet) %>%
  summarise(
    MeanIntensity = mean(BW, na.rm = TRUE),
    SEM = sd(BW, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Diet, y = MeanIntensity, fill=Diet, color = "black"),
           width=0.6, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Diet, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Diet),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = LFD_HFD_data, aes(x = Diet, y = BW, color = Diet), 
             position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8), size = 5) +
  labs(x = "", y = "Weigth (g)", title = "Body weight") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("black","grey", "#2ea2bf")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  theme(legend.title = element_blank(), aspect.ratio = 1.4/1)

#### Figure 5B right - Retro-Oribital (RO) saline/insulin injection LFD/HFD diet mice - delta BG values ####
RO_injection_data <- read.delim("LFD_HFD_Mouse_injection_BG.txt",dec=",")
RO_injection_data$Diet <- factor(RO_injection_data$Diet, levels=c("LFD","HFD"))
RO_injection_data$Injection <- factor(RO_injection_data$Injection, levels=c("Saline","Insulin"))

mean_df <- RO_injection_data %>%
  group_by(Diet,Injection) %>%
  summarise(
    MeanIntensity = mean(deltaBG, na.rm = TRUE),
    SEM = sd(deltaBG, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")


ggplot() + 
  #geom_bar(data = mean_df, aes(x = Diet, y = MeanIntensity, fill = Injection, color = Injection),
          # width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Diet, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Injection),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = RO_injection_data, aes(x = Diet, y = deltaBG, color = Injection), 
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "mmol/L \n (Delta 15-0 min)", title = "Blood glucose") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("grey", "#2ea2bf")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1) + ylim(-10,2)


#### Figure 5C - Tissue-specific response in ins-stim. Reps1 S709 phosphorylation ####
Tissue_Reps1_data <- read.delim("pREPS1_tissue_specific_response.txt",dec=",")
Tissue_Reps1_data$Diet <- factor(Tissue_Reps1_data$Diet, levels=c("LFD","HFD"))
Tissue_Reps1_data$Injection <- factor(Tissue_Reps1_data$Injection, levels=c("Saline","Insulin"))
Tissue_Reps1_data$Tissue <- factor(Tissue_Reps1_data$Tissue, levels=c("eWAT","Muscle","Liver"))
Tissue_Reps1_data$Interaction <- interaction(Tissue_Reps1_data$Diet, Tissue_Reps1_data$Injection, sep = "-")
Tissue_Reps1_data$Interaction <- factor(Tissue_Reps1_data$Interaction, levels = c("LFD-Saline", "LFD-Insulin", "HFD-Saline", "HFD-Insulin"))

mean_df <- Tissue_Reps1_data %>%
  group_by(Interaction, Diet, Tissue) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n())
  )

ggplot(data = mean_df, aes(x = Interaction, y = MeanIntensity, fill = Diet, color = Diet)) + 
  geom_bar(stat = "identity", position = "dodge", width=0.60) +
  geom_errorbar(aes(ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = Tissue_Reps1_data, 
             aes(x = Interaction, y = Intensity), 
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2) +
  labs(x = "", y = "mmol/L \n (Delta 15-0 min)", title = "p-Reps1 S709") +
  scale_fill_manual(values = c("LFD" = "white", "HFD" = "white")) +
  scale_color_manual(values = c("LFD" = "grey", "HFD" = "#2ea2bf")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"), aspect.ratio=2
  ) +
  theme(legend.title = element_blank()) +
  facet_wrap(~ Tissue, switch = "x", nrow = 1)


#### Figure 5D - Human T2D BMI and GIR ####
Human_T2D_data <- read.delim("T2D_BW_GIR_data.txt",dec=",")
Human_T2D_data$Group <- factor(Human_T2D_data$Group, levels=c("NGT","T2D"))

mean_df <- Human_T2D_data %>%
  group_by(Group) %>%
  summarise(
    MeanIntensity = mean(BW, na.rm = TRUE),
    SEM = sd(BW, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Group, y = MeanIntensity, fill=Group, color = "black"),
           width=0.6, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Group, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Group),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = Human_T2D_data, aes(x = Group, y = BW, color = Group), 
             position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8), size = 5) +
  labs(x = "", y = "kg/m2", title = "BMI") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("black","grey", "purple")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  theme(legend.title = element_blank(), aspect.ratio = 1.4/1) + ylim(0,40)                              


mean_df <- Human_T2D_data %>%
  group_by(Group) %>%
  summarise(
    MeanIntensity = mean(GIR, na.rm = TRUE),
    SEM = sd(GIR, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Group, y = MeanIntensity, fill=Group, color = "black"),
           width=0.6, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Group, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Group),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = Human_T2D_data, aes(x = Group, y = GIR, color = Group), 
             position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.8), size = 5) +
  labs(x = "", y = "μmol/kg/min", title = "Glucose infusion rate") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("black","grey", "purple")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "none") +
  theme(legend.title = element_blank(), aspect.ratio = 1.4/1)                              



#### Figure 5E - Human T2D p-REPS1 S709 WB ####
pREPS1_T2D_data <- read.delim("pREPS1_T2D_WB_Data.txt",dec=",")
pREPS1_T2D_data$Group <- factor(pREPS1_T2D_data$Group, levels=c("NGT","T2D"))
pREPS1_T2D_data$Condition <- factor(pREPS1_T2D_data$Condition, levels=c("Basal","Insulin"))
pREPS1_T2D_data$ID <- as.factor(pREPS1_T2D_data$ID)

mean_df <- pREPS1_T2D_data %>%
  group_by(Group,Condition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")


Basal_WB_data <- pREPS1_T2D_data %>% filter(Condition == "Basal")
Insulin_WB_data <- pREPS1_T2D_data %>% filter(Condition == "Insulin")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(Basal_WB_data, Insulin_WB_data, by = c("ID", "Group"), suffix = c("_Basal", "_Insulin"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("Basal" = "darkgrey", "Insulin" = "purple")  # Colors for points


# Create the plot
p <- ggplot() +
  geom_bar(data = mean_df, aes(x = Condition, y = MeanIntensity, fill = Condition, color=Condition), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  geom_point(data = pREPS1_T2D_data, aes(x = Condition, y = Intensity, color = Condition), 
             position = position_dodge(width = 0.8), size = 4) +
  geom_segment(data = line_data, aes(x = factor("Basal"), xend = factor("Insulin"), y = Intensity_Basal, yend = Intensity_Insulin, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  facet_grid(~ Group) +
  labs(x = "", y = "Fold change from NGT basal", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 3/1)

print(p)

#### Figure 5F - Human NGT/T2D GIR vs p-REPS1 S709 WB correlation ####
Corr_pREPS1_T2D_data <- read.delim("GIR_human_study_Reps1_correlation.txt",dec=",")
Corr_pREPS1_T2D_data$Group <- factor(Corr_pREPS1_T2D_data$Group, levels=c("NGT","T2D"))

NGT_data <- Corr_pREPS1_T2D_data[Corr_pREPS1_T2D_data$Group == "NGT",]
T2D_data <- Corr_pREPS1_T2D_data[Corr_pREPS1_T2D_data$Group == "T2D",]

shapiro.test(log2(T2D_data$GIR)) #chekc for normality
shapiro.test(NGT_data$GIR) #check for normality

cor.test(NGT_data$GIR,NGT_data$pReps1,method="pearson")
(0.8043437)^2
cor.test(T2D_data$GIR,T2D_data$pReps1,method="pearson")
(0.86395)^2

# ANCOVA model
model_interaction <- lm(GIR ~ Group * pReps1, data = Corr_pREPS1_T2D_data)
summary(model_interaction)

ggplot() +
  geom_point(data = NGT_data, aes(x = pReps1, y = GIR), color = "darkgrey",size=4) +
  geom_smooth(data = NGT_data, aes(x = pReps1, y = GIR), method = "lm", se = TRUE, color = "darkgrey",size=0.25,alpha=0.15) +
  geom_point(data = T2D_data, aes(x = pReps1, y = GIR), color = "purple",size=4) +
  geom_smooth(data = T2D_data, aes(x = pReps1, y = GIR), method = "lm", se = TRUE, color = "purple",size=0.25,alpha=0.15) +
  labs(title = "Insulin sensitivity association", x = "delta p-REPS1 S709 (Post-Pre)", y = "Glucose infusion rate (mg/min/m2)") +
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
        plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1, color="black"),
        axis.title = element_text(size = 18), legend.text = element_text(size=14),
        legend.title = element_text(size=16), aspect.ratio=1/1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate(geom="text",x=-1.5,y=700, label="R2 = 0.65", color = "black",size=3) +
  annotate(geom = "text",x=-1.5,y=675,label="p = 0.005",size=3) +
  annotate(geom="text",x=-1.2,y=450, label="R2 = 0.75", color = "black",size=3) +
  annotate(geom = "text",x=-1.2,y=400,label="p = 0.003",size=3) +
  annotate(geom="text",x=3,y=700, label="ANCOVA:\n Main-effect p-REPS1 S709 p < 0.001 \n Interaction p-REPS1 S709 x Group p = 0.024", color = "black",size=3)
  
  

#### Figure 4B - iRSK signaling data ####
iRSK_signaling_data <- read.delim("iRSK_signaling_data.txt",dec=",")

iRSK_signaling_data$Stimulation <- factor(iRSK_signaling_data$Stimulation, levels = c("Basal","Insulin"))
iRSK_signaling_data$Inhibition <- factor(iRSK_signaling_data$Inhibition, levels=c("DMSO","BI-D1870"))

mean_df <- iRSK_signaling_data %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
          width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iRSK_signaling_data, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-Reps1 S709") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)

#### Figure 4A - iRSK inhibition and glucose uptake ####
iRSK_GU_data <- read.delim("iRSK_GU_data.txt",dec=",")
iRSK_GU_data$Stimulation <- factor(iRSK_GU_data$Stimulation, levels = c("Basal","Insulin"))
iRSK_GU_data$Treatment <- factor(iRSK_GU_data$Treatment, levels=c("DMSO","BI-D1870"))

mean_df <- iRSK_GU_data %>%
  group_by(Stimulation, Treatment) %>%
  summarise(
    MeanIntensity = mean(GU, na.rm = TRUE),
    SEM = sd(GU, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Treatment, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Treatment, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iRSK_GU_data, aes(x = Treatment, y = GU, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "Glucose uptake") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)

#### Figure 4S1A iRSK inhibition pAkt S473 signaling ####
iRSK_signaling_data_AKT <- read.delim("pAKT_WB_C2C12_BID1870.txt",dec=",")

iRSK_signaling_data_AKT$Stimulation <- factor(iRSK_signaling_data_AKT$Stimulation, levels = c("Basal","Insulin"))
iRSK_signaling_data_AKT$Inhibition <- factor(iRSK_signaling_data_AKT$Treatment, levels=c("DMSO","BI-D1870"))

mean_df <- iRSK_signaling_data_AKT %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iRSK_signaling_data_AKT, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-Akt S473") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)

#### Figure 4S1B - pP70S6K ####
mTOR_path_signaling_data_P70S6K <- read.delim("pP70S6K_T389_WB_C2C12_rapa.txt",dec=",")

mTOR_path_signaling_data_P70S6K$Stimulation <- factor(mTOR_path_signaling_data_P70S6K$Stimulation, levels = c("Basal","Insulin"))
mTOR_path_signaling_data_P70S6K$Inhibition <- factor(mTOR_path_signaling_data_P70S6K$Treatment, levels=c("DMSO","Rapamycin","PF-4708671"))

mean_df <- mTOR_path_signaling_data_P70S6K %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = mTOR_path_signaling_data_P70S6K, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-P70S6K T389") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)


#### Figure 4S1C - pREPS1 S709 mTOR sign pathway inhibition ####
mTOR_path_signaling_data_Reps1 <- read.delim("Mtor_path_inhb_pReps1.txt",dec=",")

mTOR_path_signaling_data_Reps1$Stimulation <- factor(mTOR_path_signaling_data_Reps1$Stimulation, levels = c("Basal","Insulin"))
mTOR_path_signaling_data_Reps1$Inhibition <- factor(mTOR_path_signaling_data_Reps1$Treatment, levels=c("DMSO","Rapamycin","PF-4708671"))

mean_df <- mTOR_path_signaling_data_Reps1 %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = mTOR_path_signaling_data_Reps1, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-REPS1 S709") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)



#### Figure 4S1I - iMAPK inhibition WB Reps1 S709 ####
iErk_signaling_data_Reps1 <- read.delim("iMAPK_WB_Reps1.txt",dec=",")
iErk_signaling_data_Reps1$Stimulation <- factor(iErk_signaling_data_Reps1$Stimulation, levels = c("Basal","Insulin"))
iErk_signaling_data_Reps1$Inhibition <- factor(iErk_signaling_data_Reps1$Treatment, levels=c("DMSO","U0126"))

mean_df <- iErk_signaling_data_Reps1 %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iErk_signaling_data_Reps1, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-REPS1 S709") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)

#### Figure 4S1H - iMAPK inhibition WB ERK1/2 T202/Y204 ####
iErk_signaling_data_ERK <- read.delim("iMAPK_WB_Erk.txt",dec=",")

iErk_signaling_data_ERK$Stimulation <- factor(iErk_signaling_data_ERK$Stimulation, levels = c("Basal","Insulin"))
iErk_signaling_data_ERK$Inhibition <- factor(iErk_signaling_data_ERK$Treatment, levels=c("DMSO","U0126"))

mean_df <- iErk_signaling_data_ERK %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iErk_signaling_data_ERK, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-ERK1/2 T202/Y204") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)

#### Figure 4S1F - iPDK1 WB pREPS1 S709 ####
iPDK1_signaling_data_Reps1 <- read.delim("iPDK1_pReps1.txt",dec=",")

iPDK1_signaling_data_Reps1$Stimulation <- factor(iPDK1_signaling_data_Reps1$Stimulation, levels = c("Basal","Insulin"))
iPDK1_signaling_data_Reps1$Inhibition <- factor(iPDK1_signaling_data_Reps1$Treatment, levels=c("DMSO","MK-7"))

mean_df <- iPDK1_signaling_data_Reps1 %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iPDK1_signaling_data_Reps1, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-REPS1 S709") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)

#### Figure 4S1E - iPDK1 WB AKT T308 ####
iPDK1_signaling_data_AKT <- read.delim("New_R_TINX_data_inputs/iPDK1_pAKT_T308.txt",dec=",")

iPDK1_signaling_data_AKT$Stimulation <- factor(iPDK1_signaling_data_AKT$Stimulation, levels = c("Basal","Insulin"))
iPDK1_signaling_data_AKT$Inhibition <- factor(iPDK1_signaling_data_AKT$Treatment, levels=c("DMSO","MK-7"))

mean_df <- iPDK1_signaling_data_AKT %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iPDK1_signaling_data_AKT, aes(x = Inhibition, y = Intensity, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3) +
  labs(x = "", y = "Fold change from Basal DMSO", title = "p-AKT T308") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)


#### Figure 4G - ExVivo iRSK inhibition pAkt S473 signaling ####
pAKT_exvivo_data <- read.delim("Exvivo_pAKT_S473.txt", dec=",")
pAKT_exvivo_data$ID <- factor(pAKT_exvivo_data$ID)

# Adjusting factor levels
pAKT_exvivo_data$Stimulation <- factor(pAKT_exvivo_data$Stimulation, levels = c("Basal", "Insulin"))
pAKT_exvivo_data$Treatment <- factor(pAKT_exvivo_data$Treatment, levels = c("DMSO", "BI-D1870"))

mean_df <- pAKT_exvivo_data %>%
  group_by(Stimulation, Treatment) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")


# Split data into 'Rest' and 'Contraction'
basal_data <- pAKT_exvivo_data %>% filter(Treatment == "DMSO")
treatment_data <- pAKT_exvivo_data %>% filter(Treatment == "BI-D1870")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(basal_data, treatment_data, by = c("ID", "Stimulation"), suffix = c("_DMSO", "_BI_D1870"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("DMSO" = "darkgrey", "BI-D1870" = "darkgreen")  # Colors for points


# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Treatment, y = MeanIntensity, fill = Treatment, color=Treatment), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = pAKT_exvivo_data, aes(x = Treatment, y = Intensity, color = Treatment), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("DMSO"), xend = factor("BI-D1870"), y = Intensity_DMSO, yend = Intensity_BI_D1870, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Stimulation) +
  
  # Customizations
  labs(x = "", y = "Fold change from basal DMSO", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 2/1)

print(p)


#### Figure 4F - ExVivo iRSK inhibition pReps1 S709 / Reps1 total signaling ####
pReps1_exvivo_data <- read.delim("Exvivo_pReps1_S709.txt", dec=",")
pReps1_exvivo_data$ID <- factor(pAKT_exvivo_data$ID)

# Adjusting factor levels
pReps1_exvivo_data$Stimulation <- factor(pReps1_exvivo_data$Stimulation, levels = c("Basal", "Insulin"))
pReps1_exvivo_data$Treatment <- factor(pReps1_exvivo_data$Treatment, levels = c("DMSO", "BI-D1870"))

mean_df <- pReps1_exvivo_data %>%
  group_by(Stimulation, Treatment) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")


# Split data into 'Rest' and 'Contraction'
basal_data <- pReps1_exvivo_data %>% filter(Treatment == "DMSO")
treatment_data <- pReps1_exvivo_data %>% filter(Treatment == "BI-D1870")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(basal_data, treatment_data, by = c("ID", "Stimulation"), suffix = c("_DMSO", "_BI_D1870"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("DMSO" = "darkgrey", "BI-D1870" = "darkgreen")  # Colors for points


# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Treatment, y = MeanIntensity, fill = Treatment, color=Treatment), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = pReps1_exvivo_data, aes(x = Treatment, y = Intensity, color = Treatment), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("DMSO"), xend = factor("BI-D1870"), y = Intensity_DMSO, yend = Intensity_BI_D1870, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Stimulation) +
  
  # Customizations
  labs(x = "", y = "Fold change from basal DMSO", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 2/1)

print(p)

#### Figure 4S2A - Pilot experiment ExVivo EDL contraction delta pREPS1 ####
# Load data
delta_contraction_exvivo_data <- read.delim("EDL_Reps1_contraction.txt", dec=",")
delta_contraction_exvivo_data$ID <- factor(delta_contraction_exvivo_data$ID)

# Ensure "Rest" comes first in factor levels
delta_contraction_exvivo_data$Contraction <- factor(delta_contraction_exvivo_data$Contraction, levels = c("Rest", "Contraction"))

# Keep only Rest and Contraction conditions
filtered_data <- delta_contraction_exvivo_data %>%
  filter(Contraction %in% c("Rest", "Contraction"))

# Compute mean intensities and SEM
mean_df <- filtered_data %>%
  group_by(Contraction) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),  # Calculate SEM
    .groups = "drop"
  )

# Split for paired lines
rest_data <- filtered_data %>% filter(Contraction == "Rest")
contraction_data <- filtered_data %>% filter(Contraction == "Contraction")

# Merge based on ID to create lines
line_data <- left_join(rest_data, contraction_data, by = "ID", suffix = c("_Rest", "_Contraction"))

# Define colors
point_colors <- c("Rest" = "darkgrey", "Contraction" = "#d99f00")

# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Contraction, y = MeanIntensity, fill = Contraction, color = Contraction), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = filtered_data, aes(x = Contraction, y = Intensity, color = Contraction), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points (Rest → Contraction)
  geom_segment(data = line_data, aes(x = "Rest", xend = "Contraction", 
                                     y = Intensity_Rest, yend = Intensity_Contraction, group = ID), 
               color = "black", alpha = 0.6, size = 0.3) +
  
  # Customizations
  labs(x = "", y = "Intensity", title = "") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = point_colors) +
  
  # Adjust theme to remove legend and add custom x-axis labels
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12),  # Show x-axis labels
        axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        aspect.ratio = 2/1) +
  
  # Manually set x-axis labels
  scale_x_discrete(labels = c("Rest" = "Rest", "Contraction" = "Contraction"))

print(p)


#### Figure 4S2B - Pilot experiment ExVivo Soleus contraction delta pREPS1 ####
# Load data
delta_contraction_exvivo_data <- read.delim("Soleus_Reps1_contraction.txt", dec=",")
delta_contraction_exvivo_data$ID <- factor(delta_contraction_exvivo_data$ID)

# Ensure "Rest" comes first in factor levels
delta_contraction_exvivo_data$Contraction <- factor(delta_contraction_exvivo_data$Contraction, levels = c("Rest", "Contraction"))

# Keep only Rest and Contraction conditions
filtered_data <- delta_contraction_exvivo_data %>%
  filter(Contraction %in% c("Rest", "Contraction"))

# Compute mean intensities and SEM
mean_df <- filtered_data %>%
  group_by(Contraction) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    SEM = sd(Intensity, na.rm = TRUE) / sqrt(n()),  # Calculate SEM
    .groups = "drop"
  )

# Split for paired lines
rest_data <- filtered_data %>% filter(Contraction == "Rest")
contraction_data <- filtered_data %>% filter(Contraction == "Contraction")

# Merge based on ID to create lines
line_data <- left_join(rest_data, contraction_data, by = "ID", suffix = c("_Rest", "_Contraction"))

# Define colors
point_colors <- c("Rest" = "darkgrey", "Contraction" = "#d99f00")

# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Contraction, y = MeanIntensity, fill = Contraction, color = Contraction), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = filtered_data, aes(x = Contraction, y = Intensity, color = Contraction), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points (Rest → Contraction)
  geom_segment(data = line_data, aes(x = "Rest", xend = "Contraction", 
                                     y = Intensity_Rest, yend = Intensity_Contraction, group = ID), 
               color = "black", alpha = 0.6, size = 0.3) +
  
  # Customizations
  labs(x = "", y = "Intensity", title = "") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = point_colors) +
  
  # Adjust theme to remove legend and add custom x-axis labels
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12),  # Show x-axis labels
        axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        aspect.ratio = 2/1) +
  
  # Manually set x-axis labels
  scale_x_discrete(labels = c("Rest" = "Rest", "Contraction" = "Contraction"))

print(p)


#### Figure 4J - ExVivo EDL contraction pReps1 S709/Reps1 ####
pReps1_exvivo_data <- read.delim("ExVivo_contraction_pREPS1_REPS1.txt", dec=",")
pReps1_exvivo_data$ID <- factor(pReps1_exvivo_data$ID)

# Adjusting factor levels
pReps1_exvivo_data$Stimulation <- factor(pReps1_exvivo_data$Stimulation, levels = c("Rest", "Contraction"))
pReps1_exvivo_data$Treatment <- factor(pReps1_exvivo_data$Treatment, levels = c("DMSO", "BI-D1870"))

mean_df <- pReps1_exvivo_data %>%
  group_by(Stimulation, Treatment) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")


# Split data into 'Rest' and 'Contraction'
basal_data <- pReps1_exvivo_data %>% filter(Treatment == "DMSO")
treatment_data <- pReps1_exvivo_data %>% filter(Treatment == "BI-D1870")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(basal_data, treatment_data, by = c("ID", "Stimulation"), suffix = c("_DMSO", "_BI_D1870"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("DMSO" = "darkgrey", "BI-D1870" = "#d99f00")  # Colors for points


# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Treatment, y = MeanIntensity, fill = Treatment, color=Treatment), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = pReps1_exvivo_data, aes(x = Treatment, y = Intensity, color = Treatment), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("DMSO"), xend = factor("BI-D1870"), y = Intensity_DMSO, yend = Intensity_BI_D1870, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Stimulation) +
  
  # Customizations
  labs(x = "", y = "Fold change from rest DMSO", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 2/1)

print(p)



#### Figure 4S2D - ExVivo EDL contraction delta glucose uptake (contraction - rest) ####
delta_contraction_exvivo_data <- read.delim("Delta_contraction_rest_GU.txt", dec=",")
delta_contraction_exvivo_data$ID <- factor(delta_contraction_exvivo_data$ID)

# Ensure DMSO comes first in factor levels
delta_contraction_exvivo_data$Inhibitor <- factor(delta_contraction_exvivo_data$Inhibitor, levels = c("DMSO", "BI-D1870"))

# Keep only DMSO and BI-D1870
filtered_data <- delta_contraction_exvivo_data %>%
  filter(Inhibitor %in% c("DMSO", "BI-D1870"))

# Compute mean intensities
mean_df <- filtered_data %>%
  group_by(Inhibitor) %>%
  summarise(MeanIntensity = mean(Delta.GU, na.rm = TRUE), .groups = "drop")

# Split for paired lines
basal_data <- filtered_data %>% filter(Inhibitor == "DMSO")
treatment_data <- filtered_data %>% filter(Inhibitor == "BI-D1870")

# Merge based on ID to create lines
line_data <- left_join(basal_data, treatment_data, by = "ID", suffix = c("_DMSO", "_BI_D1870"))

# Define colors
point_colors <- c("DMSO" = "darkgrey", "BI-D1870" = "#d99f00")

# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Inhibitor, y = MeanIntensity, fill = Inhibitor, color = Inhibitor), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = filtered_data, aes(x = Inhibitor, y = Delta.GU, color = Inhibitor), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points
  geom_segment(data = line_data, aes(x = "DMSO", xend = "BI-D1870", 
                                     y = Delta.GU_DMSO, yend = Delta.GU_BI_D1870, group = ID), 
               color = "black", alpha = 0.6, size = 0.3) +
  
  # Customizations
  labs(x = "", y = "Δ GU (Contraction - Rest)", title = "") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = point_colors) +
  
  # Adjust theme to remove legend and add custom x-axis labels
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12),  # Show x-axis labels
        axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        aspect.ratio = 2/1) +
  
  # Manually set x-axis labels
  scale_x_discrete(labels = c("DMSO" = "DMSO", "BI-D1870" = "BI-D1870"))

print(p)




#### Figure 4S2C - ExVivo EDL contraction delta glucose uptake (DMSO - BI-D8170) ####
delta_contraction_exvivo_data <- read.delim("Delta_inhibitor_GU.txt", dec=",")
delta_contraction_exvivo_data$ID <- factor(delta_contraction_exvivo_data$ID)

# Ensure "Rest" comes first in factor levels
delta_contraction_exvivo_data$Contraction <- factor(delta_contraction_exvivo_data$Contraction, levels = c("Rest", "Contraction"))

# Filter data to keep only Rest and Contraction conditions
filtered_data <- delta_contraction_exvivo_data %>%
  filter(Contraction %in% c("Rest", "Contraction"))

# Compute mean and SEM
mean_df <- filtered_data %>%
  group_by(Contraction) %>%
  summarise(
    MeanIntensity = mean(Delta.GU, na.rm = TRUE),
    SEM = sd(Delta.GU, na.rm = TRUE) / sqrt(n()),  # Calculate SEM
    .groups = "drop"
  )

# Define colors
point_colors <- c("Rest" = "darkgrey", "Contraction" = "#d99f00")

# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Contraction, y = MeanIntensity, fill = Contraction, color = Contraction), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Error bars for SEM
  geom_errorbar(data = mean_df, aes(x = Contraction, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color=Contraction), 
                width = 0.2, position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = filtered_data, aes(x = Contraction, y = Delta.GU, color = Contraction), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Customizations
  labs(x = "", y = "Δ GU (Contraction - Rest)", title = "") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = point_colors) +
  
  # Adjust theme to remove legend and add custom x-axis labels
  theme_minimal() +
  theme(legend.position = "none",  # Remove legend
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12),  # Show x-axis labels
        axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        aspect.ratio = 2/1) +
  
  # Manually set x-axis labels
  scale_x_discrete(labels = c("Rest" = "Rest", "Contraction" = "Contraction"))

print(p)





#### Figure 4K - ExVivo EDL contraction pAMPK T172/AMPKa ####
pReps1_exvivo_data <- read.delim("ExVivo_contraction_pAMPK_AMPK.txt", dec=",")
pReps1_exvivo_data$ID <- factor(pReps1_exvivo_data$ID)

# Adjusting factor levels
pReps1_exvivo_data$Stimulation <- factor(pReps1_exvivo_data$Stimulation, levels = c("Rest", "Contraction"))
pReps1_exvivo_data$Treatment <- factor(pReps1_exvivo_data$Treatment, levels = c("DMSO", "BI-D1870"))

mean_df <- pReps1_exvivo_data %>%
  group_by(Stimulation, Treatment) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")


# Split data into 'Rest' and 'Contraction'
basal_data <- pReps1_exvivo_data %>% filter(Treatment == "DMSO")
treatment_data <- pReps1_exvivo_data %>% filter(Treatment == "BI-D1870")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(basal_data, treatment_data, by = c("ID", "Stimulation"), suffix = c("_DMSO", "_BI_D1870"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("DMSO" = "darkgrey", "BI-D1870" = "#d99f00")  # Colors for points


# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Treatment, y = MeanIntensity, fill = Treatment, color=Treatment), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = pReps1_exvivo_data, aes(x = Treatment, y = Intensity, color = Treatment), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("DMSO"), xend = factor("BI-D1870"), y = Intensity_DMSO, yend = Intensity_BI_D1870, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Stimulation) +
  
  # Customizations
  labs(x = "", y = "Fold change from rest DMSO", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 2/1)

print(p)


#### Figure 4E - ExVivo GU soleus muscle ####
GU_data_in_situ <- read.delim("Exvivo_ISGU_soleus.txt", dec=",")
GU_data_in_situ$ID <- factor(GU_data_in_situ$ID)

# Adjusting factor levels
GU_data_in_situ$Stimulation <- factor(GU_data_in_situ$Stimulation, levels = c("Basal", "Insulin"))
GU_data_in_situ$Treatment <- factor(GU_data_in_situ$Treatment, levels = c("DMSO", "BI-D1870"))

mean_df <- GU_data_in_situ %>%
  group_by(Stimulation, Treatment) %>%
  summarise(MeanIntensity = mean(GU, na.rm = TRUE), .groups = "drop")


# Split data into 'Rest' and 'Contraction'
basal_data <- GU_data_in_situ %>% filter(Treatment == "DMSO")
treatment_data <- GU_data_in_situ %>% filter(Treatment == "BI-D1870")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(basal_data, treatment_data, by = c("ID", "Stimulation"), suffix = c("_DMSO", "_BI_D1870"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("DMSO" = "darkgrey", "BI-D1870" = "darkgreen")  # Colors for points


# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Treatment, y = MeanIntensity, fill = Treatment, color=Treatment), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = GU_data_in_situ, aes(x = Treatment, y = GU, color = Treatment), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("DMSO"), xend = factor("BI-D1870"), y = GU_DMSO, yend = GU_BI_D1870, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Stimulation) +
  
  # Customizations
  labs(x = "", y = "umol/g/hr", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 2/1)

print(p)



#### Figure 4I - ExVivo GU EDL muscle contraction ####
GU_data_in_situ <- read.delim("ExVivo_contraction_GU.txt", dec=",")
GU_data_in_situ$ID <- factor(GU_data_in_situ$ID)

# Adjusting factor levels
GU_data_in_situ$Stimulation <- factor(GU_data_in_situ$Stimulation, levels = c("Rest", "Contraction"))
GU_data_in_situ$Treatment <- factor(GU_data_in_situ$Treatment, levels = c("DMSO", "BI-D1870"))

mean_df <- GU_data_in_situ %>%
  group_by(Stimulation, Treatment) %>%
  summarise(MeanIntensity = mean(GU, na.rm = TRUE), .groups = "drop")


# Split data into 'Rest' and 'Contraction'
basal_data <- GU_data_in_situ %>% filter(Treatment == "DMSO")
treatment_data <- GU_data_in_situ %>% filter(Treatment == "BI-D1870")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(basal_data, treatment_data, by = c("ID", "Stimulation"), suffix = c("_DMSO", "_BI_D1870"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("DMSO" = "darkgrey", "BI-D1870" = "#d99f00")  # Colors for points


# Create the plot
p <- ggplot() +
  
  # Bars for mean intensities
  geom_bar(data = mean_df, aes(x = Treatment, y = MeanIntensity, fill = Treatment, color=Treatment), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  
  # Points for individual values
  geom_point(data = GU_data_in_situ, aes(x = Treatment, y = GU, color = Treatment), 
             position = position_dodge(width = 0.8), size = 4) +
  
  # Lines connecting paired data points for each mouse
  geom_segment(data = line_data, aes(x = factor("DMSO"), xend = factor("BI-D1870"), y = GU_DMSO, yend = GU_BI_D1870, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  
  # Separate the data based on Saline and Insulin and switch the facet labels to the bottom
  facet_grid(~ Stimulation) +
  
  # Customizations
  labs(x = "", y = "umol/g/hr", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 2/1)

print(p)



#### iRSK GLUT4-translocation
iRSK_GLUT4_data_AKT <- read.delim("GLUT4_translocation_iRSK.txt",dec=",")

iRSK_GLUT4_data_AKT$Stimulation <- factor(iRSK_GLUT4_data_AKT$Stimulation, levels = c("Basal","Insulin"))
iRSK_GLUT4_data_AKT$Inhibition <- factor(iRSK_GLUT4_data_AKT$Treatment, levels=c("DMSO","BI-D1870"))

mean_df <- iRSK_GLUT4_data_AKT %>%
  group_by(Stimulation, Inhibition) %>%
  summarise(
    MeanIntensity = mean(Absorbance, na.rm = TRUE),
    SEM = sd(Absorbance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")

ggplot() + 
  geom_bar(data = mean_df, aes(x = Inhibition, y = MeanIntensity, fill = Stimulation, color = Stimulation),
           width=0.60, stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(data = mean_df, aes(x = Inhibition, ymin = MeanIntensity - SEM, ymax = MeanIntensity + SEM, color = Stimulation),
                position = position_dodge(0.8), width = 0.25) +
  geom_point(data = iRSK_GLUT4_data_AKT, aes(x = Inhibition, y = Absorbance, color = Stimulation), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 3, alpha=0.3) +
  labs(x = "", y = "Absorbance (AU)", title = "GLUT4 translocation") +
  scale_fill_manual(values = c("white", "white")) +
  scale_color_manual(values = c("darkgrey", "darkblue")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")
  ) +
  theme(legend.title = element_blank(), aspect.ratio = 1/1)



#### Figure 1S1F - Human insulin levels during exercise ####
Human_plasma_insulin_Exercise <- read.delim("Plasma_insulin_Exercise.txt",dec = ",")
Human_plasma_insulin_Exercise$ID <- as.factor(Human_plasma_insulin_Exercise$id)
Human_plasma_insulin_Exercise$Exercise <- factor(Human_plasma_insulin_Exercise$time, levels=c("pre","ex"))

mean_df <- Human_plasma_insulin_Exercise %>%
  group_by(Exercise) %>%
  summarise(MeanIntensity = mean(ins, na.rm = TRUE), .groups = "drop")

ggplot(data=Human_plasma_insulin_Exercise, aes(x=Exercise,y=ins)) +
  geom_bar(data = mean_df, aes(y = MeanIntensity, fill = Exercise, color=Exercise),width=0.60, 
           stat = "identity", position = position_dodge(width = 0.8)) + 
  geom_point(aes(color=Exercise), size=5) +
  geom_line(aes(group=ID),alpha = 0.6, size=0.3) +
  labs(x = "", y = "pmol/L", title = "Plasma insulin") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = c("black","darkblue")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 1.5/1) +
  scale_x_discrete(labels=c('Rest', 'Exercise'))

#### Figure 1S1A & B - Human cathecolamine levels during exercise ####
Human_plasma_cath_Exercise <- read.delim("Plasma_Cathe_Exercise.txt",dec = ",")

Human_plasma_cath_Exercise$ID <- as.factor(Human_plasma_cath_Exercise$id)
Human_plasma_cath_Exercise$Exercise <- factor(Human_plasma_cath_Exercise$time, levels=c("pre","post"))

mean_df <- Human_plasma_cath_Exercise %>%
  group_by(Exercise) %>%
  summarise(MeanIntensity = mean(epi, na.rm = TRUE), .groups = "drop")

ggplot(data=Human_plasma_cath_Exercise, aes(x=Exercise,y=epi)) +
  geom_bar(data = mean_df, aes(y = MeanIntensity, fill = Exercise, color=Exercise),width=0.60, 
           stat = "identity", position = position_dodge(width = 0.8)) + 
  geom_point(aes(color=Exercise), size=5) +
  geom_line(aes(group=ID),alpha = 0.6, size=0.3) +
  labs(x = "", y = "pmol/L", title = "Plasma epinephrine") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = c("black","darkblue")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 1.5/1) +
  scale_x_discrete(labels=c('Rest', 'Exercise'))

mean_df <- Human_plasma_cath_Exercise %>%
  group_by(Exercise) %>%
  summarise(MeanIntensity = mean(nor, na.rm = TRUE), .groups = "drop")

ggplot(data=Human_plasma_cath_Exercise, aes(x=Exercise,y=nor)) +
  geom_bar(data = mean_df, aes(y = MeanIntensity, fill = Exercise, color=Exercise),width=0.60, 
           stat = "identity", position = position_dodge(width = 0.8)) + 
  geom_point(aes(color=Exercise), size=5) +
  geom_line(aes(group=ID),alpha = 0.6, size=0.3) +
  labs(x = "", y = "pmol/L", title = "Plasma norepinephrine") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = c("black","darkblue")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 1.5/1) +
  scale_x_discrete(labels=c('Rest', 'Exercise'))


#### Figure 1S1E - Human plasma insulin during clamp ####
Human_plasma_Insulin_clamp <- read.delim("Plasma_insulin_Clamp.txt",dec = ",")
Human_plasma_Insulin_clamp <- Human_plasma_Insulin_clamp[!(Human_plasma_Insulin_clamp$id == "tinx03"),] #no baseline insulin measurement, therefore excluded
Human_plasma_Insulin_clamp$ID <- as.factor(Human_plasma_Insulin_clamp$id)
Human_plasma_Insulin_clamp$time <- factor(Human_plasma_Insulin_clamp$time, levels=c("0","30","60","90","120"))


ggplot(Human_plasma_Insulin_clamp, aes(x=time,y=ins)) + geom_point(color="darkblue",size=5) + geom_line(aes(group=id),color="darkblue",alpha = 0.6, size=0.3) +
  labs(x = "Time", y = "pmol/L", title = "Plasma insulin during clamp") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank())


#### Figure 1S1C & D - Human plasma cathecolamines during clamp ####
Human_plasma_cath_clamp <- read.delim("Plasma_Cathe_Insulin.txt",dec = ",")
Human_plasma_cath_clamp <- Human_plasma_cath_clamp[!(Human_plasma_cath_clamp$id == "tinx03"),] #no baseline insulin measurement, therefore excluded
Human_plasma_cath_clamp$ID <- as.factor(Human_plasma_cath_clamp$id)
Human_plasma_cath_clamp$time <- factor(Human_plasma_cath_clamp$time, levels=c("0","30","60","90","120"))

#export in 7x5 landscape
ggplot(Human_plasma_cath_clamp, aes(x=time,y=epi)) + geom_point(color="darkblue",size=5) + geom_line(aes(group=id),color="darkblue",alpha = 0.6, size=0.3) +
  labs(x = "Time", y = "pmol/L", title = "Plasma epinephrine during clamp") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank())

ggplot(Human_plasma_cath_clamp, aes(x=time,y=nor)) + geom_point(color="darkblue",size=5) + geom_line(aes(group=id),color="darkblue",alpha = 0.6, size=0.3) +
  labs(x = "Time", y = "pmol/L", title = "Plasma norepinephrine during clamp") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank())



#### Figure 3S1C & D - Human exercise GU pReps1 levels ####
Human_Exercise_GU_reps1 <- read.delim("Reps1_phos_exercise_GU.txt",dec = ",")

Human_Exercise_GU_reps1$ID <- as.factor(Human_Exercise_GU_reps1$id)
Human_Exercise_GU_reps1$Exercise <- factor(Human_Exercise_GU_reps1$time, levels=c("rest","ex"))

mean_df <- Human_Exercise_GU_reps1 %>%
  group_by(Exercise) %>%
  summarise(MeanIntensity = mean(normalised, na.rm = TRUE), .groups = "drop")

ggplot(data=Human_Exercise_GU_reps1, aes(x=Exercise,y=normalised)) +
  geom_bar(data = mean_df, aes(y = MeanIntensity, fill = Exercise, color=Exercise),width=0.60, 
           stat = "identity", position = position_dodge(width = 0.8)) + 
  geom_point(aes(color=Exercise), size=5) +
  geom_line(aes(group=ID),alpha = 0.6, size=0.3) +
  labs(x = "", y = "p-REPS1 S709/REPS1", title = "") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = c("black","darkblue")) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 1.5/1) +
  scale_x_discrete(labels=c('Rest', 'Exercise'))


Human_Exercise_GU_recovery_reps1 <- read.delim("New_R_TINX_data_inputs/GU_recovery_ex_human_reps1.txt",dec = ",")
Human_Exercise_GU_recovery_reps1$ID <- as.factor(Human_Exercise_GU_recovery_reps1$id)

cor.test(In_situ_delta_data_SALINE$delta_GU,In_situ_delta_data_SALINE$scaled_phos,method="pearson")


-(0.2535771)^2

ggplot(Human_Exercise_GU_recovery_reps1, aes(x=pReps1, y=GU_rec)) + geom_point(size=5,color="darkblue") +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.25, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1,color="black"),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16), aspect.ratio=1/1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) +
  ggtitle("Leg glucose uptake") +
  xlab("p-REPS1 S709 / REPS1") + ylab("mmol/min/kg thigh lean mass") +
  annotate(geom="text",x=1,y=1.8, label="R2 = 0.01", color = "black",size=5) +
  annotate(geom = "text",x=1,y=1.5,label="p = 0.85",size=5)





#### Figure 5S1A - LFD/HFD Ins-stimulated BG drop and p-Reps1 S709 correlation ####
LFD_HFD_BG_Reps1 <- read.delim("LFD_HFD_BG_pReps1_cor.txt",dec=",")
LFD_HFD_BG_Reps1$Diet <- factor(LFD_HFD_BG_Reps1$Diet, levels=c("LFD","HFD"))

shapiro.test(LFD_HFD_BG_Reps1$BG) ## test for normality
shapiro.test(LFD_HFD_BG_Reps1$pReps1) ## test for normality

cor.test(LFD_HFD_BG_Reps1$pReps1, LFD_HFD_BG_Reps1$BG, method="pearson")
(0.8324341)^2

ggplot(LFD_HFD_BG_Reps1, aes(x=pReps1, y=BG)) + geom_point(aes(color = Diet),size=3) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.25, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1, color="black"),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16),legend.position = "none", aspect.ratio=1/1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) + 
  scale_shape_manual(values = c(16, 21)) +
  scale_color_manual(values=c("grey","#2ea2bf")) +
  xlab("p-Reps1 S709") + ylab("ΔBG 15-0 min \n (mmol/L)") +
  annotate(geom="text",x=5,y=4, label="R2 = 0.69", color = "black",size=5) +
  annotate(geom = "text",x=5,y=3.5,label="p = 0.01",size=5)


#### Figure 3S1A & B - WB pAKT and pMTOR correlation steady-state Leg glucose uptake ####
pATK_correlation <- read.delim("pAKT_GIR_corelation_TINX.txt",dec=",")
pMTOR_correlation <- read.delim("pMTOR_GIR_correlation_TINX.txt",dec=",")

shapiro.test(pATK_correlation$pAKT) ## test for normality
shapiro.test(pMTOR_correlation$pMTOR) ## test for normality

cor.test(pATK_correlation$pAKT, pATK_correlation$AUC, method="pearson")
(0.2643036)^2

cor.test(pMTOR_correlation$pMTOR, pMTOR_correlation$AUC, method="pearson")
(0.2172571)^2

ggplot(pATK_correlation, aes(x=pAKT, y=AUC)) + geom_point(aes(color = "darkred"),size=3) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.25, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1, color="black"),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16),legend.position = "none", aspect.ratio=1/1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) + 
  scale_shape_manual(values = c(16, 21)) +
  scale_color_manual(values=c("darkred")) +
  ggtitle("Glucose uptake") +
  xlab("Δ p-AKT S473 (Post-Pre clamp)") + ylab("Steady-state AUC \n (umol/min/kg)") +
  annotate(geom="text",x=3,y=40, label="R2 = 0.07", color = "black",size=5) +
  annotate(geom = "text",x=3,y=30,label="p = 0.53",size=5)

ggplot(pMTOR_correlation, aes(x=pMTOR, y=AUC)) + geom_point(aes(color = "darkred"),size=3) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", size = 0.25, alpha=0.15, fill="#808096") +
  theme_minimal() + theme(axis.line = element_line(linewidth = 0.1, colour = "black"),
                          plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size=16, hjust = 1, color="black"),
                          axis.title = element_text(size = 18), legend.text = element_text(size=14),
                          legend.title = element_text(size=16),legend.position = "none", aspect.ratio=1/1,
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank()) + 
  scale_shape_manual(values = c(16, 21)) +
  scale_color_manual(values=c("darkred")) +
  ggtitle("Glucose uptake") +
  xlab("Δ p-MTOR S2448 (Post-Pre clamp)") + ylab("Steady-state AUC \n (umol/min/kg)") +
  annotate(geom="text",x=1,y=75, label="R2 = 0.047", color = "black",size=5) +
  annotate(geom = "text",x=1,y=65,label="p = 0.61",size=5)

#### Figure 1SG-J - WB quantification of insulin- and exercise-induced signaling in human skeletal muscle ####
Quant_human_WB_data <- read.delim("WB_quantified_Figure_S1.txt", dec=",")
Quant_human_WB_data$ID <- as.factor(Quant_human_WB_data$ID)
Quant_human_WB_data$Intervention <- factor(Quant_human_WB_data$Intervention,levels=c("Insulin","Exercise"))
Quant_human_WB_data$Timepoint <- factor(Quant_human_WB_data$Timepoint,levels=c("Pre","Post"))

#mTOR data#
mTOR_data <- Quant_human_WB_data[Quant_human_WB_data$Protein == "pMTOR_S2448",]

mean_df <- mTOR_data %>%
  group_by(Intervention,Timepoint) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

Basal_WB_data <- mTOR_data %>% filter(Timepoint == "Pre")
Insulin_WB_data <- mTOR_data %>% filter(Timepoint == "Post")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(Basal_WB_data, Insulin_WB_data, by = c("ID", "Intervention"), suffix = c("_Pre", "_Post"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("Pre" = "darkgrey", "Post" = "darkblue")  # Colors for points


# Create the plot
p <- ggplot() +
  geom_bar(data = mean_df, aes(x = Timepoint, y = MeanIntensity, fill = Timepoint, color=Timepoint), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  geom_point(data = mTOR_data, aes(x = Timepoint, y = Intensity, color = Timepoint), 
             position = position_dodge(width = 0.8), size = 4) +
  geom_segment(data = line_data, aes(x = factor("Pre"), xend = factor("Post"), y = Intensity_Pre, yend = Intensity_Post, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  facet_grid(~ Intervention) +
  labs(x = "", y = "Fold change from Pre Insulin", title = "p-mTOR S2448") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 3/1)

print(p)

#pAKT data#
pAKT_data <- Quant_human_WB_data[Quant_human_WB_data$Protein == "pAKT_S473",]

mean_df <- pAKT_data %>%
  group_by(Intervention,Timepoint) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

Basal_WB_data <- pAKT_data %>% filter(Timepoint == "Pre")
Insulin_WB_data <- pAKT_data %>% filter(Timepoint == "Post")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(Basal_WB_data, Insulin_WB_data, by = c("ID", "Intervention"), suffix = c("_Pre", "_Post"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("Pre" = "darkgrey", "Post" = "darkblue")  # Colors for points


# Create the plot
p <- ggplot() +
  geom_bar(data = mean_df, aes(x = Timepoint, y = MeanIntensity, fill = Timepoint, color=Timepoint), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  geom_point(data = pAKT_data, aes(x = Timepoint, y = Intensity, color = Timepoint), 
             position = position_dodge(width = 0.8), size = 4) +
  geom_segment(data = line_data, aes(x = factor("Pre"), xend = factor("Post"), y = Intensity_Pre, yend = Intensity_Post, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  facet_grid(~ Intervention) +
  labs(x = "", y = "Fold change from Pre Insulin", title = "p-AKT S473") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 3/1)

print(p)


#pACC data#
pACC_data <- Quant_human_WB_data[Quant_human_WB_data$Protein == "pACC",]

mean_df <- pACC_data %>%
  group_by(Intervention,Timepoint) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

Basal_WB_data <- pACC_data %>% filter(Timepoint == "Pre")
Insulin_WB_data <- pACC_data %>% filter(Timepoint == "Post")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(Basal_WB_data, Insulin_WB_data, by = c("ID", "Intervention"), suffix = c("_Pre", "_Post"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("Pre" = "darkgrey", "Post" = "darkblue")  # Colors for points


# Create the plot
p <- ggplot() +
  geom_bar(data = mean_df, aes(x = Timepoint, y = MeanIntensity, fill = Timepoint, color=Timepoint), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  geom_point(data = pACC_data, aes(x = Timepoint, y = Intensity, color = Timepoint), 
             position = position_dodge(width = 0.8), size = 4) +
  geom_segment(data = line_data, aes(x = factor("Pre"), xend = factor("Post"), y = Intensity_Pre, yend = Intensity_Post, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  facet_grid(~ Intervention) +
  labs(x = "", y = "Fold change from Pre Insulin", title = "p-ACC S221") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 3/1)

print(p)


#pAMPK data#
pAMPK_data <- Quant_human_WB_data[Quant_human_WB_data$Protein == "pAMPK_T172",]

mean_df <- pAMPK_data %>%
  group_by(Intervention,Timepoint) %>%
  summarise(MeanIntensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

Basal_WB_data <- pAMPK_data %>% filter(Timepoint == "Pre")
Insulin_WB_data <- pAMPK_data %>% filter(Timepoint == "Post")

# Combine them to form a new dataframe that will help us draw lines
line_data <- left_join(Basal_WB_data, Insulin_WB_data, by = c("ID", "Intervention"), suffix = c("_Pre", "_Post"))

# ... [previous code remains unchanged up to the plotting section]
point_colors <- c("Pre" = "darkgrey", "Post" = "darkblue")  # Colors for points


# Create the plot
p <- ggplot() +
  geom_bar(data = mean_df, aes(x = Timepoint, y = MeanIntensity, fill = Timepoint, color=Timepoint), 
           stat = "identity", position = position_dodge(width = 0.8)) +
  geom_point(data = pAMPK_data, aes(x = Timepoint, y = Intensity, color = Timepoint), 
             position = position_dodge(width = 0.8), size = 4) +
  geom_segment(data = line_data, aes(x = factor("Pre"), xend = factor("Post"), y = Intensity_Pre, yend = Intensity_Post, group = ID), 
               color = "black", alpha = 0.6, size=0.3) +
  facet_grid(~ Intervention) +
  labs(x = "", y = "Fold change from Pre Insulin", title = "p-AMPK T172") +
  scale_fill_manual(values = c("white","white")) +
  scale_color_manual(values = point_colors) +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),  # Remove x-axis labels
        axis.text = element_text(color = "black")) +
  theme(legend.title = element_blank(), aspect.ratio = 3/1)

print(p)


