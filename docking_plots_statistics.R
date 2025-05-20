###########################################################################################################################################################
###########################################################################################################################################################

# Statistics and plots for AutodockVina docking results, written by: Rodolfo Adrian Arriaga Rivera, refined and corrected by OpenAI,2025

###########################################################################################################################################################
###########################################################################################################################################################

# MIT License
# Copyright (c) 2025 Rodolfo Adrian Arriaga Rivera
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction.
# See the LICENSE file for full license text.

###########################################################################################################################################################
###########################################################################################################################################################

#A preprocess is required to adjust the data for this code, refer to bash scripts by same author: docking_10.sh, extract_vina.sh (this two scripts generate 
#the inputs necessary for this script to work)

# Set working directory
setwd("~/Documents/RNA_Genomics")

#Loading libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# First, prepare the data (that is, splitting columns and convert the file to .csv), then load it into RStudio
# For IFIT1 (won't consider ppp-RNA)
data_cap0_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_cap0/ifit1_cap0_extracted.csv", header=F)
data_cap1_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_cap1/ifit1_cap1_extracted.csv", header=F)
data_cap2_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_cap2/ifit1_cap2_extracted.csv", header=F)
data_fad_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_FAD/ifit1_FAD_extracted.csv", header=F)
data_nad_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_NAD/ifit1_NAD_extracted.csv", header=F)
# data_ppp_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_ppp/ifit1_ppp_extracted.csv", header=F)
# data_Ap4A_ifit1 = read.csv("~/Documents/RNA_Genomics/ifit1_ap4a/ifit1_extracted_columns.txt", header=F, sep = "\t", stringsAsFactors = F)
#This file is read after spliting the columns with a script modification on 26/August/2024 in extrac_vina.sh
# data_m6a_ifit1 = read.csv("~/Documents/RNA_Genomics/ifit1_m6a/ifit1_extracted_columns.txt", sep = "\t", header = F, stringsAsFactors = F)
#This file is read after spliting the columns with a script modification on 26/August/2024 in extrac_vina.sh
#data_TMG_ifit1 = read.csv("~/Documents/RNA_Genomics/ifit1_TMG/ifit1_TMG_extracted_columns.txt", sep = "\t", header = F, stringsAsFactors = F)

# Extract column 2 (energy values)
cap0_energy = select(data_cap0_ifit1, V2)
cap1_energy = select(data_cap1_ifit1, V2)
cap2_energy = select(data_cap2_ifit1, V2)
fad_energy = select(data_fad_ifit1, V2)
nad_energy = select(data_nad_ifit1, V2)
# ppp_energy = select(data_ppp_ifit1, V2)
# ap4a_energy = select(data_Ap4A_ifit1, V2)
# m6a_energy = select(data_m6a_ifit1, V2)
# TMG_energy = select(data_TMG_ifit1, V2)

#Let's create a column named source, so we can know where the data comes from
cap0_energy$source = "cap0"
cap1_energy$source = "cap1"
cap2_energy$source = "cap2"
fad_energy$source = "FAD"
nad_energy$source = "NAD"
# ppp_energy$source = "PPP-RNA"
# ap4a_energy$source = "Ap4A"
# m6a_energy$source= "m6a"
# TMG_energy$source="TMG"

combined_energies = rbind(cap0_energy, cap1_energy, cap2_energy, fad_energy, nad_energy)

# Specifying the desired order for the factors
caps_ordered_ifit1 = c("cap0", "cap1", "cap2", "NAD", "FAD")
combined_energies$source = factor(combined_energies$source, levels = caps_ordered_ifit1)
caps_ordered_ifit1 #checking order of caps

# Create the plot
ggplot(combined_energies, aes(x = source, y = V2)) +
  geom_violin() + 
  labs(title = "Violin Plot of IFIT1 binding affinities",
                       x = "Cap",
                       y = "Binding affinity (kcal/mol)") +
  geom_boxplot(width=0.1) +
  scale_color_grey() + theme_classic() +
  stat_summary(fun.y=mean, geom="point", size=1, color="orange2") 

# Extracting the average from the V2 column
average_data_ifit1 = combined_energies %>%
  group_by(source) %>%
  summarize(mean_V2 = mean(V2))
print(average_data_ifit1)

# Extracting the median from V2 column
median_data_ifit1 = combined_energies %>%
  group_by(source) %>%
  summarize(median_V2 = median(V2))
print(median_data_ifit1)

# For IFIT5 repeat the same process (will only consider NAD, FAD and ppp-RNA)
#data_cap0_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-cap0/ifit5_cap0_extracted.csv", header = F)
#data_cap1_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-cap1/ifit5_cap1_extracted.csv", header = F)
#data_cap2_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-cap2/ifit5_cap2_extracted.csv", header = F)
data_fad_ifit5 = read.csv("~/Documents/RNA_Genomics/docking_ifit5_FAD/ifit5_FAD_extracted.csv", header = F)
data_nad_ifit5 = read.csv("~/Documents/RNA_Genomics/docking_ifit5_NAD/ifit5_NAD_extracted.csv", header = F)
data_ppp_ifit5 = read.csv("~/Documents/RNA_Genomics/docking_ifit5_ppp/ifit5_ppp_extracted.csv", header = F)
data_ppp_mg_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5_mg_ppp_docking/ifit5_extracted_columns.txt", header = F, sep = "\t", stringsAsFactors = F)
#data_Ap4A_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-ap4a/ifit5_ap4a_extracted.csv", header=F)

# Selecting data from energy column V2 , and rename them for consistency (rbind requires this) instead of having V2 and V4
#cap0_ifit5_energy = select(data_cap0_ifit5, V2) 
#cap1_ifit5_energy = select(data_cap1_ifit5, V2)
#cap2_ifit5_energy = select(data_cap2_ifit5, V2)
fad_ifit5_energy = select(data_fad_ifit5, V2) %>% select(V2) %>% rename(energy = V2)
nad_ifit5_energy = select(data_nad_ifit5, V2) %>% select(V2) %>% rename(energy = V2)
ppp_ifit5_energy = select(data_ppp_ifit5, V2) %>% select(V2) %>% rename(energy = V2)
ppp_ifit5_mg_energy = select(data_ppp_mg_ifit5, V2) %>% select(V2) %>% rename(energy = V2)
#ap4a_ifit5_energy = select(data_cap0_ifit1, V2)

#Let's create a column that's called "source" so we can know where the data comes from
#cap0_ifit5_energy$source = "cap0"
#cap1_ifit5_energy$source = "cap1"
#cap2_ifit5_energy$source = "cap2"
fad_ifit5_energy$source = "FAD"
nad_ifit5_energy$source = "NAD"
ppp_ifit5_mg_energy$source = "PPP-RNA-Mg"
#ap4a_ifit5_energy$source = "Ap4A"

# Merging the energies and ordering caps by names
caps_ordered_ifit5 = c("FAD", "NAD", "PPP-RNA-Mg")
combined_energies_ifit5 = rbind(fad_ifit5_energy, nad_ifit5_energy, ppp_ifit5_mg_energy)
combined_energies_ifit5$source = factor(combined_energies_ifit5$source, levels = caps_ordered_ifit5)
caps_ordered_ifit5 #this is to check the order of the caps

# Let's create a violin plot
# A violin plot is used to represent the distribution density of the data at different values, highlighting areas where data points are more concentrated. According to R CHARTS by R CODER 2024: https://r-charts.com/distribution/violin-plot-points-ggplot2/#:~:text=If%20you%20want%20to%20create,the%20size%20of%20the%20points.&text=The%20default%20dots%20are%20of,a%20fill%20color%20inside%20geom_dotplot%20.
ggplot(combined_energies_ifit5, aes(x = source, y = energy)) +
  geom_violin() + 
  labs(title = "Violin Plot of IFIT5 binding affinities",
       x = "IFIT5",
       y = "Binding affinity (kcal/mol)") +
  geom_boxplot(width=0.1) +
  scale_color_grey() + theme_classic() +
  stat_summary(fun.y=mean, geom="point", size=1, color="blue3") 

# Extracting the average from the energy column
average_data_ifit5 = combined_energies_ifit5 %>%
  group_by(source) %>%
  summarize(mean_energy = mean(energy))
print(average_data_ifit5)

# Extracting the median from V2 column
median_data_ifit5 = combined_energies_ifit5 %>%
  group_by(source) %>%
  summarize(median_energy = median(energy))
print(median_data_ifit5)

# Combining both protein binding affinities into one violin plot
# Extracting energy values and adding source and protein columns for IFIT1
cap0_ifit1_energy = select(data_cap0_ifit1, V2) %>% mutate(source = "cap0", protein = "IFIT1")
cap1_ifit1_energy = select(data_cap1_ifit1, V2) %>% mutate(source = "cap1", protein = "IFIT1")
cap2_ifit1_energy = select(data_cap2_ifit1, V2) %>% mutate(source = "cap2", protein = "IFIT1")
fad_ifit1_energy = select(data_fad_ifit1, V2) %>% mutate(source = "FAD", protein = "IFIT1")
nad_ifit1_energy = select(data_nad_ifit1, V2) %>% mutate(source = "NAD", protein = "IFIT1")
#ppp_ifit1_energy = select(data_ppp_ifit1, V2) %>% mutate(source = "PPP-RNA", protein = "IFIT1")
#ap4a_ifit1_energy = select(data_Ap4A_ifit1, V2) %>% mutate(source = "Ap4A", protein = "IFIT1")
#m6a_ifit1_energy = select(data_m6a_ifit1, V2) %>% mutate(source = "m6a", protein = "IFIT1")

# Extracting energy values and adding source and protein columns for IFIT5
#cap0_ifit5_energy = select(data_cap0_ifit5, V2) %>% mutate(source = "cap0", protein = "IFIT5")
#cap1_ifit5_energy = select(data_cap1_ifit5, V2) %>% mutate(source = "cap1", protein = "IFIT5")
#cap2_ifit5_energy = select(data_cap2_ifit5, V2) %>% mutate(source = "cap2", protein = "IFIT5")
fad_ifit5_energy = select(data_fad_ifit5, V2) %>% mutate(source = "FAD", protein = "IFIT5")
nad_ifit5_energy = select(data_nad_ifit5, V2) %>% mutate(source = "NAD", protein = "IFIT5")
nad_ifit5_energy
ppp_ifit5_mg_energy = select(data_ppp_mg_ifit5, V2) %>% mutate(source = "PPP-MG-RNA", protein = "IFIT5")
ppp_ifit5_mg_energy
#ap4a_ifit5_energy = select(data_Ap4A_ifit5, V2) %>% mutate(source = "Ap4A", protein = "IFIT5")

# Setting caps levels in a string
caps_ifit1 = c("cap0", "cap1", "cap2", "FAD", "NAD")
caps_ifit5 = c("PPP-MG-RNA", "NAD", "FAD")

# Combining each protein's energies
combined_energies_ifit1 = bind_rows(cap0_ifit1_energy, cap1_ifit1_energy, cap2_ifit1_energy, fad_ifit1_energy, nad_ifit1_energy)
combined_energies_ifit1$source = factor(combined_energies_ifit1$source, levels = caps_ifit1)
combined_energies_ifit5 = bind_rows(ppp_ifit5_mg_energy, nad_ifit5_energy, fad_ifit5_energy)
combined_energies_ifit5$source = factor(combined_energies_ifit5$source, levels = caps_ifit5)

# Combining both protein energies
combined_energies = bind_rows(combined_energies_ifit1, combined_energies_ifit5)

# Violin plot splitted by proteins
ggplot(combined_energies, aes(x = source, y = V2)) +
  geom_violin(position = position_dodge(1)) + 
  labs(title = "Violin Plot of IFIT1 and IFIT5 binding affinities",
       x = "Cap Type",
       y = "Binding affinity (kcal/mol)") +
  geom_boxplot(width = 0.1, position = position_dodge(1)) +
  scale_color_grey() + 
  theme_classic() +
  stat_summary(fun = mean, geom = "point", size = 1, color = "orange2", position = position_dodge(0)) +
  facet_wrap(~ protein, scales = "free_x")

# Violin plot side-by-side coloured
ggplot(combined_energies, aes(x = source, y = V2, fill = protein)) +
  geom_violin(position = position_dodge(width = 0.9)) + 
  labs(title = "Violin Plot of IFIT1 and IFIT5 binding affinities",
       x = "Cap Type",
       y = "Binding affinity (kcal/mol)") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("IFIT1" = "lightblue2", "IFIT5" = "lightgreen")) + 
  theme_classic() +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red4", position = position_dodge(width = 0.9)) +
  facet_wrap(~ protein, scales = "free_x")

# Violin plot side-by-side
ggplot(combined_energies, aes(x = source, y = V2, fill = protein)) +
  geom_violin(position = position_dodge(width = 0.9)) + 
  labs(title = "Violin Plot of IFIT1 and IFIT5 binding affinities",
       x = "Cap Type",
       y = "Binding affinity (kcal/mol)") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("IFIT1" = "lightblue2", "IFIT5" = "lightgreen")) + 
  theme_classic() +
  stat_summary(fun = mean, geom = "point", size = 1, color = "red4", position = position_dodge(width = 0.9)) 

################################################################################################################################################################################
# Now, the next question to be answered should be: is my data parametric or non-parametric?
# Visual and statistical methods exist for this. Parametric refers to normal distribution and equal variances.
# Statistical tests such as Shapiro-Wilk and Kolmogrov-Smirnov (sample size >50) test the null hypothesis that the data is normally distributed. 
# Visual methods include Q-Q Plots and Histograms. Normality Tests: If the p-value from the Shapiro-Wilk test is less than 0.05, the data significantly deviates from normality.
################################################################################################################################################################################

# Installing libraries
install.packages("car")
library(car)
install.packages("nortest")
library(nortest)

#Making lists of energy data grouped for both IFIT's
groups_ifit1 = list(cap0_energy, cap1_energy, cap2_energy, fad_energy, nad_energy, ppp_energy)
group_names_ifit1 = c("PPP-RNA", "cap0", "cap1", "cap2", "FAD", "FAD")

groups_ifit5 = list(fad_ifit5_energy, nad_ifit5_energy, ppp_ifit5_energy)
group_names_ifit5 = c("PPP-RNA", "NAD", "FAD")

#Performing Shapiro-Wilk for IFIT1
cat("Shapiro-Wilk Test Results for IFIT1:\n")
for(i in 1:length(groups_ifit1)) {
  cat(group_names_ifit1[i], ":\n")
  print(shapiro.test(groups_ifit1[[i]]$V2))
  cat("\n")
}

# Defining a function for creating Q-Q plots with a main title for IFIT1, calculating its dimension to fit adequately with a rows and cols object
# In a Q-Q plot, each point represents a comparison between quantiles of two datasets. 
# The x-axis generally represents theoretical quantiles (e.g., from a normal distribution), and the y-axis represents the quantiles from your dataset.

create_qq_plots = function(groups, group_names, main_title) {
  n = length(groups)
  rows = ceiling(sqrt(n))
  cols = ceiling(n / rows)
# Adjusting layout to fit outer and inner margins in a single plot sheet
  layout(matrix(1:(rows * cols), nrow = rows, byrow = T))
  par(oma = c(0, 0, 2, 0)) # Setting outer margins to make space for the main title
  par(mar = c(4, 4, 2, 1)) # Setting inner margins for individual plots
  
  for(i in 1:n) {
    if (i > rows * cols) {
      break
    }
    qqnorm(groups[[i]]$V2, main = paste("Q-Q Plot for", group_names[i]))
    qqline(groups[[i]]$V2)
  }
  mtext(main_title, outer = T, cex = 1.5)
  par(mfrow = c(1, 1)) # Reseting the layout to default to avoid affecting other plots
}
groups_ifit1 = list(cap0_energy, cap1_energy, cap2_energy, fad_energy, nad_energy, ppp_energy, ap4a_ifit1_energy)
group_names_ifit1 = c("PPP-RNA", "cap0", "cap1", "cap2", "NAD", "FAD", "Ap4A")

create_qq_plots(groups_ifit1, group_names_ifit1, "IFIT1")

#Performing Shapiro-Wilk for IFIT5
cat("Shapiro-Wilk Test Results for IFIT5:\n")
for(i in 1:length(groups_ifit5)) {
  cat(group_names_ifit5[i], ":\n")
  print(shapiro.test(groups_ifit5[[i]]$V2))
  cat("\n")
}

# Defining a function for creating Q-Q plots with a main title for IFIT5 calculating its dimension to fit adequately with a rows and cols object
# Define a function for creating Q-Q plots with a main title
create_qq_plots = function(groups, group_names, main_title) {
  n = length(groups)
  rows = ceiling(sqrt(n))
  cols = ceiling(n / rows)
# Adjusting layout to fit outer and inner margins in a single plot sheet
  layout(matrix(1:(rows * cols), nrow = rows, byrow = T))
  par(oma = c(0, 0, 2, 0)) # Set outer margins to make space for the main title
  par(mar = c(4, 4, 2, 1)) # Set inner margins for individual plots
  
  for(i in 1:n) {
    if (i > rows * cols) {
      break
    }
    qqnorm(groups[[i]]$V2, main = paste("Q-Q Plot for", group_names[i]))
    qqline(groups[[i]]$V2)
  }
  mtext(main_title, outer = T, cex = 1.5)
  par(mfrow = c(1, 1)) # Reset the layout to default to avoid affecting other plots
}

groups_ifit5 = list(cap0_ifit5_energy, cap1_ifit5_energy, cap2_ifit5_energy, fad_ifit5_energy, nad_ifit5_energy, ppp_ifit5_energy, ap4a_ifit5_energy)
group_names_ifit5 = c("PPP-RNA", "cap0", "cap1", "cap2", "NAD", "FAD", "Ap4A")

create_qq_plots(groups_ifit5, group_names_ifit5, "IFIT5")

# Testing for variance equality (e.g.: variances are significantly different from each other)
# Levene's Test for IFIT1
cat("Levene's Test for Homogeneity of Variance for IFIT1:\n")
print(leveneTest(V2 ~ source, data = combined_energies_ifit1))

# Levene's Test for IFIT5
cat("Levene's Test for Homogeneity of Variance for IFIT5:\n")
print(leveneTest(V2 ~ source, data = combined_energies_ifit5))

######################################################################################################################################################################################################################
# Selecting the best dock by CompositeScore (considering binding affinity and both RMSD u.b & l.b values)
# Rank the results based on a composite score or use more sophisticated methods like weighted scoring.
# It is plausible to incorporate both RMSD lower bound (RMSD l.b.) and RMSD upper bound (RMSD u.b.) into a composite score for evaluating docking results, along with binding affinities. 
# This approach can provide a more nuanced assessment of docking accuracy by considering both the core stability and the overall conformational fit of the ligand, as well as the strength of the binding interaction.
######################################################################################################################################################################################################################

# In this case, a single value is desired for each cap (just to remember the data)

data_cap0_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_cap0/ifit1_cap0_extracted.csv", header=F)
data_cap1_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_cap1/ifit1_cap1_extracted.csv", header=F)
data_cap2_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_cap2/ifit1_cap2_extracted.csv", header=F)
data_fad_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_FAD/ifit1_FAD_extracted.csv", header=F)
data_nad_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_NAD/ifit1_NAD_extracted.csv", header=F)
#data_ppp_ifit1 = read.csv("~/Documents/RNA_Genomics/docking_ifit1_ppp/ifit1_ppp_extracted.csv", header=F)
#data_Ap4A_ifit1 = read.csv("~/Documents/RNA_Genomics/ifit1_ap4a/ifit1_extracted_columns.txt", header=F,  sep = "\t", stringsAsFactors = F)
#data_m6a_ifit1 = read.csv("~/Documents/RNA_Genomics/ifit1_m6a/ifit1_extracted_columns.txt", sep = "\t", header = F, stringsAsFactors = F)
#data_TMG_ifit1 = read.csv("~/Documents/RNA_Genomics/ifit1_TMG/ifit1_TMG_extracted_columns.txt", sep = "\t", header = F, stringsAsFactors = F)

# Extract column 2 (energy values)
cap0_energy = select(data_cap0_ifit1, V2)
cap1_energy = select(data_cap1_ifit1, V2)
cap2_energy = select(data_cap2_ifit1, V2)
fad_energy = select(data_fad_ifit1, V2)
nad_energy = select(data_nad_ifit1, V2)
#ppp_energy = select(data_ppp_ifit1, V2)
#ap4a_energy = select(data_Ap4A_ifit1, V2)
#m6a_energy = select(data_m6a_ifit1, V2)
#TMG_energy = select(data_TMG_ifit1, V2)

# Extracting rmsd u.b values (column V4)
cap0_ifit1_rmsd_ub = select(data_cap0_ifit1, V4)
cap1_ifit1_rmsd_ub = select(data_cap1_ifit1, V4)
cap2_ifit1_rmsd_ub = select(data_cap2_ifit1, V4)
fad_ifit1_rmsd_ub = select(data_fad_ifit1, V4)
nad_ifit1_rmsd_ub = select(data_nad_ifit1, V4)
#ppp_ifit1_rmsd_ub = select(data_ppp_ifit1, V4)
#ap4a_ifit1_rmsd_ub = select(data_Ap4A_ifit1, V4)
#m6a_ifit1_rmsd_ub = select(data_m6a_ifit1, V4)
#TMG_ifit1_rmsd_ub = select(data_TMG_ifit1, V4)

# Extracting rmsd l.b values (column V4)
cap0_ifit1_rmsd_lb = select(data_cap0_ifit1, V3)
cap1_ifit1_rmsd_lb = select(data_cap1_ifit1, V3)
cap2_ifit1_rmsd_lb = select(data_cap2_ifit1, V3)
fad_ifit1_rmsd_lb = select(data_fad_ifit1, V3)
nad_ifit1_rmsd_lb = select(data_nad_ifit1, V3)
#ppp_ifit1_rmsd_lb = select(data_ppp_ifit1, V3)
#ap4a_ifit1_rmsd_lb = select(data_Ap4A_ifit1, V3)
#m6a_ifit1_rmsd_lb = select(data_m6a_ifit1, V3)
#TMG_ifit1_rmsd_lb = select(data_TMG_ifit1, V3)

# Normalizing data for each cap
# Normalized values of rmsd and binding affinities cap0-IFIT1
data_cap0_ifit1 = data_cap0_ifit1 %>%
  mutate(NormalizedBindingAffinity_ifit1_cap0 = (cap0_energy - min(cap0_energy)) / (max(cap0_energy) - min(cap0_energy)),
         NormalizedRMSD_lb_ifit1_cap0 = (cap0_ifit1_rmsd_lb - min(cap0_ifit1_rmsd_lb)) / (max(cap0_ifit1_rmsd_lb) - min(cap0_ifit1_rmsd_lb)),
         NormalizedRMSD_ub_ifit1_cap0 = (cap0_ifit1_rmsd_ub - min(cap0_ifit1_rmsd_ub)) / (max(cap0_ifit1_rmsd_ub) - min(cap0_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities cap1-IFIT1
data_cap1_ifit1 = data_cap1_ifit1 %>%
  mutate(NormalizedBindingAffinity_ifit1_cap1 = (cap1_energy - min(cap1_energy)) / (max(cap1_energy) - min(cap1_energy)),
         NormalizedRMSD_lb_ifit1_cap1 = (cap1_ifit1_rmsd_lb - min(cap1_ifit1_rmsd_lb)) / (max(cap1_ifit1_rmsd_lb) - min(cap1_ifit1_rmsd_lb)),
         NormalizedRMSD_ub_ifit1_cap1 = (cap1_ifit1_rmsd_ub - min(cap1_ifit1_rmsd_ub)) / (max(cap1_ifit1_rmsd_ub) - min(cap1_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities cap2-IFIT1
data_cap2_ifit1 = data_cap2_ifit1 %>%
  mutate(NormalizedBindingAffinity_ifit1_cap2 = (cap2_energy - min(cap2_energy)) / (max(cap2_energy) - min(cap2_energy)),
         NormalizedRMSD_lb_ifit1_cap2 = (cap2_ifit1_rmsd_lb - min(cap2_ifit1_rmsd_lb)) / (max(cap2_ifit1_rmsd_lb) - min(cap2_ifit1_rmsd_lb)),
         NormalizedRMSD_ub_ifit1_cap2 = (cap2_ifit1_rmsd_ub - min(cap2_ifit1_rmsd_ub)) / (max(cap2_ifit1_rmsd_ub) - min(cap2_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities FAD-IFIT1
data_fad_ifit1 = data_fad_ifit1 %>%
  mutate(NormalizedBindingAffinity_ifit1_fad = (fad_energy - min(fad_energy)) / (max(fad_energy) - min(fad_energy)),
         NormalizedRMSD_lb_ifit1_fad = (fad_ifit1_rmsd_lb - min(fad_ifit1_rmsd_lb)) / (max(fad_ifit1_rmsd_lb) - min(fad_ifit1_rmsd_lb)),
         NormalizedRMSD_ub_ifit1_fad = (fad_ifit1_rmsd_ub - min(fad_ifit1_rmsd_ub)) / (max(fad_ifit1_rmsd_ub) - min(fad_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities NAD-IFIT1
data_nad_ifit1 = data_nad_ifit1 %>%
  mutate(NormalizedBindingAffinity_ifit1_nad = (nad_energy - min(nad_energy)) / (max(nad_energy) - min(nad_energy)),
         NormalizedRMSD_lb_ifit1_nad = (nad_ifit1_rmsd_lb - min(nad_ifit1_rmsd_lb)) / (max(nad_ifit1_rmsd_lb) - min(nad_ifit1_rmsd_lb)),
         NormalizedRMSD_ub_ifit1_nad = (nad_ifit1_rmsd_ub - min(nad_ifit1_rmsd_ub)) / (max(nad_ifit1_rmsd_ub) - min(nad_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities ppp-IFIT1
#data_ppp_ifit1 = data_ppp_ifit1 %>%
  #mutate(NormalizedBindingAffinity_ifit1_ppp = (ppp_energy - min(ppp_energy)) / (max(ppp_energy) - min(ppp_energy)),
         #NormalizedRMSD_lb_ifit1_ppp = (ppp_ifit1_rmsd_lb - min(ppp_ifit1_rmsd_lb)) / (max(ppp_ifit1_rmsd_lb) - min(ppp_ifit1_rmsd_lb)),
         #NormalizedRMSD_ub_ifit1_ppp = (ppp_ifit1_rmsd_ub - min(ppp_ifit1_rmsd_ub)) / (max(ppp_ifit1_rmsd_ub) - min(ppp_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities ap4a-IFIT1
#data_Ap4A_ifit1 = data_Ap4A_ifit1 %>%
 # mutate(NormalizedBindingAffinity_ifit1_ap4a = (ap4a_energy - min(ap4a_energy)) / (max(ap4a_energy) - min(ap4a_energy)),
  #       NormalizedRMSD_lb_ifit1_ap4a = (ap4a_ifit1_rmsd_lb - min(ap4a_ifit1_rmsd_lb)) / (max(ap4a_ifit1_rmsd_lb) - min(ap4a_ifit1_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit1_ap4a = (ap4a_ifit1_rmsd_ub - min(ap4a_ifit1_rmsd_ub)) / (max(ap4a_ifit1_rmsd_ub) - min(ap4a_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities m6a-IFIT1
#data_m6a_ifit1 = data_m6a_ifit1 %>%
 # mutate(NormalizedBindingAffinity_ifit1_m6a = (m6a_energy - min(m6a_energy)) / (max(m6a_energy) - min(m6a_energy)),
  #       NormalizedRMSD_lb_ifit1_m6a = (m6a_ifit1_rmsd_lb - min(m6a_ifit1_rmsd_lb)) / (max(m6a_ifit1_rmsd_lb) - min(m6a_ifit1_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit1_m6a = (m6a_ifit1_rmsd_ub - min(m6a_ifit1_rmsd_ub)) / (max(m6a_ifit1_rmsd_ub) - min(m6a_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities
# Normalized values of rmsd and binding affinities TMG-IFIT1
#data_TMG_ifit1 = data_TMG_ifit1 %>%
 # mutate(NormalizedBindingAffinity_ifit1_TMG = (TMG_energy - min(TMG_energy)) / (max(TMG_energy) - min(TMG_energy)),
  #       NormalizedRMSD_lb_ifit1_TMG = (TMG_ifit1_rmsd_lb - min(TMG_ifit1_rmsd_lb)) / (max(TMG_ifit1_rmsd_lb) - min(TMG_ifit1_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit1_TMG = (TMG_ifit1_rmsd_ub - min(TMG_ifit1_rmsd_ub)) / (max(TMG_ifit1_rmsd_ub) - min(TMG_ifit1_rmsd_ub))) # normalized values of rmsd and binding affinities

# Assigning weights
weight_binding_affinity = 0.4
weight_lb = 0.3
weight_ub = 0.3

# Calculating CompositeScore for each cap
# CompositeScore for cap0-IFIT1
data_cap0_ifit1 = data_cap0_ifit1 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_cap0 +
           weight_lb * NormalizedRMSD_lb_ifit1_cap0 +
           weight_ub * NormalizedRMSD_ub_ifit1_cap0)
head(data_cap0_ifit1)
# CompositeScore for cap1-IFIT1
data_cap1_ifit1 = data_cap1_ifit1 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_cap1 +
           weight_lb * NormalizedRMSD_lb_ifit1_cap1 +
           weight_ub * NormalizedRMSD_ub_ifit1_cap1)
head(data_cap1_ifit1)
# CompositeScore for cap2-IFIT1
data_cap2_ifit1 = data_cap2_ifit1 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_cap2 +
           weight_lb * NormalizedRMSD_lb_ifit1_cap2 +
           weight_ub * NormalizedRMSD_ub_ifit1_cap2)
head(data_cap2_ifit1)
# CompositeScore for fad-IFIT1
data_fad_ifit1 = data_fad_ifit1 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_fad +
           weight_lb * NormalizedRMSD_lb_ifit1_fad +
           weight_ub * NormalizedRMSD_ub_ifit1_fad)
head(data_fad_ifit1)
# CompositeScore for nad-IFIT1
data_nad_ifit1 = data_nad_ifit1 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_nad +
           weight_lb * NormalizedRMSD_lb_ifit1_nad +
           weight_ub * NormalizedRMSD_ub_ifit1_nad)
head(data_nad_ifit1)
# CompositeScore for ppp-IFIT1
#data_ppp_ifit1 = data_ppp_ifit1 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_ppp +
  #         weight_lb * NormalizedRMSD_lb_ifit1_ppp +
   #        weight_ub * NormalizedRMSD_ub_ifit1_ppp)
#head(data_ppp_ifit1)
# CompositeScore for Ap4A-IFIT1
#data_Ap4A_ifit1 = data_Ap4A_ifit1 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_ap4a +
  #         weight_lb * NormalizedRMSD_lb_ifit1_ap4a +
   #        weight_ub * NormalizedRMSD_ub_ifit1_ap4a)
#head(data_Ap4A_ifit1)
# CompositeScore for m6a-IFIT1
#data_m6a_ifit1 = data_m6a_ifit1 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_m6a +
  #         weight_lb * NormalizedRMSD_lb_ifit1_m6a +
   #        weight_ub * NormalizedRMSD_ub_ifit1_m6a)
#head(data_m6a_ifit1)
#Composite score for TMG-IFIT1
#data_TMG_ifit1 = data_TMG_ifit1 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit1_TMG +
  #         weight_lb * NormalizedRMSD_lb_ifit1_TMG +
   #        weight_ub * NormalizedRMSD_ub_ifit1_TMG)
#head(data_TMG_ifit1)

#Picking best dock for each cap
# A result with RMSD = 0 usually indicates that the pose is identical to the reference pose, meaning there was no movement or change in the ligand's position during docking.
# Best dock cap0-IFIT1
ranked_docking_results_ifit1_cap0 = data_cap0_ifit1 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit1_cap0)
filtered_results_ifit1_cap0 = ranked_docking_results_ifit1_cap0 %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit1_cap0) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit1_cap0, file = "filtered_results_ifit1_cap0.csv", row.names = F)

# Best dock cap1-IFIT1
ranked_docking_results_ifit1_cap1 = data_cap1_ifit1 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit1_cap1)
filtered_results_ifit1_cap1 = ranked_docking_results_ifit1_cap1 %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit1_cap1) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit1_cap1, file = "filtered_results_ifit1_cap1.csv")

# Best dock cap2-IFIT1
ranked_docking_results_ifit1_cap2 = data_cap2_ifit1 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit1_cap2)
filtered_results_ifit1_cap2 = ranked_docking_results_ifit1_cap2 %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit1_cap2) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit1_cap2, file = "filtered_results_ifit1_cap2.csv")

# Best dock FAD-IFIT1
ranked_docking_results_ifit1_fad = data_fad_ifit1 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit1_fad)
filtered_results_ifit1_fad = ranked_docking_results_ifit1_fad %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit1_fad) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit1_fad, file = "filtered_results_ifit1_fad.csv")

# Best dock NAD-IFIT1
ranked_docking_results_ifit1_nad = data_nad_ifit1 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit1_nad)
filtered_results_ifit1_nad = ranked_docking_results_ifit1_nad %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit1_nad) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit1_nad, file = "filtered_results_ifit1_nad.csv")

# Best dock ppp-IFIT1
#ranked_docking_results_ifit1_ppp = data_ppp_ifit1 %>%
  #arrange(CompositeScore)
#View(ranked_docking_results_ifit1_ppp)
#filtered_results_ifit1_ppp = ranked_docking_results_ifit1_ppp %>%
  #filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit1_ppp) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit1_ppp, file = "filtered_results_ifit1_ppp.csv")

# Best dock ap4a-IFIT1
#ranked_docking_results_ifit1_ap4a = data_Ap4A_ifit1 %>%
 # arrange(CompositeScore)
#View(ranked_docking_results_ifit1_ap4a)
#filtered_results_ifit1_ap4a = ranked_docking_results_ifit1_ap4a %>%
 #filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit1_ap4a) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit1_ap4a, file = "filtered_results_ifit1_ap4a.csv")

# Best dock m6a-IFIT1
#ranked_docking_results_ifit1_m6a = data_m6a_ifit1 %>%
 # arrange(CompositeScore)
#View(ranked_docking_results_ifit1_m6a)
#filtered_results_ifit1_m6a = ranked_docking_results_ifit1_m6a %>%
#  filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit1_m6a) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit1_m6a, file = "filtered_results_ifit1_m6a.csv", row.names = F)

#Best dock TMG_ifit1
#ranked_docking_results_ifit1_TMG = data_TMG_ifit1 %>%
  #arrange(CompositeScore)
#View(ranked_docking_results_ifit1_TMG)
#filtered_results_ifit1_TMG = ranked_docking_results_ifit1_TMG %>%
  #filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit1_TMG) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit1_TMG, file = "filtered_results_ifit1_TMG.csv", row.names = F)

######################################## Now, same process applied to IFIT5 dockings #################################################

# Data frames for IFIT5
#data_cap0_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-cap0/ifit5_cap0_extracted.csv", header = F)
#data_cap1_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-cap1/ifit5_cap1_extracted.csv", header = F)
#data_cap2_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-cap2/ifit5_cap2_extracted.csv", header = F)
data_fad_ifit5 = read.csv("~/Documents/RNA_Genomics/docking_ifit5_FAD/ifit5_FAD_extracted.csv", header = F)
data_nad_ifit5 = read.csv("~/Documents/RNA_Genomics/docking_ifit5_NAD/ifit5_NAD_extracted.csv", header = F)
#data_ppp_ifit5 = read.csv("~/Documents/RNA_Genomics/docking_ifit5_ppp/ifit5_ppp_extracted.csv", header = F)
data_ppp_mg_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5_mg_ppp_docking/ifit5_extracted_columns.txt", header = F, sep = "\t", stringsAsFactors = F)
#data_Ap4A_ifit5 = read.csv("~/Documents/RNA_Genomics/ifit5-ap4a/ifit5_extracted_columns.txt", sep = "\t", header = F, stringsAsFactors = F)

# Extracting energy data from column V2
#cap0_ifit5_energy = select(data_cap0_ifit5, V2)
#cap1_ifit5_energy = select(data_cap1_ifit5, V2)
#cap2_ifit5_energy = select(data_cap2_ifit5, V2)
fad_ifit5_energy = select(data_fad_ifit5, V2)
nad_ifit5_energy = select(data_nad_ifit5, V2)
#ppp_ifit5_energy = select(data_ppp_ifit5, V2)
ppp_mg_energy = select(data_ppp_mg_ifit5, V2)
#ap4a_ifit5_energy = select(data_cap0_ifit1, V2)

# Extracting rmsd u.b values (column V4)
#cap0_ifit5_rmsd_ub = select(data_cap0_ifit5, V4)
#cap1_ifit5_rmsd_ub = select(data_cap1_ifit5, V4)
#cap2_ifit5_rmsd_ub = select(data_cap2_ifit5, V4)
fad_ifit5_rmsd_ub = select(data_fad_ifit5, V4)
nad_ifit5_rmsd_ub = select(data_nad_ifit5, V4)
#ppp_ifit5_rmsd_ub = select(data_ppp_ifit5, V4)
ppp_mg_ifit5_rmsd_ub = select(data_ppp_mg_ifit5, V4)
#ap4a_ifit5_rmsd_ub = select(data_Ap4A_ifit5, V4)

# Extracting rmsd l.b values (column V4)
#cap0_ifit5_rmsd_lb = select(data_cap0_ifit5, V3)
#cap1_ifit5_rmsd_lb = select(data_cap1_ifit5, V3)
#cap2_ifit5_rmsd_lb = select(data_cap2_ifit5, V3)
fad_ifit5_rmsd_lb = select(data_fad_ifit5, V3)
nad_ifit5_rmsd_lb = select(data_nad_ifit5, V3)
#ppp_ifit5_rmsd_lb = select(data_ppp_ifit5, V3)
ppp_mg_ifit5_rmsd_lb = select(data_ppp_mg_ifit5, V3)
#ap4a_ifit5_rmsd_lb = select(data_Ap4A_ifit5, V3)

# Normalizing data for each cap
# Normalized values of rmsd and binding affinities cap0-IFIT5
#data_cap0_ifit5 = data_cap0_ifit5 %>%
 # mutate(NormalizedBindingAffinity_ifit5_cap0 = (cap0_ifit5_energy - min(cap0_ifit5_energy)) / (max(cap0_ifit5_energy) - min(cap0_ifit5_energy)),
  #       NormalizedRMSD_lb_ifit5_cap0 = (cap0_ifit5_rmsd_lb - min(cap0_ifit5_rmsd_lb)) / (max(cap0_ifit5_rmsd_lb) - min(cap0_ifit5_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit5_cap0 = (cap0_ifit5_rmsd_ub - min(cap0_ifit5_rmsd_ub)) / (max(cap0_ifit5_rmsd_ub) - min(cap0_ifit5_rmsd_ub))) 
# Normalized values of rmsd and binding affinities cap1-IFIT5
#data_cap1_ifit5 = data_cap1_ifit5 %>%
 # mutate(NormalizedBindingAffinity_ifit5_cap1 = (cap1_ifit5_energy - min(cap1_ifit5_energy)) / (max(cap1_ifit5_energy) - min(cap1_ifit5_energy)),
  #       NormalizedRMSD_lb_ifit5_cap1 = (cap1_ifit5_rmsd_lb - min(cap1_ifit5_rmsd_lb)) / (max(cap1_ifit5_rmsd_lb) - min(cap1_ifit5_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit5_cap1 = (cap1_ifit5_rmsd_ub - min(cap1_ifit5_rmsd_ub)) / (max(cap1_ifit5_rmsd_ub) - min(cap1_ifit5_rmsd_ub))) 
# Normalized values of rmsd and binding affinities cap2-IFIT5
#data_cap2_ifit5 = data_cap2_ifit5 %>%
 # mutate(NormalizedBindingAffinity_ifit5_cap2 = (cap2_ifit5_energy - min(cap2_ifit5_energy)) / (max(cap2_ifit5_energy) - min(cap2_ifit5_energy)),
  #       NormalizedRMSD_lb_ifit5_cap2 = (cap2_ifit5_rmsd_lb - min(cap2_ifit5_rmsd_lb)) / (max(cap2_ifit5_rmsd_lb) - min(cap2_ifit5_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit5_cap2 = (cap2_ifit5_rmsd_ub - min(cap2_ifit5_rmsd_ub)) / (max(cap2_ifit5_rmsd_ub) - min(cap2_ifit5_rmsd_ub))) 
# Normalized values of rmsd and binding affinities FAD-IFIT5
data_fad_ifit5 = data_fad_ifit5 %>%
  mutate(NormalizedBindingAffinity_ifit5_fad = (fad_ifit5_energy - min(fad_ifit5_energy)) / (max(fad_ifit5_energy) - min(fad_ifit5_energy)),
         NormalizedRMSD_lb_ifit5_fad = (fad_ifit5_rmsd_lb - min(fad_ifit5_rmsd_lb)) / (max(fad_ifit5_rmsd_lb) - min(fad_ifit5_rmsd_lb)),
         NormalizedRMSD_ub_ifit5_fad = (fad_ifit5_rmsd_ub - min(fad_ifit5_rmsd_ub)) / (max(fad_ifit5_rmsd_ub) - min(fad_ifit5_rmsd_ub))) 
# Normalized values of rmsd and binding affinities NAD-IFIT5
data_nad_ifit5 = data_nad_ifit5 %>%
  mutate(NormalizedBindingAffinity_ifit5_nad = (nad_ifit5_energy - min(nad_ifit5_energy)) / (max(nad_ifit5_energy) - min(nad_ifit5_energy)),
         NormalizedRMSD_lb_ifit5_nad = (nad_ifit5_rmsd_lb - min(nad_ifit5_rmsd_lb)) / (max(nad_ifit5_rmsd_lb) - min(nad_ifit5_rmsd_lb)),
         NormalizedRMSD_ub_ifit5_nad = (nad_ifit5_rmsd_ub - min(nad_ifit5_rmsd_ub)) / (max(nad_ifit5_rmsd_ub) - min(nad_ifit5_rmsd_ub))) 
# Normalized values of rmsd and binding affinities PPP-IFIT5
#data_ppp_ifit5 = data_ppp_ifit5 %>%
 # mutate(NormalizedBindingAffinity_ifit5_ppp = (ppp_ifit5_energy - min(ppp_ifit5_energy)) / (max(ppp_ifit5_energy) - min(ppp_ifit5_energy)),
  #       NormalizedRMSD_lb_ifit5_ppp = (ppp_ifit5_rmsd_lb - min(ppp_ifit5_rmsd_lb)) / (max(ppp_ifit5_rmsd_lb) - min(ppp_ifit5_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit5_ppp = (ppp_ifit5_rmsd_ub - min(ppp_ifit5_rmsd_ub)) / (max(ppp_ifit5_rmsd_ub) - min(ppp_ifit5_rmsd_ub))) 
# Normalized values of rmsd and binding affinities PPP-IFIT5-MG
data_ppp_mg_ifit5 = data_ppp_mg_ifit5 %>%
  mutate(NormalizedBindingAffinity_ifit5_ppp_mg = (ppp_mg_energy - min(ppp_mg_energy)) / (max(ppp_mg_energy) - min(ppp_mg_energy)),
         NormalizedRMSD_lb_ifit5_ppp_mg = (ppp_mg_ifit5_rmsd_lb - min(ppp_mg_ifit5_rmsd_lb)) / (max(ppp_mg_ifit5_rmsd_lb) - min(ppp_mg_ifit5_rmsd_lb)),
         NormalizedRMSD_ub_ifit5_ppp_mg = (ppp_mg_ifit5_rmsd_ub - min(ppp_mg_ifit5_rmsd_ub)) / (max(ppp_mg_ifit5_rmsd_ub) - min(ppp_mg_ifit5_rmsd_ub))) 

# Normalized values of rmsd and binding affinities Ap4A-IFIT5
#data_Ap4A_ifit5 = data_Ap4A_ifit5 %>%
 # mutate(NormalizedBindingAffinity_ifit5_ap4a = (ap4a_ifit5_energy - min(ap4a_ifit5_energy)) / (max(ap4a_ifit5_energy) - min(ap4a_ifit5_energy)),
  #       NormalizedRMSD_lb_ifit5_ap4a = (ap4a_ifit5_rmsd_lb - min(ap4a_ifit5_rmsd_lb)) / (max(ap4a_ifit5_rmsd_lb) - min(ap4a_ifit5_rmsd_lb)),
   #      NormalizedRMSD_ub_ifit5_ap4a = (ap4a_ifit5_rmsd_ub - min(ap4a_ifit5_rmsd_ub)) / (max(ap4a_ifit5_rmsd_ub) - min(ap4a_ifit5_rmsd_ub))) 

# Calculating CompositeScore for each cap
# CompositeScore for cap0-IFIT5
#data_cap0_ifit5 = data_cap0_ifit5 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_cap0 +
  #         weight_lb * NormalizedRMSD_lb_ifit5_cap0 +
   #        weight_ub * NormalizedRMSD_ub_ifit5_cap0)
#head(data_cap0_ifit5)
# CompositeScore for cap1-IFIT5
#data_cap1_ifit5 = data_cap1_ifit5 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_cap1 +
  #         weight_lb * NormalizedRMSD_lb_ifit5_cap1 +
   #        weight_ub * NormalizedRMSD_ub_ifit5_cap1)
#head(data_cap1_ifit5)
# CompositeScore for cap2-IFIT5
#data_cap2_ifit5 = data_cap2_ifit5 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_cap2 +
  #         weight_lb * NormalizedRMSD_lb_ifit5_cap2 +
   #        weight_ub * NormalizedRMSD_ub_ifit5_cap2)
#head(data_cap2_ifit5)
# CompositeScore for FAD-IFIT5
data_fad_ifit5 = data_fad_ifit5 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_fad +
           weight_lb * NormalizedRMSD_lb_ifit5_fad +
           weight_ub * NormalizedRMSD_ub_ifit5_fad)
head(data_fad_ifit5)
# CompositeScore for NAD-IFIT5
data_nad_ifit5 = data_nad_ifit5 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_nad +
           weight_lb * NormalizedRMSD_lb_ifit5_nad +
           weight_ub * NormalizedRMSD_ub_ifit5_nad)
head(data_nad_ifit5)
# CompositeScore for ppp-IFIT5
#data_ppp_ifit5 = data_ppp_ifit5 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_ppp +
  #         weight_lb * NormalizedRMSD_lb_ifit5_ppp +
   #        weight_ub * NormalizedRMSD_ub_ifit5_ppp)
#head(data_ppp_ifit5)
# CompositeScore for ppp-IFIT5-MG
data_ppp_mg_ifit5 = data_ppp_mg_ifit5 %>%
  mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_ppp_mg +
           weight_lb * NormalizedRMSD_lb_ifit5_ppp_mg +
           weight_ub * NormalizedRMSD_ub_ifit5_ppp_mg)
head(data_ppp_mg_ifit5)
# CompositeScore for Ap4A-IFIT5
#data_Ap4A_ifit5 = data_Ap4A_ifit5 %>%
 # mutate(CompositeScore = weight_binding_affinity * NormalizedBindingAffinity_ifit5_ap4a +
  #         weight_lb * NormalizedRMSD_lb_ifit5_ap4a +
   #        weight_ub * NormalizedRMSD_ub_ifit5_ap4a)
#head(data_Ap4A_ifit5)

#Picking best dock
# A result with RMSD = 0 usually indicates that the pose is identical to the reference pose, meaning there was no movement or change in the ligand's position during docking.
# Best dock cap0-IFIT5
#ranked_docking_results_ifit5_cap0 = data_cap0_ifit5 %>%
 # arrange(CompositeScore)
#View(ranked_docking_results_ifit5_cap0)
#filtered_results_ifit5_cap0 = ranked_docking_results_ifit5_cap0 %>%
 # filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit5_cap0) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit5_cap0, file = "filtered_results_ifit5_cap0.csv")

# Best dock cap1-IFIT5
#ranked_docking_results_ifit5_cap1 = data_cap1_ifit5 %>%
 # arrange(CompositeScore)
#View(ranked_docking_results_ifit5_cap1)
#filtered_results_ifit5_cap1 = ranked_docking_results_ifit5_cap1 %>%
 # filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit5_cap1) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit5_cap1, file = "filtered_results_ifit5_cap1.csv")

# Best dock cap2-IFIT5
#ranked_docking_results_ifit5_cap2 = data_cap2_ifit5 %>%
 # arrange(CompositeScore)
#View(ranked_docking_results_ifit5_cap2)
#filtered_results_ifit5_cap2 = ranked_docking_results_ifit5_cap2 %>%
 # filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit5_cap2) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit5_cap2, file = "filtered_results_ifit5_cap2.csv")

# Best dock FAD-IFIT5
ranked_docking_results_ifit5_fad = data_fad_ifit5 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit5_fad)
filtered_results_ifit5_fad = ranked_docking_results_ifit5_fad %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit5_fad) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit5_fad, file = "filtered_results_ifit5_fad.csv")

# Best dock NAD-IFIT5
ranked_docking_results_ifit5_nad = data_nad_ifit5 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit5_nad)
filtered_results_ifit5_nad = ranked_docking_results_ifit5_nad %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit5_nad) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit5_nad, file = "filtered_results_ifit5_nad.csv")

# Best dock ppp-IFIT5
ranked_docking_results_ifit5_ppp = data_ppp_ifit5 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit5_ppp)
filtered_results_ifit5_ppp = ranked_docking_results_ifit5_ppp %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit5_ppp) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit5_ppp, file = "filtered_results_ifit5_ppp.csv")

# Best dock ppp-IFIT5-MG
ranked_docking_results_ifit5_ppp_mg = data_ppp_mg_ifit5 %>%
  arrange(CompositeScore)
View(ranked_docking_results_ifit5_ppp_mg)
filtered_results_ifit5_ppp_mg = ranked_docking_results_ifit5_ppp_mg %>%
  filter(V3 != 0 & V4 != 0)
print(filtered_results_ifit5_ppp_mg) #Choosing the best docking, exlcuding RMSD's = 0
write.csv(filtered_results_ifit5_ppp_mg, file = "filtered_results_ifit5_ppp_mg.csv")

# Best dock Ap4A-IFIT5
#ranked_docking_results_ifit5_ap4a = data_Ap4A_ifit5 %>%
 # arrange(CompositeScore)
#View(ranked_docking_results_ifit5_ap4a)
#filtered_results_ifit5_ap4a = ranked_docking_results_ifit5_ap4a %>%
 # filter(V3 != 0 & V4 != 0)
#print(filtered_results_ifit5_ap4a) #Choosing the best docking, exlcuding RMSD's = 0
#write.csv(filtered_results_ifit5_ap4a, file = "filtered_results_ifit5_ap4a.csv")

#Creating a Heatmap for all best docks binding energies (according to R CHARTS by R CODER 2024: https://r-charts.com/correlation/heat-map-ggplot2/)
#First, let's place all the values in a dataframe
bd_IFIT5_energies = c(-9.9, -9.2, -10.4) #these values correspond to the selected docking for each case
bd_IFIT5_caps = c("FAD", "PPP", "NAD")
bd_IFIT1_energies = c(-13.0, -10.5, -11.5, -10.5, -10.4) #these values correspond to the selected docking for each cap
bd_IFIT1_caps = c("FAD", "NAD", "Cap1", "Cap2", "Cap0")
all_energies_caps_IFIT5 = data.frame(Energy = bd_IFIT5_energies, Cap = bd_IFIT5_caps)
all_energies_caps_IFIT1 = data.frame(Energy = bd_IFIT1_energies, Cap = bd_IFIT1_caps)
combined_data_energies = rbind(all_energies_caps_IFIT5, all_energies_caps_IFIT1)
combined_data_energies$Protein = rep(c("IFIT5", "IFIT1"), times = c(nrow(all_energies_caps_IFIT5), nrow(all_energies_caps_IFIT1)))

#Now, the heatmap that fits both best docks binding energies from IFIT5 and IFIT1
ggplot(combined_data_energies, aes(x = Protein, y = Cap, fill = Energy)) +
  geom_tile(color = "black") + ggtitle("Best dockings's binding energies") +
  scale_fill_gradientn(colors = hcl.colors(20, "Cork")) +
  coord_fixed() + geom_text(aes(label = Energy), color = "white", size = 4.0) +
  guides(fill = guide_colourbar(title = "Energy (kcal/mol)"))

#Let's do the same for the RMSD's for the best dockings! Because this could give an insight of how much the ligand had to change it's conformation from the original pose
bd_IFIT5_rmsd_lb = c(4.409, 2.352, 3.195)
bd_IFIT5_rmsd_lb_caps = c("PPP", "NAD", "FAD")
bd_IFIT5_rmsd_ub = c(8.850, 4.362, 10.833)
bd_IFIT1_rmsd_lb = c(1.920, 2.480, 1.700, 1.302, 1.783)
bd_IFIT1_rmsd_lb_caps = c("Cap0", "Cap1", "Cap2", "FAD", "NAD")
bd_IFIT1_rmsd_ub = c(1.320, 1.526, 2.469, 2.129, 2.851)
bd_IFIT1_rmsd_ub_caps = c("Cap0", "Cap1", "Cap2", "FAD", "NAD")

all_rmsd_lb_IFIT1_caps = data.frame(RMSD_l.b = bd_IFIT1_rmsd_lb, Caps= bd_IFIT1_rmsd_lb_caps)
all_rmsd_ub_IFIT1_caps = data.frame(RMSD_u.b = bd_IFIT1_rmsd_ub, Caps= bd_IFIT1_rmsd_lb_caps)
all_rmsd_lb_IFIT5_caps = data.frame(RMSD_l.b = bd_IFIT5_rmsd_lb, Caps = bd_IFIT5_rmsd_lb_caps)
all_rmsd_ub_IFIT5_caps = data.frame(RMSD_u.b = bd_IFIT5_rmsd_ub, Caps = bd_IFIT5_rmsd_lb_caps)

#Now, the heatmap for RMSD's lb and ub for IFIT1's best dockings
combined_rmsd_lb = rbind(all_rmsd_lb_IFIT1_caps, all_rmsd_lb_IFIT5_caps)
combined_rmsd_lb$Protein = rep(c("IFIT1", "IFIT5"), times = c(nrow(all_rmsd_lb_IFIT1_caps), nrow(all_rmsd_lb_IFIT5_caps)))

ggplot(combined_rmsd_lb, aes(x = Protein, y = Caps, fill = RMSD_l.b)) +
  geom_tile(color = "black") + ggtitle("Best dockings's RMSD's_l.b.") +
  scale_fill_gradientn(colors = hcl.colors(20, "Heat")) +
  coord_fixed() + geom_text(aes(label = RMSD_l.b), color = "black", size = 4.0) +
  guides(fill = guide_colourbar(title = "RMSD_l.b"))

combined_rmsd_ub = rbind(all_rmsd_ub_IFIT1_caps, all_rmsd_ub_IFIT5_caps)
combined_rmsd_ub$Protein = rep(c("IFIT1", "IFIT5"), times = c(nrow(all_rmsd_ub_IFIT1_caps), nrow(all_rmsd_ub_IFIT5_caps)))

ggplot(combined_rmsd_ub, aes(x = Protein, y = Caps, fill = RMSD_u.b)) +
  geom_tile(color = "black") + ggtitle("Best dockings's RMSD's_u.b.") +
  scale_fill_gradientn(colors = hcl.colors(20, "Heat")) +
  coord_fixed() + geom_text(aes(label = RMSD_u.b), color = "black", size = 4.0) +
  guides(fill = guide_colourbar(title = "RMSD_u.b"))

# Referece of language model used to correct, refine and optimize this code:
# (OpenAI. (2025). ChatGPT (Jan 20 version) [Large language model]. Retrieved from https://chat.openai.com)

# Wilcoxon test or Mann-Whitney test????
# Used for comparing the means of two independent samples (x & y), in this case, binding affinity means of each cap from both IFIT's. it's comparing two independent groups of samples. 
# It’s used when your data are not normally distributed.
# alternative: the alternative hypothesis. Allowed value is one of “two.sided” (default), “greater” or “less”.
# wilcox.test(x, y, alternative = "two.sided")
