################################################################################

                              #Plots for MMPBSA analysis

################################################################################
# MIT License
# Copyright (c) 2025 Rodolfo Adrian Arriaga Rivera
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction.
# See the LICENSE file for full license text.
################################################################################

#Setting working directory
setwd("~/Documents/RNA_Genomics")

#calling ggplot2 for plotting and some other tools for df processing
library(dplyr)
library(tidyr)
library(ggplot2)

#First, we need to input the data into R, from an Excel spreasheet (I know, don't judge)
#in which I placed the DeltaTotal Energy Contributions of the MMPBSA analysis of each cap
gmx_ifit1_data= readxl::read_xlsx("data_gmxMMPBSA_ifit1.xlsx")

#Grouping data into objects called ifit1_cap1_data and ifit1_cap2_data (and so on, for each cap)
ifit1_cap0_data=select(gmx_ifit1_data, `Frame #`, `DeltaInicialC0`, `DeltaFinalC0`)
ifit1_cap0_data
ifit1_cap1_data=select(gmx_ifit1_data, `Frame #`, `DeltaInicialC1`, `DeltaFinalC1`)
ifit1_cap1_data
ifit1_cap2_data=select(gmx_ifit1_data, `Frame #`, `DeltaInicialC2`, `DeltaFinalC2`)
ifit1_cap2_data
ifit1_fad_data=select(gmx_ifit1_data, `Frame #`, `DeltaInicialFAD`, `DeltaFinalFAD`)
ifit1_fad_data
ifit1_nad_data=select(gmx_ifit1_data, `Frame #`, `DeltaInicialNAD`, `DeltaFinalNAD`)
ifit1_nad_data
ifit1_fad_ex_data=select(gmx_ifit1_data, `Frame #`, `DeltaInicialFADex`, `DeltaFinalFADex`)
ifit1_fad_ex_data

# Convert to long format for plotting Cap0 using pivot_longer (this merges DeltaInicial and DeltaFinal columns into one), as suggested by ChatGPT, referencing: https://tidyr.tidyverse.org/reference/pivot_longer.html
ifit1_cap0_long = ifit1_cap0_data %>%
  pivot_longer(cols = c(`DeltaInicialC0`, `DeltaFinalC0`), 
               names_to = "Initial/Final", #this makes both columns have the name Initial/Final
               values_to = "Delta") #And this takes the actual data into another column called Delta
head(ifit1_cap0_long) #check the data

# Calculate mean ΔG for each group for plotting for Cap0
mean_values_cap0 = ifit1_cap0_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) #calculates mean of Delta column and removes any NA spaces
mean_values_cap0

# Plot Cap0 energies (suggested by ChatGPT), also referencing:  ggplot2 line plot : Quick start guide - R software and data visualization (https://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization)
ggplot(ifit1_cap0_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_cap0, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT1-Cap0", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaInicialC0" = "Initial Binding", 
                                  "DeltaFinalC0" = "Final Binding")) + 
  theme_light()

# Convert to long format for plotting Cap1 using pivot_longer (this merges DeltaInicial and DeltaFinal columns into one), as suggested by ChatGPT, referencing: https://tidyr.tidyverse.org/reference/pivot_longer.html
ifit1_cap1_long = ifit1_cap1_data %>%
  pivot_longer(cols = c(`DeltaInicialC1`, `DeltaFinalC1`), 
               names_to = "Initial/Final", #this makes both columns have the name Initial/Final
               values_to = "Delta") #And this takes the actual data into another column called Delta
head(ifit1_cap1_long) #check the data

# Calculate mean ΔG for each group for plotting for Cap1
mean_values_cap1 = ifit1_cap1_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) #calculates mean of Delta column and removes any NA spaces
mean_values_cap1

# Plot Cap1 energies (suggested by ChatGPT), also referencing:  ggplot2 line plot : Quick start guide - R software and data visualization (https://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization)
ggplot(ifit1_cap1_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_cap1, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT1-Cap1", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaInicialC1" = "Initial Binding", 
                                  "DeltaFinalC1" = "Final Binding")) + 
  theme_light()

#Now, same deal, but for Cap2 (and so on, for each cap)
# Convert to long format for plotting Cap2 using pivot_longer (this merges DeltaInicial and DeltaFinal columns into one), as suggested by ChatGPT, referencing: https://tidyr.tidyverse.org/reference/pivot_longer.html
ifit1_cap2_long = ifit1_cap2_data %>%
  pivot_longer(cols = c(`DeltaInicialC2`, `DeltaFinalC2`), 
               names_to = "Initial/Final", #this makes both columns have the name Initial/Final
               values_to = "Delta") #And this makes the actual data into another column called Delta
head(ifit1_cap2_long) #check the data

# Calculate mean ΔG for each group for plotting for Cap2
mean_values_cap2 = ifit1_cap2_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = TRUE)) #calculates mean of Delta column and removes any NA spaces
mean_values_cap2

# Plot Cap2 energies (suggested by ChatGPT) also referencing:  ggplot2 line plot : Quick start guide - R software and data visualization 
# (https://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization), 
ggplot(ifit1_cap2_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_cap2, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Dashed mean lines
  labs(title = "Binding Free Energy Over Time - IFIT1-Cap2", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaInicialC2" = "Initial Binding", 
                                  "DeltaFinalC2" = "Final Binding")) + 
  theme_light()

# Convert to long format for plotting FAD using pivot_longer
ifit1_fad_long = ifit1_fad_data %>%
  pivot_longer(cols = c(`DeltaInicialFAD`, `DeltaFinalFAD`), 
               names_to = "Initial/Final", #this makes both columns have the name Initial/Final
               values_to = "Delta") #And this takes the actual data into another column called Delta
head(ifit1_fad_long) #check the data

# Calculate mean ΔG for each group for plotting for FAD
mean_values_fad = ifit1_fad_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) #calculates mean of Delta column and removes any NA spaces
mean_values_fad

# Plot FAD energies (suggested by ChatGPT), also referencing:  ggplot2 line plot : Quick start guide - R software and data visualization 
# (https://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization)
ggplot(ifit1_fad_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_fad, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT1-FAD", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaInicialFAD" = "Initial Binding",
                                  "DeltaFinalFAD" = "Final Binding")) + 
  theme_light()

# Convert to long format for plotting NAD using pivot_longer
ifit1_nad_long = ifit1_nad_data %>%
  pivot_longer(cols = c(`DeltaInicialNAD`,`DeltaFinalNAD`), 
               names_to = "Initial/Final", #this makes both columns have the name Initial/Final
               values_to = "Delta") #And this takes the actual data into another column called Delta
head(ifit1_nad_long) #check the data

# Calculate mean ΔG for each group for plotting for NAD
mean_values_nad = ifit1_nad_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) #calculates mean of Delta column and removes any NA spaces
mean_values_nad

# Plot NAD energies (suggested by ChatGPT), also referencing:  ggplot2 line plot : Quick start guide - R software and data visualization 
# (https://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization)
ggplot(ifit1_nad_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_nad, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT1-NAD", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaFinalNAD" = "Final Binding",
                                  "DeltaInicialNAD" = "Initial Binding")) + 
  theme_light()

##################################################################################################################

                                           #IFIT5 MMPBSA analysis

##################################################################################################################

#Now, same deal but for IFIT5-cap data

#Importing data from Excel spreadsheet
gmx_ifit5_data= readxl::read_xlsx("data_gmxMMPBSA_ifit5.xlsx")

# Extracting data for each cap
ifit5_NAD_data=select(gmx_ifit5_data, `Frame #`, `DeltaInicialNAD`, `DeltaFinalNAD`)
ifit5_NAD_data
ifit5_FAD_data=select(gmx_ifit5_data, `Frame #`, `DeltaInicialFAD`, `DeltaFinalFAD`)
ifit5_FAD_data
ifit5_ppp_mg_data=select(gmx_ifit5_data, `Frame #`, `DeltaInicialPPP`, `DeltaFinalPPP`)

# Convert to long format for plotting NAD using pivot_longer 
ifit5_NAD_long = ifit5_NAD_data %>%
  pivot_longer(cols = c(`DeltaInicialNAD`, `DeltaFinalNAD`), 
               names_to = "Initial/Final",
               values_to = "Delta")
head(ifit5_NAD_long)

# Calculate mean ΔG for each group for plotting for NAD
mean_values_ifit5_NAD = ifit5_NAD_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) 
mean_values_ifit5_NAD

# Plot NAD energies 
ggplot(ifit5_NAD_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_ifit5_NAD, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT5-NAD", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaInicialNAD" = "Initial Binding", 
                                  "DeltaFinalNAD" = "Final Binding")) + 
  theme_light()

# Convert to long format for plotting FAD using pivot_longer 
ifit5_FAD_long = ifit5_FAD_data %>%
  pivot_longer(cols = c(`DeltaInicialFAD`, `DeltaFinalFAD`), 
               names_to = "Initial/Final",
               values_to = "Delta")
head(ifit5_FAD_long)

# Calculate mean ΔG for each group for plotting for FAD
mean_values_ifit5_FAD = ifit5_FAD_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) 
mean_values_ifit5_FAD

# Plot FAD energies 
ggplot(ifit5_FAD_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_ifit5_FAD, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT5-FAD", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(labels = c("DeltaFinalFAD" = "Final Binding",
                                  "DeltaInicialFAD" = "Initial Binding")) + 
  theme_light()

# Convert to long format for plotting PPP using pivot_longer 
ifit5_PPP_long = ifit5_ppp_mg_data %>%
  pivot_longer(cols = c(`DeltaInicialPPP`, `DeltaFinalPPP`), 
               names_to = "Initial/Final",
               values_to = "Delta")
head(ifit5_PPP_long)

# Calculate mean ΔG for each group for plotting for PPP
mean_values_ifit5_PPP = ifit5_PPP_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) 
mean_values_ifit5_PPP

# Plot PPP energies 
ggplot(ifit5_PPP_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_ifit5_PPP, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT5-PPP", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaInicialPPP" = "Initial Binding", 
                                  "DeltaFinalPPP" = "Final Binding")) + 
  theme_light()


# Convert to long format for plotting FAD_ex using pivot_longer (this merges DeltaInicial and DeltaFinal columns into one), as suggested by ChatGPT, referencing: https://tidyr.tidyverse.org/reference/pivot_longer.html
ifit1_fadex_long = ifit1_fad_ex_data %>%
  pivot_longer(cols = c(`DeltaInicialFADex`, `DeltaFinalFADex`), 
               names_to = "Initial/Final", #this makes both columns have the name Initial/Final
               values_to = "Delta") #And this takes the actual data into another column called Delta
head(ifit1_fadex_long) #check the data

# Calculate mean ΔG for each group for plotting for FAD_ex
mean_values_fadex = ifit1_fadex_long %>%
  group_by(`Initial/Final`) %>%
  summarize(mean_Delta = mean(Delta, na.rm = T)) #calculates mean of Delta column and removes any NA spaces
mean_values_fadex

# Plot FAD_ex energies (suggested by ChatGPT), also referencing:  ggplot2 line plot : Quick start guide - R software and data visualization (https://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization)
ggplot(ifit1_fadex_long, aes(x = `Frame #`, y = Delta, color = `Initial/Final`)) +
  geom_line() + 
  geom_hline(data = mean_values_fadex, aes(yintercept = mean_Delta, color = `Initial/Final`), 
             linetype = "dashed", size = 0.3, color="black") +  # Adds mean lines
  labs(title = "Binding Free Energy Over Time - IFIT1-FAD_extended", 
       x = "Frame #", 
       y = "ΔG (kcal/mol)") +
  scale_color_discrete(name = "Initial/Final", 
                       labels = c("DeltaInicialFADex" = "Initial Binding", 
                                  "DeltaFinalFADex" = "Final Binding")) + 
  theme_light()
