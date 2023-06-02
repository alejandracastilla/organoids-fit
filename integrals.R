#-------------- Intro -------------------
# username: alecastillab
# Name: María Alejandra Castilla Bolaños
# Affiliation: Ph.D. Student in Medical Biophysics, Near's lab, University of Toronto
# Physical Sciences, Sunnybrook
#-------------- Declaring integral vector -------------------
integrals_scan4 <- c(1013681.83077109,	1740065.72646358,	3158945.37129312,	492051.309594797,	555472.774266667,	358637.817561697,	937387.584418775,	1259140.66491945,	995873.107287070)

integrals_scan3 <- c(6528424.54934505,	5314296.16059844,	24837035.5865943,	1099090.19830726,	2479930.32369288,	1209053.96561301,	2333201.41742437,	1036917.02000616,	2345479.47206918)
#-------------- Saving integrals into a data frame------------
df_integrals_scan3 <- as.data.frame(integrals_scan3)
df_integrals_scan3$m_org <- c(0.0113, 0.0158, 0.0091, 0.0039, 0.0098, 0.0043, 0.0076, 0.0089, 0.0038)
# clean data, by deleting third row
cleanoutput <- df_integrals_scan3[-c(3),]

plot(df_integrals_scan3$m_org, df_integrals_scan3$integrals_scan3)

plot(cleanoutput$m_org, cleanoutput$integrals_scan3)
#-------------- Fitting a model m_org~integrals ------------
lm_mass = lm(m_org~integrals_scan3, data = df_integrals_scan3) #Create the linear regression
lm_mass = lm(m_org~integrals_scan3, data = cleanoutput) #Create the linear regression for cleanoutput
summary(lm_mass) #Review the results

#----------------plotting ----------------------------------
library(tidyverse)
library(ggplot2)
#install.packages("ggpubr") #to put the equation in the plot
library(ggpubr)


regression <- ggplot(cleanoutput, aes(x = m_org, y = integrals_scan3)) + 
  geom_point()
# adding the line to the plot
regression <- regression + geom_smooth(method="lm", col="black")

# adding the equation to the plot
regression <- regression + stat_regline_equation(label.x = 0.01, label.y = 8e6)
regression


