#load spant library
library(spant)


### Run from making metabolites first
## until fit rest not found
#now run the firstr part until line 68

################################################################################
#reference sample fit
#volume = 30ul
#concentration = 20mM
ref_list <- list(ala2)
ref_basis <- sim_basis(ref_list, pul_seq = seq_pulse_acquire, acq_paras = def_acq_paras(N = 8124, fs = 10000, ft = 5.0028e8, ref = 4.65, nuc="1H"))
refname <- "~/Documents/data/studies/organoid_NMR/hESC_LT/alanineRef/ala29_lcm"
ref_data <- read_mrs(refname, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref = 4.65)
opts <- abfit_opts(ppm_left=4.8)
?opts
ref_res <- fit_mrs(ref_data, ref_basis, opts=opts)
plot(ref_res)

################################################################################
#fitting data

#load metabolites and simulate basis set
met_list <- list(ace2, lac2, ala2, get_mol_paras("cho"), get_mol_paras("glc"), gln2, gly2, get_mol_paras("gpc"), get_mol_paras("ins"), get_mol_paras("naa"), get_mol_paras("pch"), get_mol_paras("gaba"), cr2, etoh, htau, get_mol_paras("glu"))
basis <- sim_basis(met_list, pul_seq = seq_pulse_acquire, acq_paras = def_acq_paras(N = 8124, fs = 10000, ft = 5.0028e8, ref = 4.65, nuc="1H"))
stackplot(basis, xlim = c(4.8, 0.2), y_offset = 50, labels = basis$names)

#load data files/spectra
b <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/B/liqph_lcm"
c <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/C/liqph_lcm"
d <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/D/liqph_lcm"
e <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/E/liqph_lcm"
f <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/F/liqph_lcm"
g <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/G/liqph_lcm"
h <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/H/liqph_lcm"
i <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/I/liqph_lcm"
k <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/K/liqph_lcm"
l <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/L/liqph_lcm"
s <- "~/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/S/liqph_lcm"
t <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/T/liqph_lcm"
u <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/U/liqph_lcm"
aa <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/AA/liqph_lcm"
bb <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/BB/liqph_lcm"
cc <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/CC/liqph_lcm"
dd <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/DD/liqph_lcm"
ee <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/EE/liqph_lcm"
ff <- "/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/FF/liqph_lcm"


#read files
files <- c(b, c, d, e, f, g, h, i, k, l, s, t, u, aa, bb, cc, dd, ee, ff)
data <- read_mrs(files, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref=4.65)

#optional: check data by plotting

for (spectra in data) {
  plot(spectra)}

#set/customize options
opts <- abfit_opts(ppm_left=4.8)
#to see options:
?abfit_opts

#fit analyis

fit_res <- fit_mrs(data, basis, opts=opts)

#plot fits

for (fit in fit_res) {
  plot(fit)}

################################################################################
# manually fit/anaylze dataset

# #B
# met_list <- list(get_mol_paras("lac"), ala2, get_mol_paras("cho"), get_mol_paras("glc"), gln2, gly2, get_mol_paras("gpc"), get_mol_paras("ins"), get_mol_paras("naa"), get_mol_paras("pch"), get_mol_paras("gaba"), cr2, custom_mol, htau)
# basis <- sim_basis(met_list, pul_seq = seq_pulse_acquire, acq_paras = def_acq_paras(N = 8124, fs = 10000, ft = 5.0028e8, ref = 4.65, nuc="1H"))
# stackplot(basis, xlim = c(4.8, 0.2), y_offset = 50, labels = basis$names)
# b <- "~/Documents/spant/B"
# b_data <- read_mrs(b, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref=4.65)
# opts <- abfit_opts(ppm_left=4.8)
# fit_res_b <- fit_mrs(b_data, basis, opts=opts)
# plot(fit_res_b)

################################################################################
#making metabolites

#ethanol

n_etoh <- rep("1H", 5)
chem_shift_etoh <- c(1.16455, 1.16455, 1.16455, 3.67, 3.67)

#make matrix with J coupling values

j_coupling_mat_etoh <- matrix(0, 5, 5)
j_coupling_mat_etoh[1, 4] <- 7.1
j_coupling_mat_etoh[1,5] <- 7.1
j_coupling_mat_etoh[2,4] <- 7.1
j_coupling_mat_etoh[2,5] <- 7.1
j_coupling_mat_etoh[3,4] <- 7.1
j_coupling_mat_etoh[3,5] <- 7.1

spin_group_etoh <- list(nucleus = n_etoh, chem_shift = chem_shift_etoh, 
                        j_coupling_mat = j_coupling_mat_etoh, scale_factor = 1,
                        lw = 2, lg = 0)

source <- "test"

etoh <- list(spin_groups = list(spin_group_etoh), name = "EtOH",
                   source = source, full_name = "Ethanol")

etoh |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(etoh) <- "mol_parameters"

#hypotaurine
n_htau <- rep("1H", 4)
chem_shift_htau <- c(3.32, 3.32, 2.626, 2.626)
j_coupling_mat_htau <- matrix(0, 4, 4)
j_coupling_mat_htau[1,3] <- 6.742
j_coupling_mat_htau[1,4] <- 6.464
j_coupling_mat_htau[2,3] <- 6.403
j_coupling_mat_htau[2,4] <- 6.792

spin_group_htau <- list(nucleus = n_htau, chem_shift = chem_shift_htau, 
                        j_coupling_mat = j_coupling_mat_htau, scale_factor = 1,
                        lw = 2, lg = 0)

source <- "test"
htau <- list(spin_groups = list(spin_group_htau), name = "hTau",
             source = source, full_name = "Hypotaurine")


htau |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(htau) <- "mol_parameters"

#alanine
n_ala2 <- rep("1H", 4)
chem_shift_ala2 <- c(3.74, 1.455, 1.455, 1.455)
j_coupling_mat_ala2 <- matrix(0, 4, 4)
j_coupling_mat_ala2[2, 1] <- 7.234
j_coupling_mat_ala2[3, 1] <- 7.234
j_coupling_mat_ala2[4, 1] <- 7.234
j_coupling_mat_ala2[3, 2] <- -14.366
j_coupling_mat_ala2[4, 2] <- -14.366
j_coupling_mat_ala2[4, 3] <- -14.366

spin_group_ala2 <- list(nucleus = n_ala2, chem_shift = chem_shift_ala2, 
                        j_coupling_mat = j_coupling_mat_ala2, scale_factor = 1,
                        lw = 2, lg = 0)
source <- "test"
#called it Alanine2 to differentiate original versus one with edited chemical shifts
ala2 <- list(spin_groups = list(spin_group_ala2), name = "Ala",
             source = source, full_name = "Alanine2")
ala2 |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(ala2) <- "mol_parameters"

#creatine

n_cr2 <- rep("1H", 6)
chem_shift_cr2 <- c(6.649, 3.913, 3.913, 3.01, 3.01, 3.01)
j_coupling_mat_cr2 <- matrix(0, 6, 6)

spin_group_cr2 <- list(nucleus = n_cr2, chem_shift = chem_shift_cr2, 
                       j_coupling_mat = j_coupling_mat_cr2, scale_factor = 1,
                       lw = 2, lg = 0)

source <- "test"
cr2 <- list(spin_groups = list(spin_group_cr2), name = "Cr",
            source = source, full_name = "Creatine2")
cr2 |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(cr2) <- "mol_parameters"

#glycine

n_gly2 <- rep("1H", 2)
chem_shift_gly2 <- c(3.529, 3.529)
j_coupling_mat_gly2 <- matrix(0, 2, 2)

spin_group_gly2 <- list(nucleus = n_gly2, chem_shift = chem_shift_gly2, 
                        j_coupling_mat = j_coupling_mat_gly2, scale_factor = 1,
                        lw = 2, lg = 0)

source <- "test"
gly2 <- list(spin_groups = list(spin_group_gly2), name = "Gly",
             source = source, full_name = "Glycine2")
gly2 |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(gly2) <- "mol_parameters"

#glutamine

n_gln2 <- rep("1H", 5)
chem_shift_gln2 <- c(3.753, 2.129, 2.109, 2.417, 2.417)
j_coupling_mat_gln2 <- matrix(0, 5, 5)
j_coupling_mat_gln2[1, 2] <- 5.847 
j_coupling_mat_gln2[1, 3] <- 6.5
j_coupling_mat_gln2[2, 3] <- -14.504
j_coupling_mat_gln2[2, 4] <- 9.165
j_coupling_mat_gln2[2, 5] <- 6.347
j_coupling_mat_gln2[3, 4] <- 6.324
j_coupling_mat_gln2[3, 5] <- 9.209
j_coupling_mat_gln2[4, 5] <- -15.371

spin_group_gln2 <- list(nucleus = n_gln2, chem_shift = chem_shift_gln2, 
                        j_coupling_mat = j_coupling_mat_gln2, scale_factor = 1,
                        lw = 2, lg = 0)

source <- "test"
gln2 <- list(spin_groups = list(spin_group_gln2), name = "Gln",
             source = source, full_name = "Glutamine2")
gln2 |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(gln2) <- "mol_parameters"

#acetate
n_ace2 <- rep("1H")
chem_shift_ace2 <- c(1.9004)
j_coupling_mat_ace2 <- matrix(0, 1, 1)

spin_group_ace2 <- list(nucleus = n_ace2, chem_shift = chem_shift_ace2, 
                        j_coupling_mat = j_coupling_mat_ace2, scale_factor = 3,
                        lw = 2, lg = 0)

source <- "test"
ace2 <- list(spin_groups = list(spin_group_ace2), name = "Ace",
             source = source, full_name = "Acetate")


ace2 |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(ace2) <- "mol_parameters"

#glucose (test: 2 spin groups)
# n_a_glc2 <- rep("1H", 7)
# chem_shift_a_glc2 <- c(5.216, 3.519, 3.698, 3.395, 3.822, 3.826, 3.749)
# j_coupling_mat_a_glc2 <- matrix(0, 7, 7)
# j_coupling_mat_a_glc2[1, 2] <- 3.8
# j_coupling_mat_a_glc2[2, 3] <- 9.6
# j_coupling_mat_a_glc2[3, 4] <- 9.4
# j_coupling_mat_a_glc2[4, 5] <- 9.9
# j_coupling_mat_a_glc2[5, 6] <- 1.5
# j_coupling_mat_a_glc2[5, 7] <- 6
# j_coupling_mat_a_glc2[6, 7] <- -12.1
# 
# spin_group_a_glc2 <- list(nucleus = n_a_glc2, chem_shift = chem_shift_a_glc2, 
#                         j_coupling_mat = j_coupling_mat_a_glc2, scale_factor = 0.36,
#                         lw = 2, lg = 0)
# 
# n_b_glc2 <- rep("1H", 7)
# chem_shift_b_glc2 <- c(4.63, 3.23, 3.473, 3.387, 3.45, 3.882, 3.707)
# 
# j_coupling_mat_b_glc2 <- matrix(0, 7, 7)
# j_coupling_mat_b_glc2[1, 2] <- 8
# j_coupling_mat_b_glc2[2, 3] <- 9.1
# j_coupling_mat_b_glc2[3, 4] <- 9.4
# j_coupling_mat_b_glc2[4, 5] <- 8.9
# j_coupling_mat_b_glc2[5, 6] <- 1.6
# j_coupling_mat_b_glc2[5, 7] <- 5.4
# j_coupling_mat_b_glc2[6, 7] <- -12.3
# 
# spin_group_b_glc2 <- list(nucleus = n_b_glc2, chem_shift = chem_shift_b_glc2, 
#                           j_coupling_mat = j_coupling_mat_b_glc2, scale_factor = 0.64,
#                           lw = 2, lg = 0)
# source <- "test"
# glc2 <- list(spin_groups = list(spin_group_a_glc2, spin_group_b_glc2), name = "glc2", source=source, full_name="glucose2")
# glc2 |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
# class(glc2) <- "mol_parameters"

#lactate

n_lac2 <- rep("1H", 4)
chem_shift_lac2 <- c(4.0974, 1.309, 1.309, 1.309)
j_coupling_mat_lac2 <- matrix(0, 4, 4)
j_coupling_mat_lac2[1,2] <- 6.933
j_coupling_mat_lac2[1,3] <- 6.933
j_coupling_mat_lac2[1, 4] <- 6.933

spin_group_lac2 <- list(nucleus = n_lac2, chem_shift = chem_shift_lac2, 
                        j_coupling_mat = j_coupling_mat_lac2, scale_factor = 1,
                        lw = 2, lg = 0)

source <- "test"
lac2 <- list(spin_groups = list(spin_group_lac2), name = "Lac",
             source = source, full_name = "Lactate2")


lac2 |> sim_mol() |> lb(2) |> zf() |> plot(xlim = c(4.4, 0.5))
class(lac2) <- "mol_parameters"

################################################################################
#organoid areas in mm^3

d = c(2.1, 2, 1.5, 1.5, 1.5, 1.5)
vol_mm = (4/3)*pi*(d/2)^3
print(vol_mm)

#conversion to litres

vol_l = vol_mm*(0.001)*(1/1000)
print(vol_l)

################################################################################
#quantify concentrations using for loop 

output <- matrix(nrow = length(fit_res), ncol = 16)
  
index=1
for (i in 1:length(fit_res)) {
  metabs <- fit_res[[i]]$res_tab[, 6:21]
  for (j in 1:length(metabs)) {
    conc = (0.02)*(metabs[[j]]/120629.2)*(0.00003/vol_l[[index]])
    output[i, j] <- conc
  }
  index=(index+1)
}

# clean data, by deleting third row
cleanoutput <- output[-c(3),]

#write data frame to csv file
write.csv(cleanoutput,"~/Documents/test.csv", row.names = FALSE)

#write vector as data frame
df <- as.data.frame(cleanoutput)

#add age column
df$Age <- c(103, 103, 121, 121)

#rewrite/rename column names
colnames(df) <- c("Ace", "Lac", "Ala", "Cho", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")

################################################################################
#plot the fit for a single metabolite?

plot(fit_res[[1]]$fits$`1.fit`$PPMScale, fit_res[[1]]$fits$`1.fit`$Gln, type="l")

################################################################################
#mean standard deviation for all metabolites

for (fit in fit_res) {
  meansd = sum(fit$res_tab[,22:37])/length(fit$res_tab[,22:37])
  print(meansd)
}

################################################################################
#making graphs using ggplot
library(ggplot2)

#creating a list where the plots will be stored
plot_list = list()

for (i in 1:16) {
  df <- as.data.frame(cleanoutput)
  df$Age <- c(103, 103, 121, 121)
  colnames(df) <- c("Ace", "Lac", "Ala", "Cho", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")
  xvar <- df$Age
  yvar <- df[, i]
  tempdata <- data.frame(xvar, yvar)
  plot_list[[i]] <- ggplot(tempdata, aes(x=xvar, y=yvar)) + 
            geom_point() +
            geom_smooth(method="lm", formula=y~x, se=FALSE, color = "#A3A7E2") +
            labs(x = "Age (days)", y = "Concentration (M)") +
            ylim(0, max(df[,i])*1.1) +
            ggtitle(names(df)[i]) +
            theme(panel.background = element_rect(fill = "white"),
                  panel.grid = element_line(color = "#EEEEEE"))
}

#note: i created a for loop that makes a new temporary df each time since only the
#last metabolite was being plotted every time without doing so

#plot graphs side by side

library(gridExtra)

out <- do.call(grid.arrange, plot_list)
plot(out)

#do.call: calls a function and a list and passes the elements of the list 
#as arguments for the function
#in matlab, i think the "equivalent" would be listing the elements in curly brackets
#but i couldn't do that in r

################################################################################
#statistical analysis

#dataset with  metabolites with 0 concentration removed (choline)

stats_output <- cleanoutput[,-c(4)]

for (i in 1:15) {
  statsdf <- as.data.frame(stats_output)
  statsdf$Age <- c(103, 103, 121, 121)
  colnames(statsdf) <- c("Ace", "Lac", "Ala", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")
  Age <- statsdf$Age
  Concentration <- statsdf[, i]
  tempdf <- data.frame(Age, Concentration)
  ttest = t.test(Concentration ~ Age, data = tempdf)
  print(ttest)
}

################################################################################

