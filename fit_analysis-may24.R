#-------------- Intro -------------------
# username: alecastillab
# Name: María Alejandra Castilla Bolaños
# Affiliation: Ph.D. Student in Medical Biophysics, Near's lab, University of Toronto
# Physical Sciences, Sunnybrook
# ------------- installing packages ------------------------------------------------------------------
# install.packages("devtools")
# install.packages("plotly")
# install.packages("ggplotly") #not available for this verion of R?
# devtools::install_git("https://github.com/cbarbu/R-package-zoom",subdir="zoom")
# install.packages("tidyverse")
# install.packages("htmlwidgets")
# install.packages("ggplotly")
#install.packages("ggpubr") #to put the equation in the plot
#-------------- Loading packages -------------------
# Signal processing library
library(spant)
# Graphs libraries
library(plotly)
library(ggplot2)
library(gridExtra)
# for plots w ggplot
library(tidyverse)
library(htmlwidgets)
library(ggplot2)
library(ggpubr)
# library(zoom)
#-------------- Setting paths --------------------

# to go to the scripts folder
# In Mac
setwd("/Users/alecastillab/Library/CloudStorage/OneDrive-UniversityofToronto/Research/organoids/hESC_LT/scripts")
# In SRI Desktop
setwd("C:/Users/jnear/OneDrive - University of Toronto/Research/organoids/hESC_LT/scripts")
#load the R.Data in Desktop
# load("C:/Users/jnear/OneDrive - University of Toronto/Research/organoids/hESC_LT/scripts/.RData")
# In surface
setwd("C:/Users/macas/OneDrive - University of Toronto/Research/organoids/hESC_LT/scripts")
# ------------- making metabolites----------------------------------------------------------------------------------
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
chem_shift_ala2 <- c(3.746, 1.453, 1.453, 1.453)
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
#-------------- organoid volumes----------------------------------------------- 
#organoid volumes in mm^3
d <- c(2.1, 2, 1.5, 1.5, 1.5, 1.8, 1.9, 2, 2.1, 2, 2, 2.5, 2, 2, 2.2, 2, 2, 2, 2)
vol_mm = (4/3)*pi*(d/2)^3
#print(vol_mm)

#conversion to litres
vol_l = vol_mm*(1/1000000)
#print(vol_l)

#-------------- reference fit ala28_lcm & ala29_lcm-----
# we use the same concentration and value of ala29 for ala28
#reference sample fit
#volume = 30ul
#concentration = 20mM

ref_list <- list(ala2)
ref_basis <- sim_basis(ref_list, pul_seq = seq_pulse_acquire, acq_paras = def_acq_paras(N = 8124, fs = 10000, ft = 5.0028e8, ref = 4.65, nuc="1H"))

#refname28 <- "../alanineRef/ala28_lcm"
#ref_data_28 <- read_mrs(refname28, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref = 4.65)
refname29 <- "../alanineRef/ala29_lcm"
ref_data_29 <- read_mrs(refname29, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref = 4.65)

#-------------- abfit_opts default parameters ----------------------------------
# opts <- abfit_opts(
#                    pre_fit_bl_ed_pppm = 1, #allows the fit to move
#                    max_bl_ed_pppm = 7, #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                    min_bl_ed_pppm = NULL, #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                    auto_bl_flex_n = 20, #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                    algo_pre = "NLOPT_LN_NELDERMEAD", #NLOPT_LN_SBPLX !!!!!!!!!!!
#                    auto_bl_flex=TRUE,
#                    init_damping=5, # no deja 25
#-------------- other abfit_opts default parameters ----------------------------
                   # max_shift=0.078,
                   # max_damping=15,
                   # max_phase=360,
                   # lambda=NULL, #no deja ningun numero, es un smoothing parammetter
                   # ppm_left=4.8, #son parte de xlim
                   # ppm_right=0.2, #son parte de xlim
                   # maxiters=1024, #se revienta con el doble e incluso con 1800
                   # zp=TRUE,
                   # bl_ed_ppm=2, # baseline smoothness parameter (ED per ppm). want to be like 0
                   # bl_comps_pppm=2,
                   # auto_bl_flex_n = TRUE,
                   # export_sp_fit = FALSE,
                   # max_asym = 0.25,
                   # max_basis_shift = 0.0078,
                   # max_basis_damping = 2,
                   # maxiters_pre = 1000,
                   # remove_lip_mm_prefit = FALSE,
                   # pre_align = TRUE,
                   # max_pre_align_shift = 0.1,
                   # pre_align_ref_freqs = c(2.01, 3.03, 3.22),
                   # noise_region = c(-0.5, -2.5),
                   # optimal_smooth_criterion = "maic",
                   # aic_smoothing_factor = 5,
                   # anal_jac = TRUE,
                   # pre_fit_ppm_left = 4,
                   # pre_fit_ppm_right = 1.8,
                   # phi1_optim = FALSE,
                   # phi1_init = 0,
                   # max_dphi1 = 0.2,
                   # max_basis_shift_broad = 0.0078,
                   # max_basis_damping_broad = 2,
                   # ahat_calc_method = "lh_pnnls",
                   # prefit_phase_search = TRUE,
                   # freq_reg = NULL,
                   # output_all_paras = FALSE)
#default
#opts <- abfit_opts(ppm_left=4.8, pre_fit_bl_ed_pppm = 1, max_bl_ed_pppm = 7, algo_pre = "NLOPT_LN_NELDERMEAD", bl_ed_ppm=2,auto_bl_flex=TRUE, prefit_phase_search = TRUE,max_phase=360)
#                    
#testing
#opts <- abfit_opts(ppm_left=4.8, pre_fit_bl_ed_pppm = 6.8, max_bl_ed_pppm = 14, auto_bl_flex_n = 5, algo_pre = "NLOPT_LN_SBPLX", bl_ed_ppm=0)
#opts <- abfit_opts(ppm_left=4.8, pre_fit_bl_ed_pppm = 6.8, max_bl_ed_pppm = 14, auto_bl_flex_n = 5, algo_pre = "NLOPT_LN_SBPLX", auto_bl_flex=FALSE, prefit_phase_search = FALSE) #GOOD STARTING POINT THE BASELINE IS FLAT you dont want the fit to fit the noise
#-------------- abfit opts for ala28--------------------------------------------
#fit for ala28
opts28 <- abfit_opts(pre_fit_bl_ed_pppm = 2, algo_pre = "NLOPT_GN_ISRES")

ref_res28 <- fit_mrs(ref_data_28, ref_basis, opts=opts28)
plot(ref_res28, main= 'ref_res28')
plot(ref_res28, xlim = c(1.4,1.5), main= 'ref_res28' )
plot(ref_res28, xlim = c(3.7,3.8), main= 'ref_res28' )

#-------------- abfit opts for ala29 -------------------------------------------
opts29 <- abfit_opts(pre_fit_bl_ed_pppm = 6.8, max_bl_ed_pppm = 1, auto_bl_flex_n = 3, max_shift=0.1, algo_pre = "NLOPT_LN_COBYLA")

ref_res29 <- fit_mrs(ref_data_29, ref_basis, opts=opts29)

plot(ref_res29, main= 'ref_res29')
plot(ref_res29, xlim = c(1.4,1.55), main= 'ref_res29')
plot(ref_res29, xlim = c(3.7,3.8) )

# to change: integrate ggplot in ggplotly
ggplot(data = ref_res29$res_tab, aes (x = data ))
p <- ggplot(data = ref_res29, aes (x = data, fill = clarity) )
ggplotly(p)
#---------------loading organoid data--------------------------------------------------
#load metabolites and simulate basis set
met_list <- list(ace2, lac2, ala2, get_mol_paras("cho"), get_mol_paras("glc"), gln2, gly2, get_mol_paras("gpc"), get_mol_paras("ins"), get_mol_paras("naa"), get_mol_paras("pch"), get_mol_paras("gaba"), cr2, etoh, htau, get_mol_paras("glu"))
basis <- sim_basis(met_list, pul_seq = seq_pulse_acquire, acq_paras = def_acq_paras(N = 8124, fs = 10000, ft = 5.0028e8, ref = 4.65, nuc="1H"))
stackplot(basis, xlim = c(4.2, 1), y_offset = 50, labels = basis$names)
#load data files/spectra
folders_path <- "../" 
b <- paste(folders_path, "B/liqph_lcm", sep="")
c <- paste(folders_path, "C/liqph_lcm", sep="")
d <- paste(folders_path, "D/liqph_lcm", sep="")
e <- paste(folders_path, "E/liqph_lcm", sep="")
f <- paste(folders_path, "F/liqph_lcm", sep="")
g <- paste(folders_path, "G/liqph_lcm", sep="")
h <- paste(folders_path, "H/liqph_lcm", sep="")
i <- paste(folders_path, "I/liqph_lcm", sep="")
k <- paste(folders_path, "K/liqph_lcm", sep="")
l <- paste(folders_path, "L/liqph_lcm", sep="")
s <- paste(folders_path, "S/liqph_lcm", sep="")
t <- paste(folders_path, "T/liqph_lcm", sep="")
u <- paste(folders_path, "U/liqph_lcm", sep="")
aa <- paste(folders_path, "AA/liqph_lcm", sep="")
bb <- paste(folders_path, "BB/liqph_lcm", sep="")
cc <- paste(folders_path, "CC/liqph_lcm", sep="")
dd <- paste(folders_path, "DD/liqph_lcm", sep="")
ee <- paste(folders_path, "EE/liqph_lcm", sep="")
ff <- paste(folders_path, "FF/liqph_lcm", sep="")

#read files
files <- c(b, c, d, e, f, h, i, k, l, s, t, u, aa, bb, cc, dd, ee, ff)
data <- read_mrs(files, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref=4.65)

#optional: check data by plotting

for (spectra in data) {
  plot(spectra)}

#----------------------------------
#--------------- fit analyis of the organoid data -----------------------------
#----------------------------------
opts <- abfit_opts(pre_fit_bl_ed_pppm = 6.8, max_bl_ed_pppm = 1, auto_bl_flex_n = 3, 
                   max_shift=0.1, algo_pre = "NLOPT_LN_COBYLA")
fit_res <- fit_mrs(data, basis, opts=opts)

#plot fits 
for (fit in fit_res) {
  plot(fit)
}

# As an option: it is possible to make the fit individually
# load libraries, and make metabolites in line 37 and load organoid data in line 336 
opts_test <- abfit_opts(pre_fit_bl_ed_pppm = 6.8, max_bl_ed_pppm = 1, auto_bl_flex_n = 3, 
                   max_shift=0.1, algo_pre = "NLOPT_LN_COBYLA")
organoide <- ff
n_organoide <- str_sub(organoide, 4,5)
data_test <- read_mrs(organoide, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref = 4.65)
fit_res_test <- fit_mrs(data_test, basis, opts=opts_test)

plot(fit_res_test, main= paste('fit',  n_organoide, 'with', opts_test$algo_pre) )

# para cada organoide
#correr el fit en un loop que le pase el vector de los strings de diferentes NLOPT algoritmo 

#--------------- saving fits in a workspace -----------------------------
save(file = "../workspace_May18.RData")
save()

# you can run the fit or load it from this file:
#load("../hESC_LT/scripts/workspace_may17_surface.RData")

#------------------manually fit/anaylze dataset ----------------
# #B
# met_list <- list(get_mol_paras("lac"), ala2, get_mol_paras("cho"), get_mol_paras("glc"), gln2, gly2, get_mol_paras("gpc"), get_mol_paras("ins"), get_mol_paras("naa"), get_mol_paras("pch"), get_mol_paras("gaba"), cr2, custom_mol, htau)
# basis <- sim_basis(met_list, pul_seq = seq_pulse_acquire, acq_paras = def_acq_paras(N = 8124, fs = 10000, ft = 5.0028e8, ref = 4.65, nuc="1H"))
# stackplot(basis, xlim = c(4.8, 0.2), y_offset = 50, labels = basis$names)
# b <- "~/Documents/spant/B"
# b_data <- read_mrs(b, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref=4.65)
# opts <- abfit_opts(ppm_left=4.8)
# fit_res_b <- fit_mrs(b_data, basis, opts=opts)
# plot(fit_res_b)

# input for fit mrs: lcm file of each organoid
# organoid_b <- "../B/liqph_lcm"
# ref_data_b <- read_mrs(organoid_b, format = "lcm_raw", ft=5.0028e8, fs = 10000, ref = 4.65)
# 
# opts_b <- abfit_opts(pre_fit_bl_ed_pppm = 15, algo_pre = "NLOPT_LN_NEWUOA")
# ref_res_b <- fit_mrs(ref_data_b, ref_basis, opts=opts_b2)
# plot(ref_res_b, main= 'organoid_b')

# opts_b2 <- abfit_opts(pre_fit_bl_ed_pppm = 2, 
#                      max_bl_ed_pppm = 3,
#                      algo_pre = "NLOPT_LN_NEWUOA")
# 
# ref_res_b2 <- fit_mrs(ref_data_b, ref_basis, opts=opts_b2)
# plot(ref_res_b, main= 'organoid_b: pre2-max3')
# 
# plot(ref_res_b, xlim = c(1.4,1.5), main= 'ref_res_b' )
# plot(ref_res_b, xlim = c(3.7,3.8), main= 'ref_res_b' )
#---------------quantify concentrations ala29 and the volume of organoids-------
output <- matrix(nrow = length(fit_res), ncol = 16)
# attempt 1
#index=1
#for (i in seq_along(fit_res)) {
#  metabs <- fit_res[[i]]$res_tab[, 6:21]
#  for (j in seq_along(metabs)) {
#    conc <- (0.02)*(metabs[[j]]/120629.2)*(0.00003/vol_l[[index]])
#    output[i, j] <- conc
#  }
#  index <- (index+1)
#}
ref_res <- ref_res29
index=1
for (i in seq_along(fit_res)) {
  metabs <- fit_res[[i]]$res_tab[, 6:21]
  for (j in seq_along(metabs)) {
    #conc <- (S_metabolite/S_alaref)*(conc_alaref*vol_alaref)/(Mass_org)
    # Mass_org is suposing a density of 1 
    conc <- (metabs[[j]]/ref_res$res_tab$Ala) * (20*0.00003/vol_l[[index]])
    
    output[i, j] <- conc
  }
  index <- (index+1)
}
#---------------saving table of [metabs]  by volume ----------- -----------
#write vector as data frame
df <- as.data.frame(output)
#add age column
df$Age <- c(103, 104, 121, 122, 123, 137, 138, 185, 184, 311, 312, 313, 85, 86, 87, 144, 145, 146)
#rewrite/rename column names
colnames(df) <- c("Ace", "Lac", "Ala", "Cho", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")
rownames(df)  <- c("b", "c", "d", "e", "f", "h", "i", "k", "l", "s", "t", "u", "aa", "bb", "cc", "dd", "ee", "ff")
#write data frame to csv file
# in mMols!!!
write.csv(df,"../csv/conc-byVol-May18_wAla29.csv", row.names = TRUE)

#---------------quantify concentrations by mass of organoids-------
# mass of organoids in grames
# m_org <- c(0.0113, 0.0158, 0.0091, 0.0039, 0.0098, 0.0043, 0.0076, 0.0089, 0.0038)
# output_mass <- matrix(nrow = length(fit_res[10:18]), ncol = 16)

Kg_org <- c(0.0000113, 0.0000158, 0.0000091, 0.000039, 0.0000098, 0.0000043, 0.0000076, 0.0000089, 0.0000038)
output_Kmass <- matrix(nrow = length(fit_res[10:18]), ncol = 16)
ref_res <- ref_res29

index=1
for (i in seq_along(fit_res[10:18])) {
  metabs <- fit_res[[i]]$res_tab[, 6:21]
  for (j in seq_along(metabs)) {
    #conc <- (S_metabolite/S_alaref)*(conc_alaref*vol_alaref)/(Mass_org)
    # we should make sure that 30 uL is in L!!! then we can cancel
    conc <- (metabs[[j]]/ref_res$res_tab$Ala) * (20*0.00003/Kg_org[[index]])
    
    output_Kmass[i, j] <- conc
  }
  index <- (index+1)
}
#write vector as data frame
df_mass <- as.data.frame(output_Kmass)
#add age column
df_mass$Age <- c(311, 312, 313, 85, 86, 87, 144, 145, 146)
#rewrite/rename column names
colnames(df_mass) <- c("Ace", "Lac", "Ala", "Cho", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")
rownames(df_mass)  <- c("s", "t", "u", "aa", "bb", "cc", "dd", "ee", "ff")
#write data frame to csv file
# in mMols!!!
write.csv(df_mass,"../csv/conc-bymass-May24_wAla29-Kmass.csv", row.names = TRUE)

#---------------plots -------------------------------------------------------------------------------------------------------------------
#plot the fit for a single metabolite
plot(fit_res[[1]]$fits$`1.fit`$PPMScale, fit_res[[1]]$fits$`1.fit`$Gln, type="l")
#---------------mean standard deviation for all metabolites --------------------------
for (fit in fit_res) {
  meansd = sum(fit$res_tab[,22:37])/length(fit$res_tab[,22:37])
  print(meansd)
}
#creating a list where the plots will be stored
plot_list = list()
# df <- as.data.frame(cleanoutput)
# df$Age <- c(103, 104, 121, 122, 123, 136, 137, 138, 185, 184, 311, 312, 313, 85, 86, 87, 144, 145, 146)
# colnames(df) <- c("Ace", "Lac", "Ala", "Cho", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")
xvar <- df$Age

for (i in 1:16) {
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

#plot graphs side by side with library gridExtra
out <- do.call(grid.arrange, plot_list)
plot(out)

#do.call: calls a function and a list and passes the elements of the list 
#as arguments for the function
#in matlab, i think the "equivalent" would be listing the elements in curly brackets
#but i couldn't do that in r

#--------------statistical analysis -------------------------------------------- -----

#dataset with  metabolites with 0 concentration removed (choline)

# stats_output <- cleanoutput[,-c(4)]
stats_output <- df
statsdf <- as.data.frame(stats_output)
statsdf$Age <- c(103, 104, 121, 122, 123, 137, 138, 185, 184, 311, 312, 313, 85, 86, 87, 144, 145, 146)
colnames(statsdf) <- c("Ace", "Lac", "Ala", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")
Age <- statsdf$Age

tempdf <- data.frame(Age, Concentration)
for (i in 1:15) {
  Concentration <- statsdf[, i]
  ttest = t.test(Concentration ~ Age, data = tempdf)
  print(ttest)
}

##-----------------------------------------------------------------------------
#-------------- Declaring integral vector -------------------
integrals_scan3 <- c(6528424.54934505,	5314296.16059844,	24837035.5865943,	1099090.19830726,	2479930.32369288,	1209053.96561301,	2333201.41742437,	1036917.02000616,	2345479.47206918)
#-------------- Making a data frame of integrals------------
df_integrals_scan3 <- as.data.frame(integrals_scan3)
# adding mass column
df_integrals_scan3$m_org <- c(0.0113, 0.0158, 0.0091, 0.0039, 0.0098, 0.0043, 0.0076, 0.0089, 0.0038)

# plotting results w outlier
plot(x = df_integrals_scan3$m_org, y= df_integrals_scan3$integrals_scan3)

# clean data, by deleting third row (it seems like an outlier)
cleanoutput <- df_integrals_scan3[-c(3),]

# plotting results w/o outlier plot(x,y)
plot(cleanoutput$m_org, cleanoutput$integrals_scan3)

#-------------- Fitting a model m_org~integrals ------------
lm_mass = lm(m_org~integrals_scan3, data = cleanoutput) #Create the linear regression for cleanoutput
summary(lm_mass) #Review the results
# Note: the higher the t-statistic (and the lower the p-value), the more significant the predictor. 
#----------------plotting ----------------------------------
regression <- ggplot(cleanoutput, aes(x = m_org, y = integrals_scan3)) + 
  geom_point()
# adding the line to the plot
regression <- regression + geom_smooth(method="lm", col="black")

# adding the equation to the plot
regression <- regression + stat_regline_equation(label.x = 0.01, label.y = 8e6)
regression

#-------------- Making predictions -----------------------
model <- lm_mass
predictions <- model %>% predict(cleanoutput)
# # Model performance
# # (a) Prediction error, RMSE
# RMSE(predictions, cleanoutput$m_org)
# # (b) R-square
# R2(predictions, cleanoutput$m_org)

integrals_scan3_missingwight <- c(3361433.67172262,	1035671.79148143,	1094674.48911844,	1034696.27058579,	1286596.44318902,	1468438.17280932,	1284491.58874730,	2171514.84000100,	1287818.60650576);
new_integrals_scan3 <- data.frame(integrals_scan3_missingwight)
rownames(new_integrals_scan3) <- c("B", "C", "D", "E", "F", "H", "I", "K", "L")

# Predicts the mass values
new_mass <- predict(model, newdata = new_integrals_scan3)

# Adding the weight on side of the integrals
new_integrals_scan3$Weight <- new_mass

# ------------- Adding new m_org to the vector -------------------------
mass_g <- c(df_integrals_scan3$m_org, new_integrals_scan3$Weight)
mass_Kg <- mass_g/1000
# ------------- quantify concentrations by Kg of organoids, including new weights -------------------------
# if it were like just 9 conc_by_weightKg <- matrix(nrow = length(fit_res[10:18]), ncol = 16)
# for (i in seq_along(fit_res[10:18])) { ... 
conc_by_weightKg <- matrix(nrow = length(fit_res), ncol = 16)

#ref_res <- ref_res29
index=1
for (i in seq_along(fit_res)) {
  metabs <- fit_res[[i]]$res_tab[, 6:21]
  for (j in seq_along(metabs)) {
    #conc <- (S_metabolite/S_alaref)*(conc_alaref*vol_alaref)/(Mass_org)
    # Mass_org is suposing a density of 1 
    conc <- (metabs[[j]]/ref_res$res_tab$Ala) * (20*0.00003/mass_Kg[[index]])
    
    conc_by_weightKg[i, j] <- conc
  }
  index <- (index+1)
}
#write vector as data frame
df_Kg <- as.data.frame(conc_by_weightKg)
#add age column
df_Kg$Age <- c(103, 104, 121, 122, 123, 137, 138, 185, 184, 311, 312, 313, 85, 86, 87, 144, 145, 146)
#rewrite/rename column names
colnames(df_Kg) <- c("Ace", "Lac", "Ala", "Cho", "Glc", "Gln", "Gly", "GPC", "Ins", "NAA", "PCh", "GABA", "Cr", "EtOH", "hTau", "tCho", "Age")
rownames(df_Kg)  <- c("b", "c", "d", "e", "f", "h", "i", "k", "l", "s", "t", "u", "aa", "bb", "cc", "dd", "ee", "ff")
#write data frame to csv file # in mMols!!!
write.csv(df_Kg,"../csv/conc-byKg-May24_wAla29_18organoids.csv", row.names = TRUE)
