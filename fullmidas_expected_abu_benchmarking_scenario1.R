


pacman::p_load(
  tidyverse, # data cleaning
  janitor, # data cleaning
  rio, # importing data
  here, # relative file paths
  Hmisc, # data exploration
  skimr, # data exploration
  VIM, # missing data check
  gtsummary, # tables
  flextable, # tables 
  RColorBrewer, # ggplot Colour Brewer
  cowplot, # plots together
  patchwork,
  ggrepel,
  glue) # pasting estimates


rm(list = ls())



#######
# Script for estimating expected DDD of antibiotics 
# Total DID using benchmarking then A/W/R


# Full MIDAS data - scenario 1 - Primary analysis  
#######

# Total DDD
# First step is to define benchmark by cluster 
# Then using benchmark DID get ratio of DID to total infection incidence using matrix of total infe ction & imputations (where relevant)
# Then taking this ratio matrix and multiplying by matrix of total infection for all countries in each benchmark group 
# That gives you matrix of expected DID then need to convert to expected DDD using population total 
# Estimate Reserve matrix
# Estimate Watch matrix 
# Using matrix of total DDD expected, can subtract matrices of Reserve, Watch to get Access 
#####


# Scenario 1 details: this looks to establish the primary analysis scenario of expected total and AWaRe precsribing
# Cases account for adjusted UTI and cellulitis and uses dysentary as the primary diarrhoea cases + a few additional infections for HAP and IAI
# Total infection incidence uses these adjusted case counts 
# Use Watch scenario for typhoid where it is inflated >50 per 100k, 20.7% upper UTI, 15% severe LRI in HIC + additional Watch for HAI (for HAP, IAI and NF-SSTI)
# this is the main dataset & same benchmark countries (as scenarios 2 and 3) from MIDAS
# Typhoid is inflated 10x in the primary scenario

# *** UPDATE - 01 APRIL 2026 - updated for response to reviewers for LPH; removed XDR TB from Watch ***

#-------------------------------#
# Read in datasets #####
#-------------------------------#

# import full MIDAS cluster assignments 
clusters <- import(here("04_final_datasets", "clustering", "fullmidas", "fullmidas.cluster.assign.tidy.RDATA"))

# import covariates (post processing)
covars <- import(here("04_final_datasets", "post_imputation", "covariates.postimp.RDATA"))

# select descriptive vars
covars <- covars %>% dplyr::select(iso3_code, country_name, who_reg, midas_country_flag, glass_country_flag, pop_totl, gni_cap_usdat, income)

# Read in Full MIDAS data post imputation - complete imputed dataset - select DID draws for benchmark countries 
fullmidas <- import(here("04_final_datasets", "post_imputation", "midasfull.postimp.RDATA"))

# analysis countries 
analysis.countries.df <- import(here("04_final_datasets", "analysis.countries.RDATA"))
analysis.countries <- analysis.countries.df$iso3_code # want as a vector of values

# Read in analysis functions ####
source(here("05_analysis", "fair_share_paper", "abu_estimate_functions.R"))



#---------------------------------------------# 
# Scenario specific datasets & settings ####
#---------------------------------------------#

# infection incidence draws when including the estimates of HAIs that weren't in GBD 
total.inf.inc <- import(here("04_final_datasets", "infection_draws", "total.inc.per100k.draws.scen1.RDATA"))


# Watch infection case counts 

# read in case counts of variables needed for watch DDD - cases as is and adjusted cases 
dysentery.cases <- import(here("04_final_datasets", "infection_draws", "dys.cases.df.RDATA"))
cholera.cases <- import(here("04_final_datasets", "infection_draws", "cholera.cases.draws.RDATA"))
chlamydia.cases <- import(here("04_final_datasets", "infection_draws", "chlamydia.cases.df.RDATA"))
gonorrhoea.cases <- import(here("04_final_datasets", "infection_draws", "gonorrhoea.cases.df.RDATA"))
xdr.tb.cases <- import(here("04_final_datasets", "infection_draws", "xdr_tb.cases.df.RDATA"))
mdr.tb.cases <- import(here("04_final_datasets", "infection_draws", "mdr_tb.cases.df.RDATA"))
sepsis.cases <- import(here("04_final_datasets", "infection_draws", "sepsis.cases.sim.df.RDATA")) %>% rename(iso3_code = country) %>% rename_with(~ sub("^V", "X", .x), .cols = starts_with("V"))
uti.upper.cases <- import(here("04_final_datasets", "infection_draws", "uti.predicted.upper.cases.RDATA")) %>% dplyr::select(-pop_totl)
severe.lri.cases <- import(here("04_final_datasets", "infection_draws", "severe.lri.cases.RDATA"))
typhoid.cases <- import(here("04_final_datasets", "infection_draws", "typhoid.10xsuspected.50per100k.cases.RDATA")) # changed to teh 10x inflated 

ssti.nf.cases <- import(here("04_final_datasets", "infection_draws", "ssti.nf.cases.RDATA")) %>% rename_with(~ sub("^V", "X", .x), .cols = starts_with("V"))
iai.cases <- import(here("04_final_datasets", "infection_draws", "iai.cases.RDATA")) %>% rename_with(~ sub("^V", "X", .x), .cols = starts_with("V"))
hap.cases <- import(here("04_final_datasets", "infection_draws", "hap.cases.lriratio.RDATA")) #%>% rename_with(~ sub("^V", "X", .x), .cols = starts_with("V"))


# specify scenario 
# scenario.folder <- "scenario_1"
# scenario.file <- "scen1"
scenario.folder <- "scenario_1_revision" # make new folder for the revised scenario w/o XDR TB
scenario.file <- "scen1"

# specify benchmark countries 
benchmark.country.c4 <- "CHE"
benchmark.country.c123 <- "MAR"



# ----------# 
# Imputed data processing into matrix ###
#-----------#

# === Filter to MIDAS countries only ==== # 

fullmidas <- fullmidas %>% filter(midas_country_flag == "midas" & .imp > 0)

# select total DID 

totaldid.midas <- fullmidas %>% dplyr::select(c(iso3_code, .imp, did_total))


# === Convert imputed dataset to wide format for total DID ==== # 

totaldid.midas.w <- totaldid.midas %>%
  pivot_wider(id_cols = "iso3_code",
              names_from = c(".imp"),
              values_from = c("did_total")) %>%
  rename_with(~ paste0("X", .x), .cols = -c(iso3_code))




# --------------- #
# Define benchmark countries & select imputations matrix ####
#-----------------# 

# Benchmark coutnries selected - see table / text
# Cluster 4 (HIC) = Switzerland
# Cluster 3 = Morocco
# Cluster 1 & 2 - no suitable benchmark so using cluster 3 benchmark = Morocco


# Benchmark countries
benchmark.did.c4 <- totaldid.midas.w %>% filter(iso3_code == benchmark.country.c4)
benchmark.did.c123 <- totaldid.midas.w %>% filter(iso3_code == benchmark.country.c123)



# select infection incidence matrices for benchmark countries separately 
total.inf.c4 <- total.inf.inc %>% filter(iso3_code == benchmark.country.c4)
total.inf.c123 <- total.inf.inc %>% filter(iso3_code == benchmark.country.c123)



#-------------------------------#
# calculate ratio of DID to infection incidence #####
#-------------------------------#
# do for all 1000 draws for both benchmark countries 
# do this element wise (e.g. x1 to x1)


# Run function using the infection incidence (withouth HAI)

ratio.c4 <- divide_draws_f(benchmark.did.c4, total.inf.c4)

ratio.c123 <- divide_draws_f(benchmark.did.c123, total.inf.c123)



#-----------------------#
# join cluster assignments to the matrix of total infection incidence & then split 
#------------------------#

total.inf.inc <- total.inf.inc %>% left_join(clusters %>% dplyr::select(c(iso3_code, cluster))) %>%
  relocate(cluster, .after = iso3_code)


# === split into two datasets for two difference benchmarks === #

total.inf.inc.c4 <- total.inf.inc %>% filter(cluster == "Cluster 4") %>% dplyr::select(-c(cluster))

total.inf.inc.c123 <- total.inf.inc %>% filter(cluster != "Cluster 4") %>% dplyr::select(-c(cluster))




#--------------------------------#
# multiply relevant ratio matrix by total infection incidence matrix for each country by cluster assignment 
#---------------------------------#

# === Multiply for Cluster 4 ==== # 
# for switzerland should return the original DID 

exp.total.c4 <- multiply_draws_matrix_f(total.inf.inc.c4, ratio.c4)


# === Multiply for Clusters 1, 2, 3 === # 
# for morocco should return the original did per row 
exp.total.c123 <- multiply_draws_matrix_f(total.inf.inc.c123, ratio.c123)


#------------------------------------------#
# convert matrix of total estimated DID to DDD 
#------------------------------------------#

# combine both groups
exp.total <- bind_rows(exp.total.c4, exp.total.c123) %>% arrange(iso3_code)


# join with population data 
exp.total <- exp.total %>% left_join(covars %>% dplyr::select(c(iso3_code, pop_totl))) %>% relocate(pop_totl, .after = "iso3_code")

# convert to DDD
exp.total.ddd <- exp.total %>%
  mutate(across(starts_with("X"),  ~ (. * (pop_totl * 365)) / 1000)) %>%
  dplyr::select(-c(pop_totl))


#-------------------------------------#
# Estimate Reserve DDD needed ####
#-------------------------------------#

# ==== Read in datasets ==== # 

# Reserve DDD from benchmark countries

reserve.benchmark <- import(here("04_final_datasets", "reserve.ddd.benchmark.RDATA"))

# Read in summed CRO, VRO and TB cases 
cro.cases <- import(here("04_final_datasets", "resistance_draws", "cro.cases.df.RDATA"))
vro.cases <- import(here("04_final_datasets", "resistance_draws", "vro.cases.df.RDATA"))
mdrtb.cases <- import(here("04_final_datasets", "infection_draws", "mdr_tb.cases.df.RDATA"))
xdrtb.cases <- import(here("04_final_datasets", "infection_draws", "xdr_tb.cases.df.RDATA"))

# actual case counts and ddds not cases/100k and DID 

# Filter cases to just the analysis countries 
cro.cases <- cro.cases %>% filter(iso3_code %in% analysis.countries)
vro.cases <- vro.cases %>% filter(iso3_code %in% analysis.countries)
mdrtb.cases <- mdrtb.cases %>% filter(iso3_code %in% analysis.countries)
xdrtb.cases <- xdrtb.cases %>% filter(iso3_code %in% analysis.countries)

# === Sum up datasets for VRO and MDR/XDR TB === # 

# List to store summed data
#summed_cases_midas <- list()

# List of the datasets to iterate over
reserve_midas_cases_datasets <- list(vro.cases, mdrtb.cases, xdrtb.cases) # not CRO 


#prep data to be summed - ensure correct order etc 
summed_cases_midas <- summing_data_prep_f(datasets = reserve_midas_cases_datasets, 
                                          analysis_countries= analysis.countries, 
                                          dataset_names = NULL)



# Sum all numeric columns across datasets
summed_midas_cases <- Reduce("+", summed_cases_midas)

# Return a final dataframe with summed draws and iso3_code column
summed.gpos.cases <- data.frame(iso3_code = analysis.countries, summed_midas_cases)



# ---------------------------------------#
# Estimate Reserve using benchmark ###
#----------------------------------------#
# Do this in two parts - GPB/TB and GNB DDDs then sum 


# === CRO cases === # 

# CRO cases from Switzerland
midas.cro.benchmark.country.c4 <- cro.cases %>% filter(iso3_code == benchmark.country.c4)

# get GNB DDD from Switzerland
midas.gnb.ddd.benchmark.country.c4 <- reserve.benchmark$reserve_ddd[reserve.benchmark$dataset == "full_midas" & reserve.benchmark$abx_cov_group == "gram_neg" & reserve.benchmark$iso3_code == benchmark.country.c4]

# Ratio of CRO to GNB coverage 
midas.gnb.ratio.benchmark.country.c4 <- midas.cro.benchmark.country.c4 %>%
  mutate(across(-iso3_code, ~ midas.gnb.ddd.benchmark.country.c4 / .))

# use multiply matrix function 
reserve.gnb.ddd.exp <- multiply_draws_matrix_f(cro.cases, midas.gnb.ratio.benchmark.country.c4)


# === VRO & TB === #

# VRO & TB cases from switzerland
midas.vrotb.benchmark.country.c4 <- summed.gpos.cases %>% filter(iso3_code == benchmark.country.c4)

# get GPB & TB DDD from switzerland
midas.vrotb.ddd.benchmark.country.c4 <- reserve.benchmark$reserve_ddd[reserve.benchmark$dataset == "full_midas" & reserve.benchmark$abx_cov_group == "gram_pos_tb" & reserve.benchmark$iso3_code == benchmark.country.c4] 

# Ratio of VRO to GPB coverage 
midas.vrotb.ratio.benchmark.country.c4 <- midas.vrotb.benchmark.country.c4 %>%
  mutate(across(-iso3_code, ~ midas.vrotb.ddd.benchmark.country.c4 / .))


# use multiply matrix function
reserve.gpb.ddd.exp <- multiply_draws_matrix_f(summed.gpos.cases, midas.vrotb.ratio.benchmark.country.c4)


# == Sum up both expected DDD === # 

# make sure both are arranged in iso3_code order 
reserve.gpb.ddd.exp <- reserve.gpb.ddd.exp %>% arrange(iso3_code)
reserve.gnb.ddd.exp <- reserve.gnb.ddd.exp %>% arrange(iso3_code)

# check order is actually matching
if (!all(reserve.gnb.ddd.exp$iso3_code == reserve.gpb.ddd.exp$iso3_code)) {
  stop("Ordering mismatch after sorting — check for duplicates or NA iso3_codes.")
}

# draws columns 
draws.cols <- setdiff(names(reserve.gnb.ddd.exp), "iso3_code")

# Sum only numeric columns (excluding iso3_code)
reserve.ddd.exp <- data.frame(
  iso3_code = reserve.gnb.ddd.exp$iso3_code,
  reserve.gnb.ddd.exp[draws.cols] + reserve.gpb.ddd.exp[draws.cols]
)

names(reserve.ddd.exp)




#------------------------------------#
# Watch antibiotic use estimates ####
#-------------------------------------#

# this is for the primary scenario which includes: 
# - inflated typhoid cases 10x (treating suspected typhoid) in countries with >50 per 100k cases 
# - growth-model adjusted UTI and cellulitis cases 
# - 20% UTI cases are upper UTI 
# - 15% LRI cases are severe LRI 
# - Dysentery estimates 
# Use the additional HAI cases that have been estimated and DDDs for these 

# Update 01 April 2026 - excluded XDR TB from scenario

# ==== Calculate Watch DDD for each dataset ==== # 

# watch DDD per treatment course from AWaRe Book using more conservative estimate 

dysentery.ddd <- dysentery.cases %>% mutate(across(-iso3_code, ~ . * 4.2))
cholera.ddd <- cholera.cases %>% mutate(across(-iso3_code, ~ . * 3.3))
chlamydia.ddd <- chlamydia.cases %>% mutate(across(-iso3_code, ~ . *3.3))
gonorrhoea.ddd <- gonorrhoea.cases %>% mutate(across(-iso3_code, ~ . * 3.425))
#xdr.tb.ddd <- xdr.tb.cases %>% mutate(across(-iso3_code, ~ . * 1218))
mdr.tb.ddd <- mdr.tb.cases %>% mutate(across(-iso3_code, ~ . * 364))
sepsis.ddd <- sepsis.cases %>% mutate(across(-iso3_code, ~ . * 10.5))

typhoid.ddd <- typhoid.cases %>% mutate(across(-iso3_code, ~ . * 10)) # cipro and ceftriaxone are this 
upper.uti.ddd <- uti.upper.cases %>% mutate(across(-iso3_code, ~ . * 7))
severe.lri.ddd <- severe.lri.cases %>% mutate(across(-iso3_code, ~ . * 12.5)) # this is for CAP Cefotax + Clarithromycin

ssti.nf.ddd <- ssti.nf.cases %>% mutate(across(-iso3_code, ~ . * 27))
iai.ddd <- iai.cases %>% mutate(across(-iso3_code, ~ . * 7.5))
hap.ddd <- hap.cases %>% mutate(across(-iso3_code, ~ . * 10.5))


# ==== Restrict all to the analysis countries & order ==== # 

# === For Watch w/ HAI extra === # 

#specify datasets to sum/filter
watch.datasets <- list(dysentery.ddd, typhoid.ddd, cholera.ddd, chlamydia.ddd, gonorrhoea.ddd, #xdr.tb.ddd, 
                       mdr.tb.ddd, sepsis.ddd, upper.uti.ddd, severe.lri.ddd, ssti.nf.ddd, iai.ddd, hap.ddd)


watch.sum.list <- summing_data_prep_f(datasets = watch.datasets,
                                      analysis_countries= analysis.countries, 
                                      dataset_names = NULL)


# Sum up datasets element wise 
watch.exp.ddd <- Reduce(`+`, watch.sum.list)

# Add back iso3_code from list
watch.exp.ddd <- cbind(iso3_code = analysis.countries, watch.exp.ddd)


# Convert this amount to DID 
watch.exp.did <- watch.exp.ddd %>%
  left_join(covars %>% dplyr::select(c(iso3_code, pop_totl))) %>%
  mutate(across(starts_with("X"), ~((. * 1000)/(pop_totl * 365))))



#-----------------------------------#
# Estimate Access DDD ####
#-----------------------------------#

# subtract Watch and Reserve from Total DDD to get Access DDD - element wise again 

# examine datasets 
str(exp.total.ddd)
str(watch.exp.ddd)
str(reserve.ddd.exp)

# Make sure they are in the same order by iso3_code and have same number of columns

exp.total.ddd <- exp.total.ddd %>% arrange(match(iso3_code, analysis.countries)) 
watch.exp.ddd <- watch.exp.ddd %>% arrange(match(iso3_code, analysis.countries))
reserve.ddd.exp <- reserve.ddd.exp %>% arrange(match(iso3_code, analysis.countries))


# specify column names for draws columns 
draws.cols <- setdiff(names(exp.total.ddd), c("iso3_code"))


# Subtract Watch and Reserve from total using draws 
# using case counts for Watch with the HAI 
access.exp.ddd <- data.frame(
  iso3_code = exp.total.ddd$iso3_code,
  exp.total.ddd[draws.cols] - reserve.ddd.exp[draws.cols] - watch.exp.ddd[draws.cols]
)



exp.total.ddd[1,2] - reserve.ddd.exp[1,2] - watch.exp.ddd[1,2]
access.exp.ddd[1,2]

exp.total.ddd[150,100] - reserve.ddd.exp[150,100] - watch.exp.ddd[150,100]
access.exp.ddd[150,100]



#-----------------------------#
# Calculate access percent DDD by country ####
#-----------------------------#

# calculate access % - no excess Watch 
draws.cols <- setdiff(names(access.exp.ddd), "iso3_code")

access.ddd.perc <- data.frame(
  iso3_code = exp.total.ddd$iso3_code,
  access.exp.ddd[draws.cols] / exp.total.ddd[draws.cols])




#----------------------------------------#
# Create dataframe of global totals DDDs ####
#-----------------------------------------#

# total ddd

global.total.ddd <- summarise_global_stat_f(data = exp.total.ddd, stat_fn = sum, "total_ddd")

# reserve ddd
global.reserve.ddd <- summarise_global_stat_f(reserve.ddd.exp, sum, "reserve_ddd")

# watch ddd 
global.watch.ddd <- summarise_global_stat_f(watch.exp.ddd, sum, "watch_ddd")

# access ddd
global.access.ddd <- summarise_global_stat_f(access.exp.ddd, sum, "access_ddd")

# global access % 
total.cols <- setdiff(names(global.access.ddd), "sum_type")

global.access.perc <- data.frame(
  sum_type = "access_ddd_perc",
  global.access.ddd[total.cols] / global.total.ddd[total.cols]
)


# check global access %
global.access.perc %>% 
  rowwise() %>%
  mutate(draw_median = median(c_across(-sum_type)),
         draw_lower = quantile(c_across(-sum_type), probs = 0.025),
         draw_upper = quantile(c_across(-sum_type), probs = 0.975)) %>%
  ungroup() %>%
  dplyr::select(draw_median, draw_lower, draw_upper)




# Bind rows together into one global dataset 

global.exp.ddd <- bind_rows(global.total.ddd, global.reserve.ddd, global.watch.ddd, global.access.ddd, global.access.perc)




#----------------------------#
# Save DDDs draws datasets ####
#----------------------------#



export(global.exp.ddd, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("global.exp.ddd.", scenario.file, ".RDATA")))
export(exp.total.ddd, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("total.exp.ddd.", scenario.file, ".RDATA")))
export(reserve.ddd.exp, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("reserve.exp.ddd.", scenario.file, ".RDATA")))
export(watch.exp.ddd, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("watch.exp.ddd.", scenario.file, ".RDATA")))
export(access.exp.ddd, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("access.exp.ddd.", scenario.file, ".RDATA")))
export(access.ddd.perc, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("access.ddd.perc.", scenario.file, ".RDATA")))





#------------------------------#
# Convert country-level DDDs to DIDs ####
#------------------------------#

# total did already calculated

total.exp.did <- exp.total %>% arrange(match(iso3_code, analysis.countries)) %>% relocate(pop_totl, .after = X1000)

# AWaRe
reserve.exp.did <- did_convert_f(data = reserve.ddd.exp, covariate_data = covars, pop_column = "pop_totl")

watch.exp.did <- did_convert_f(watch.exp.ddd, covars, "pop_totl")

access.exp.did <- did_convert_f(access.exp.ddd, covars, "pop_totl")

# spot check totals 
# ethiopia
total.exp.did[56, 2]
reserve.exp.did[56, 2] + watch.exp.did[56,2] + access.exp.did[56,2] 


# switzerland
total.exp.did[160,3]
reserve.exp.did[160, 2] + watch.exp.did[160, 2] + access.exp.did[160,2]


# === Export DID datasets ==== # 

export(total.exp.did, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("total.exp.did.", scenario.file, ".RDATA")))
export(reserve.exp.did, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("reserve.exp.did.", scenario.file, ".RDATA")))
export(watch.exp.did, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("watch.exp.did.", scenario.file, ".RDATA")))
export(access.exp.did, here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0("access.exp.did.", scenario.file, ".RDATA")))




# === Create summary datasets === # 

assign(
  paste0("total.summary.", scenario.file),
  summarise_draws_f(
    data = total.exp.did,
    estimate_label = "Total DID",
    scenario_label = scenario.folder
  ),
  envir = .GlobalEnv
)


assign(paste0("total.summary.", scenario.file), summarise_draws_f(data = total.exp.did, estimate_label = "Total DID", scenario_label = scenario.folder), envir = .GlobalEnv)

assign(paste0("access.summary.", scenario.file), summarise_draws_f(data = access.exp.did, estimate_label = "Access DID", scenario_label = scenario.folder), envir = .GlobalEnv)

assign(paste0("watch.summary.", scenario.file), summarise_draws_f(data = watch.exp.did, estimate_label = "Watch DID", scenario_label = scenario.folder), envir = .GlobalEnv)

assign(paste0("reserve.summary.", scenario.file), summarise_draws_f(data = reserve.exp.did, estimate_label = "Reserve DID", scenario_label = scenario.folder), envir = .GlobalEnv)


# bind together 

assign(paste0(scenario.folder, ".summary"), bind_rows(get(paste0("total.summary.", scenario.file)), get(paste0("access.summary.", scenario.file)), get(paste0("watch.summary.", scenario.file)), get(paste0("reserve.summary.", scenario.file))), envir = .GlobalEnv)


# create pasted together summary columns & add cluster
assign(paste0(scenario.folder, ".summary"), get(paste0(scenario.folder, ".summary")) %>%
         mutate(median_CI = glue("{round(draw_median, 1)} ({round(draw_lower, 1)}, {round(draw_upper, 1)})")) %>%
         left_join(clusters %>% dplyr::select(c(iso3_code, cluster = ordered_imp_4group_assign_f))) %>%
         relocate(cluster, .after = iso3_code), envir = .GlobalEnv)


# save

export(get(paste0(scenario.folder, ".summary")), here("06_results", "fair_share_paper", "abu_estimates", scenario.folder, paste0(scenario.file, ".summary.RDATA")))
