## Propensity Score Matching with MatchIt
## by Jonas Brändli & Nelly Blindenbacher
## June 2023

# Set path to file path
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Install packages
install.packages("MatchIt")
install.packages("readxl")
install.packages("ggpubr")
install.packages("VIM")

# Load Libraries
library(dplyr)
library(ggplot2)
library(MatchIt)
library(readxl)
library(ggpubr)
library(kableExtra)
library(VIM)

# Source Files
source("analysisPlot.R") # distribution plot

# Load data
clinicaldata_unprocessed <- read_xlsx("Originaldata_incl_imaging.xlsx")

# Limit data set to patients with CTA and CTP
clinicaldata_CT <- clinicaldata_unprocessed %>%
  filter(CTA == "1") %>%
  filter(CTP == "1")

# Set variable type
clinicaldata <- clinicaldata_CT %>%
  mutate(IAT = as.factor(i_iatrt),
         sex = as.factor(r_gender),
         age = as.numeric(age),
         IVT = as.factor(r_ivtrom),
         INR = as.numeric(linr_abl), 
         serum_creatinin = as.numeric(lkreat_abl), 
         systolic_blood_pressure = as.numeric(rrsyst_abl), 
         diastolic_blood_pressure = as.numeric(rrdias_abl), 
         previous_stroke = as.factor(b_pvstr), 
         diabetes_mellitus = as.factor(b_pvdm),
         hypertension = as.factor(b_pvrr),
         atrial_fibrillation = as.factor(b_pvaf),
         prestroke_mrs = as.factor(premrs),
         symptom_onset_to_door = as.numeric(dur_oa),
         NIHSS_baseline = as.factor(nihsco_abl_c),
         occlusion_site = as.factor(loc_cta_abl),
         collateral_score = as.factor(cgsc_cta_abl_c))

# Analyse missing data
variables_of_interest <- c("sex", "age", "IVT", "INR", "serum_creatinin", "systolic_blood_pressure", "diastolic_blood_pressure", "previous_stroke", "diabetes_mellitus", "hypertension", "atrial_fibrillation", "prestroke_mrs", "symptom_onset_to_door", "NIHSS_baseline", "occlusion_site", "collateral_score")
subset_data_of_interest <- clinicaldata[, variables_of_interest]

missing_values <- colSums(is.na(subset_data_of_interest))

variables_with_missing <- names(missing_values[missing_values > 0])
variables_with_missing

subset_missing_data <- clinicaldata[,variables_with_missing]

a <- aggr(subset_missing_data, plot = FALSE)
plot(a, numbers = TRUE, prop = FALSE, cex.axis = 0.4)


# Visualize Data Distribution
plotColumns0 <- c("sex", "age", "IVT", "INR", "serum_creatinin", 
                  "systolic_blood_pressure", "diastolic_blood_pressure", 
                  "previous_stroke", "diabetes_mellitus", "hypertension", 
                  "atrial_fibrillation", "prestroke_mrs", 
                  "symptom_onset_to_door", "NIHSS_baseline", 
                  "occlusion_site", "collateral_score")
p0 <- lapply(plotColumns0, analysisPlot, df = clinicaldata, therapy="IAT")
ggarrange(plotlist = p0, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")


# PSM Matching
set.seed(1)
psmMatch <- matchit(
  IAT ~ sex + age + IVT + systolic_blood_pressure + 
    diastolic_blood_pressure + previous_stroke + 
    diabetes_mellitus + hypertension + 
    atrial_fibrillation + prestroke_mrs + 
    NIHSS_baseline,
  data = clinicaldata, 
  distance = "glm", # logistic regression
  link = "logit", # logit link function
  method = "nearest", # matches the two nearest patients in treated and untreated group
  m.order = "random",
  caliper = 0.2) # exludes "bad" matches

# Pre Matching
preMatchDF <- clinicaldata
preMatchDF$distance <- psmMatch$distance

# Post Matching
summary(psmMatch)
psmMatchDF <- match.data(psmMatch)
# 92 patients in the Control Group are unmatched and therefore omitted

# Propensity Score Plot
ggarrange(analysisPlot(preMatchDF, "distance", "IAT"),
          analysisPlot(psmMatchDF, "distance", "IAT"),
          nrow = 2)

## Plot Distribution (first column before matching, second after matching)
plotColumns <- c("distance", "sex", "age", "IVT", "INR", "serum_creatinin", 
                 "systolic_blood_pressure", "diastolic_blood_pressure", 
                 "previous_stroke", "diabetes_mellitus", "hypertension", 
                 "atrial_fibrillation", "prestroke_mrs", 
                 "symptom_onset_to_door", "NIHSS_baseline", 
                 "occlusion_site", "collateral_score")
p1 <- lapply(plotColumns, analysisPlot, df = preMatchDF, therapy = "IAT")
p2 <- lapply(plotColumns, analysisPlot, df = psmMatchDF, therapy = "IAT")
newOrder <- as.vector(matrix(1:(length(plotColumns)*2), nrow = 2, byrow = T))
ggarrange(plotlist = c(p1,p2)[newOrder], nrow = 7, ncol = 2,  
          common.legend = TRUE, legend = "bottom")

## Reduce data to 20 patients per group

# Extract the 20 treated patients with the smallest difference in propensity score
pairs1 <- psmMatchDF %>%
  group_by(subclass) %>%
  summarize(diff=abs(distance[1]-distance[2])) %>%
  slice_min(diff, n=20) %>%
  select(., subclass)

psmMatchDF_short1 <- psmMatchDF %>%
  filter(., subclass %in% unlist(pairs1))

# or select randomly 20 pairs
pairs2 <- sample(unique(psmMatchDF$subclass), 20)

psmMatchDF_short2 <- psmMatchDF %>%
  filter(., subclass %in% pairs2)


## Plot Distributions 
# first col before matching, 
# second col after matching, 
# third col reduced to 20 pairs, with smallest diff
# 4th col reduced to 20 pairs, random
p3 <- lapply(plotColumns, analysisPlot, df = psmMatchDF_short1, therapy = "IAT")
p4 <- lapply(plotColumns, analysisPlot, df = psmMatchDF_short2, therapy = "IAT")
newOrder <- as.vector(matrix(1:(length(plotColumns)*4), nrow = 4, byrow = T))
ggarrange(plotlist = c(p1,p2,p3,p4)[newOrder], nrow = 4, ncol = 4,  
          common.legend = TRUE, legend = "bottom")


# Display the 20 pairs with the smallest difference in propensity score
clinicaldata_selection <- psmMatchDF_short1 %>%
  select(`IAT` = i_iatrt,
         `age` = age,
         `sex` = r_gender,
         `IVT`= r_ivtrom,
         `INR` = linr_abl,
         `serum_creatinin` = lkreat_abl,
         `systolic_blood_pressure` = rrsyst_abl,
         `diastolic_blood_pressure` = rrdias_abl,
         `previous_stroke` = b_pvstr,
         `diabetes_mellitus` = b_pvdm,
         `hypertension` = b_pvrr,
         `atrial_fibrillation` = b_pvaf,
         `prestroke_mrs` = premrs,
         `symptom_onset_to_door` = dur_oa,
         `NIHSS_baseline` = nihsco_abl_c,
         `occlusion_site` = loc_cta_abl,
         `collateral_score` = cgsc_cta_abl_c) %>%
  mutate(`age` = round(`age`),
         `sex` = ifelse(`sex` == 0, "M", "F"),
         `IAT` = ifelse(`IAT` == 0, "No", "Yes"),
         `IVT` = ifelse(`IVT` == 0, "No", "Yes"),
         `previous_stroke`= ifelse(`previous_stroke` == 0, "No", "Yes"),
         `diabetes_mellitus`= ifelse(`diabetes_mellitus` == 0, "No", "Yes"),
         `hypertension`= ifelse(`hypertension` == 0, "No", "Yes"),
         `atrial_fibrillation`= ifelse(`atrial_fibrillation` == 0, "No", "Yes"), 
         `occlusion_site` = as.character(`occlusion_site`),
         `occlusion_site` = case_when(
           `occlusion_site` == 0 ~ "No occlusion present",
           `occlusion_site` == 1 ~ "ICA",
           `occlusion_site` == 2 ~ "ICA-T",
           `occlusion_site` == 3 ~ "M1",
           `occlusion_site` == 4 ~ "M2",
           `occlusion_site` == 5 ~ "A1",
           `occlusion_site` == 6 ~ "A2",
           TRUE ~ `occlusion_site`))

# Final table of 40 selected patients
kable(clinicaldata_selection, format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)
