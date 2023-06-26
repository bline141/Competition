## Propensity Score Matching with MatchIt
## by Jonas Brändli & Nelly Blindenbacher
## June 2023


## Set path to file path-------------------------------------------------------
getwd()
setwd("Projekte/9 CT Model vs. Neurologist/Competition/")


## General Dependencies--------------------------------------------------------
install.packages(c("dplyr", "ggplot", "MatchIt", "readxl", "ggpubr", "kableExtra"))
library(dplyr)
library(ggplot2)
library(MatchIt)
library(readxl)
library(ggpubr)
library(kableExtra)


## Load, filter, set type and select clinical data-----------------------------
clinicaldata <- read_xlsx("Originaldata_incl_imaging.xlsx") %>%
  filter(CTA == "1", CTP == "1") %>%
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
         NIHSS_baseline = as.numeric(nihsco_abl_c),
         occlusion_site = as.factor(loc_cta_abl),
         collateral_score = as.factor(cgsc_cta_abl_c),
         followid = as.numeric(followid)) %>%
  select(c("IAT", "sex", "age", "IVT", "INR", "serum_creatinin", "systolic_blood_pressure", "diastolic_blood_pressure", "previous_stroke", 
           "diabetes_mellitus", "hypertension", "atrial_fibrillation", "prestroke_mrs", "symptom_onset_to_door", "NIHSS_baseline", "occlusion_site", "collateral_score", "followid"))


## Analyse missing data---------------------------------------------------------
# Dependencies
install.packages("VIM")
library(VIM)

# Find missing values
anyNA(clinicaldata)
missing_values <- colSums(is.na(clinicaldata))
print(missing_values)

# Show missing values and check for systematic missings
clinicaldata_missing <- clinicaldata[,c("INR", "serum_creatinin", "symptom_onset_to_door", "collateral_score")]
na <- aggr(clinicaldata_missing, plot = FALSE)
plot(na, numbers = TRUE, prop = FALSE, cex.axis = 0.3) # data is not systematically missing


## Impute missing data with missForest------------------------------------------
# Dependencies 
install.packages("missForest")
library(missForest)

# Impute missing values
set.seed(7)
clinicaldata <- as.data.frame(clinicaldata)
variables_with_na <- c("INR", "serum_creatinin", "symptom_onset_to_door", "collateral_score")
clinicaldata_imp <- missForest(clinicaldata[, variables_with_na], 
                               xtrue = clinicaldata[, variables_with_na], 
                               verbose = TRUE, maxiter = 50, ntree = 500, mtry = 3)
clinicaldata_imp$OOBerror

# Save, round and import imputed values
clinicaldata_imputed <- clinicaldata_imp$ximp #Save as dataframe
numeric_vars <- sapply(clinicaldata_imputed, is.numeric) # Identify numeric variables
clinicaldata_imputed[numeric_vars] <- round(clinicaldata_imputed[numeric_vars], 0) # Round numeric variables

# Import imputed values into clinical data set
clinicaldata$INR <- clinicaldata_imputed$INR
clinicaldata$serum_creatinin <- clinicaldata_imputed$serum_creatinin
clinicaldata$symptom_onset_to_door <- clinicaldata_imputed$symptom_onset_to_door
clinicaldata$collateral_score <- clinicaldata_imputed$collateral_score

## Source Files---------------------------------------------------------------
source("analysisPlot.R") # distribution plot

# Visualize Data Distribution---------------------------------------------------
plotColumns0 <- c("IAT", "sex", "age", "IVT", "INR", "serum_creatinin", "systolic_blood_pressure", "diastolic_blood_pressure", "previous_stroke", 
                  "diabetes_mellitus", "hypertension", "atrial_fibrillation", "prestroke_mrs", "symptom_onset_to_door", "NIHSS_baseline", "occlusion_site", "collateral_score")
p0 <- lapply(plotColumns0, analysisPlot, df = clinicaldata, therapy="IAT")
ggarrange(plotlist = p0, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")


# PSM Matching------------------------------------------------------------------
set.seed(1)
psmMatch <- matchit(
  IAT ~ sex + age + IVT + INR + serum_creatinin + systolic_blood_pressure + 
    diastolic_blood_pressure + previous_stroke + 
    diabetes_mellitus + hypertension + 
    atrial_fibrillation + prestroke_mrs + symptom_onset_to_door +
    NIHSS_baseline + occlusion_site + collateral_score,
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
# 82 patients are unmatched and therefore omitted

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


#Export
install.packages("openxlsx")
library(openxlsx)

write.xlsx(clinicaldata_selection, "Clinical Data Selection.xlsx")
