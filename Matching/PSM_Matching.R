## Propensity Score Matching with MatchIt
## by Jonas Brändli & Nelly Blindenbacher
## June 2023


## Set path to file path-------------------------------------------------------
getwd()
setwd("\\\\fs-home/bline$/Documents/Projekte/9 CT Model vs. Neurologist/2 R Versions/3 Week 27")

## Load data-------------------------------------------------------------------
clinicaldata_cleaned <- read_xlsx("Data_cleaned.xlsx") %>%
  mutate(IAT = as.factor(IAT))

## Source Files-----------------------------------------------------------------
source("analysisPlot.R") # distribution plot

## Dependencies-----------------------------------------------------------------
install.packages(c("dplyr", "ggplot", "MatchIt", "readxl", "ggpubr", "kableExtra"))
library(dplyr)
library(ggplot2)
library(MatchIt)
library(readxl)
library(ggpubr)
library(kableExtra)

# Visualize Data Distribution---------------------------------------------------
plotColumns0 <- c("IAT", "sex", "age", "IVT", "INR", "serum_creatinin", "systolic_blood_pressure", "diastolic_blood_pressure", "previous_stroke", 
                  "diabetes_mellitus", "hypertension", "atrial_fibrillation", "prestroke_mrs", "symptom_onset_to_door", "NIHSS_baseline", "occlusion_site", "collateral_score")
p0 <- lapply(plotColumns0, analysisPlot, df = clinicaldata_cleaned, therapy="IAT")
ggarrange(plotlist = p0, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")


## PSM Matching------------------------------------------------------------------
set.seed(1)
psmMatch <- matchit(
  IAT ~ sex + age + IVT + INR + serum_creatinin + systolic_blood_pressure + 
    diastolic_blood_pressure + previous_stroke + 
    diabetes_mellitus + hypertension + 
    atrial_fibrillation + prestroke_mrs + symptom_onset_to_door +
    NIHSS_baseline + occlusion_site + collateral_score,
  data = clinicaldata_cleaned, 
  distance = "glm", # logistic regression
  link = "logit", # logit link function
  method = "nearest", # matches the two nearest patients in treated and untreated group
  m.order = "random",
  caliper = 0.2) # exludes "bad" matches

# Pre Matching
preMatchDF <- clinicaldata_cleaned
preMatchDF$distance <- psmMatch$distance

# Post Matching
summary(psmMatch)
psmMatchDF <- match.data(psmMatch)
# 73 patients are unmatched and therefore omitted

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

## Reduce data to 20 patients per group----------------------------------------

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
clinicaldata_selection <- psmMatchDF_short1

# Final table of 40 selected patients----------------------------------------------
kable(clinicaldata_selection, format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

#Export excel of 40 selected patients 
install.packages("openxlsx")
library(openxlsx)
write.xlsx(clinicaldata_selection, "Data_matched.xlsx")
