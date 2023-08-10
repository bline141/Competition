## Selection & Propensity Score Matching with MatchIt
## by J.B. & N.B.
## June 2023

# Set path to file path

# Dependencies
library(readxl)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(MatchIt)
library(ggpubr)
library(openxlsx)
library(VIM)
library(openxlsx)

# Source Files
source("analysisPlot.R") # distribution plot
source("pairPlot.R")
source("get_pval.R")
source("get_sdiff.R")

# Load data
clinicaldata_cleaned <- read_xlsx("Data_cleaned.xlsx") %>%
  mutate(IAT = as.factor(IAT))

## Selection for competition

# Select relevant baseline variables
clinicaldata_variables <- clinicaldata_cleaned %>%
  select(c("IAT", "sex", "age", "IVT", "INR", "serum_creatinin", "systolic_blood_pressure", "diastolic_blood_pressure", "previous_stroke", 
           "diabetes_mellitus", "hypertension", "atrial_fibrillation", "prestroke_mrs", "symptom_onset_to_door", "NIHSS_baseline", "occlusion_site", "collateral_score", "followid"))

# Visualize Missing data
anyNA(clinicaldata_variables)
clinicaldata_missing <- clinicaldata_variables[sapply(clinicaldata_variables, anyNA)]
aggr(clinicaldata_missing, 
     plot = TRUE, numbers = TRUE, prop = FALSE,
     labels = names(clinicaldata_missing), cex.axis = .9, oma = c(10,5,5,3))


# Exclude patients with missing BL data
missing_values <- colSums(is.na(clinicaldata_variables))
print(missing_values) # data is not systematically missing
clinicaldata_selection <- clinicaldata_variables %>% 
  filter(!is.na(INR), !is.na(serum_creatinin), !is.na(symptom_onset_to_door), !is.na(collateral_score))

# Imputation with missForest is not necessary, because patients with missing data are excluded


## Matching

# Visualize Data Distribution
plotColumns0 <- c("IAT", "sex", "age", "IVT", "systolic_blood_pressure", "diastolic_blood_pressure", "previous_stroke", 
                  "diabetes_mellitus", "hypertension", "atrial_fibrillation", "prestroke_mrs", "NIHSS_baseline", "occlusion_site")
p0 <- lapply(plotColumns0, analysisPlot, df = clinicaldata_selection, therapy="IAT")
ggarrange(plotlist = p0, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")


## PSM Matching
set.seed(1)
psmMatch <- matchit(
  IAT ~ sex + age + IVT + systolic_blood_pressure + 
    diastolic_blood_pressure + previous_stroke + 
    diabetes_mellitus + hypertension + 
    atrial_fibrillation + prestroke_mrs +
    NIHSS_baseline + occlusion_site,
  data = clinicaldata_selection, 
  distance = "glm", # logistic regression
  link = "logit", # logit link function
  method = "nearest", # matches the two nearest patients in treated and untreated group
  m.order = "random",
  caliper = 0.2) # exludes "bad" matches

# Pre Matching
preMatchDF <- clinicaldata_selection
preMatchDF$distance <- psmMatch$distance

# Post Matching
summary(psmMatch)
psmMatchDF <- match.data(psmMatch)
# 58 patients are unmatched and therefore omitted

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


## calculate test:
#         Paired    No                  Yes
# Values           
#    Continuous      Wilcoxon rank-sum   Wilcoxon signed-rank
#    Categorical     Chi-squared         McNemar

all_pvals <- rbind(
  sapply(plotColumns, 
         function(var, dat){get_pval(variable = dat[[var]], group = dat$IAT)},
         dat = preMatchDF),
  sapply(plotColumns, 
         function(var, dat){get_pval(variable = dat[[var]], group = dat$IAT, pair_id = dat$subclass)},
         dat = psmMatchDF),
  sapply(plotColumns, 
         function(var, dat){get_pval(variable = dat[[var]], group = dat$IAT, pair_id = dat$subclass)},
         dat = psmMatchDF_short1),
  sapply(plotColumns, 
         function(var, dat){get_pval(variable = dat[[var]], group = dat$IAT, pair_id = dat$subclass)},
         dat = psmMatchDF_short2)
)
all_pvals <- t(all_pvals)
colnames(all_pvals) <- c("PreMatch", "Matched", "Reduced_1", "Reduced_2")
round(all_pvals, 3)


# Plot standardized Differences
# numerical: standardized difference between groups
# categorical: standardized difference of proportions between groups
sdiff_short1 <- sapply(plotColumns, 
               function(var, dat){get_sdiff(variable = dat[[var]], group = dat$IAT)},
               dat = psmMatchDF_short1)
sdiff_short2<- sapply(plotColumns, 
                       function(var, dat){get_sdiff(variable = dat[[var]], group = dat$IAT)},
                       dat = psmMatchDF_short2)

par(mar = c(12,4,2,2))
plot(x = sdiff_short1, type = "b", col = 1, ylim = c(-0.5, 0.6), 
     xaxt = "n", xlab = "", ylab = "standardized difference")
axis(1, at=1:length(plotColumns), labels=plotColumns, las = 2)
lines(x = sdiff_short2, type = "b", col = 2)
abline(h = 0, col = 3, lwd = 2)

# Plot Differences after Matching 
pair_diff_psm <- pairDiffCalc(psmMatchDF)
pair_diff_short1 <- pairDiffCalc(psmMatchDF_short1)
pair_diff_short2 <- pairDiffCalc(psmMatchDF_short2)

# 1 col all matched differences, 
# 2 col differences after reducing to 20 pairs, with smallest diff
# 3 col differences after reducing to 20 pairs, random
p20 <- lapply(plotColumns, pairPlot, df = pair_diff_psm)
p21 <- lapply(plotColumns, pairPlot, df = pair_diff_short1)
p22 <- lapply(plotColumns, pairPlot, df = pair_diff_short2)
newOrder <- as.vector(matrix(1:(length(plotColumns)*3), nrow = 3, byrow = T))
ggarrange(plotlist = c(p20,p21,p22)[newOrder], nrow = 4, ncol = 3,  
          common.legend = TRUE, legend = "bottom")

## Adjustments and Statistical Tests
# Save psmMatchDF_short1 as a new dataframe to preserve the original data
psmMatchDF_short1_adj <- psmMatchDF_short1

# Remove pair 45, 65 and 73 
psmMatchDF_short1_adj <- filter(psmMatchDF_short1_adj, !(subclass %in% c(45, 65, 73)))

# Create data frames for subclass 54, 63 and 66
subclass_54_obs <- filter(psmMatchDF, subclass == 54)
subclass_63_obs <- filter(psmMatchDF, subclass == 63)
subclass_66_obs <- filter(psmMatchDF, subclass == 66)

# Add the rows for subclass 54, 63 and 66 to psmMatchDF_short1_adj
psmMatchDF_short1_adj <- bind_rows(psmMatchDF_short1_adj, subclass_54_obs, subclass_63_obs, subclass_66_obs)


# Plot standardized Differences
# numerical: standardized difference between groups
# categorical: standardized difference of proportions between groups
sdiff_short1 <- sapply(plotColumns, 
                       function(var, dat){get_sdiff(variable = dat[[var]], group = dat$IAT)},
                       dat = psmMatchDF_short1)
sdiff_short2<- sapply(plotColumns, 
                      function(var, dat){get_sdiff(variable = dat[[var]], group = dat$IAT)},
                      dat = psmMatchDF_short1_adj)

par(mar = c(12,4,2,2))
plot(x = sdiff_short1, type = "b", col = 1, ylim = c(-0.5, 0.6), 
     xaxt = "n", xlab = "", ylab = "standardized difference")
axis(1, at=1:length(plotColumns), labels=plotColumns, las = 2)
lines(x = sdiff_short2, type = "b", col = 2)
abline(h = 0, col = 3, lwd = 2)

# Plot Differences after Adjustments 
pair_diff_adj <- pairDiffCalc(psmMatchDF_short1_adj)

# 1 col differences after reducing to 20 pairs, with smallest diff
# 2 adjusted reducing
p23 <- lapply(plotColumns, pairPlot, df = pair_diff_adj)
newOrder <- as.vector(matrix(1:(length(plotColumns)*2), nrow = 2, byrow = T))
ggarrange(plotlist = c(p21,p23)[newOrder], nrow = 4, ncol = 2,  
          common.legend = TRUE, legend = "bottom")

sel <- which(plotColumns %in% c("prestroke_mrs", "collateral_score"))
ggarrange(plotlist = c(p21[sel],p23[sel])[c(1,3,2,4)], nrow = 2, ncol = 2,  
          common.legend = TRUE, legend = "bottom")

## Export Patient Data for Competition

# Display the 20 pairs with the smallest difference in propensity score
clinicaldata_competition <- psmMatchDF_short1_adj

# Final table of 40 selected patients
kable(clinicaldata_competition, format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

#Export excel of 40 selected patients 
write.xlsx(clinicaldata_competition, "Data_matched.xlsx")

