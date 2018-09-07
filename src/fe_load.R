#-----------------------------------------------------------------------
# First episode data
# Load data
# PH, 8/24/2017
#-----------------------------------------------------------------------
source("../src/fe_func.R")
#-----------------------------------------------------------------------
# Declare variables
#-----------------------------------------------------------------------
#atlas     <- "destrieux"
atlas     <- "dk"
#nrois     <- 68
if (atlas == "dk") {
  atlasdf <- dk
  nrois     <- 68
} else {
  atlasdf <- destrieux
  nrois     <- 148
}
#-----------------------------------------------------------------------
# Load CIDAR
fef       <- read.csv("../data/fe_cidar.csv")

# load cidar prior meds
fefmeds <- read.csv("../preproc/clinical/fe_cidar_meds.csv")
fef <- fef %>% left_join(fefmeds)


# load dosage data
fefdose <- read.csv("../preproc/clinical/fe_cidar_dosage.csv") %>%
  dplyr::select(grid, modal.dose.mg)

fef <- fef %>% left_join(fefdose)

# load diagnosis data
fefd <- read.csv("../preproc/clinical/fe_cidar_diagnoses.csv")
fef <- fef %>% left_join(fefd)


# Load OMEGA3
feo       <- read.csv("../data/fe_omega3.csv")

# Exclude 10870 from Omega3 for missing thickness data in one ROI
# when using Destrieux atlas
if (atlas == "destrieux") {
  feo <- feo %>% filter(grid!="10870")
}

# load dosage data
feodose <- read.csv("../preproc/clinical/fe_omega3_dosage.csv") %>%
  dplyr::select(grid, modal.dose.mg)

feo <- feo %>% left_join(feodose)

# load omega3 prior meds
feomeds <- read.csv("../preproc/clinical/fe_omega3_meds.csv")
feo <- feo %>% left_join(feomeds)

# load diagnosis data
feod <- read.csv("../preproc/clinical/fe_omega3_diagnoses.csv")
feo <- feo %>% left_join(feod)

# Load HC data
hcm       <- read.csv("../data/fe_hcm.csv")

visits <- c(0, 1, 2, 3, 4, 6, 8, 10, 12)
slopes1 <- compute_slopes(df=fef %>% filter(visit %in% visits),
                         fm=formula("bprs_td ~ 
                                    log(day+1) + (1 + log(day+1)|grid)"))



slopes2 <- compute_slopes(df=feo %>% filter(visit %in% visits),
                         fm=formula("bprs_td ~ 
                                    log(day+1) + (1 + log(day+1)|grid)"))

# ROI names
rois      <- read.table("../output/tables/fe_freesurfer_roinames.csv",
                        header=TRUE)


slopes1 <- slopes1 %>% bind_cols(fef %>% filter(day==0) %>%
                                 dplyr::select(iq, mccb, dup, age,
                                               sex, race, icv,
                                               med,
                                               bmi,
                                               mccb_reasoning,
                                               mccb_wm,
                                               mccb_speed,
                                               mccb_socialcog,
                                               mccb_attnvig,
                                               mccb_verbalmem,
                                               mccb_visualmem,
                                               bprs_tdbl, hasmri))

slopes2 <- slopes2 %>% bind_cols(feo %>% filter(day==0) %>%
                                 dplyr::select(iq, mccb, dup, age,
                                               sex, race, icv,
                                               med,
                                               bmi,
                                               mccb_reasoning,
                                               mccb_wm,
                                               mccb_speed,
                                               mccb_socialcog,
                                               mccb_attnvig,
                                               mccb_verbalmem,
                                               mccb_visualmem,
                                               bprs_tdbl, hasmri))

# create fe full data frame
# make character types from factors to
# suppress warning messages when joining data frames
slopes1$med <- as.character(slopes1$med)
slopes2$med <- as.character(slopes2$med)
fef$med <- as.character(fef$med)
fef$site <- as.character(fef$site)
fef$study <- as.character(fef$study)
feo$med <- as.character(feo$med)
feo$study <- as.character(feo$study)
fe_full <- rbind(
  slopes1 %>% dplyr::select(grid, slope),
  slopes2 %>% dplyr::select(grid, slope)) %>%
  left_join(rbind(fef, feo), by="grid") 
  # filter out bipolars
  # filter(is.na(bipolar))
    
# create fe data frame with all mri cases
fe <- fe_full %>% filter(hasmri==1)

# create baseline data frame
feb <- fe %>% filter(day==0)


