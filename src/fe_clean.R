# First episode data
# Clean data
# PH, 8/24/2017
#
#library("plyr")
#library("tidyr")
#library("dplyr")
#library("magrittr")
#library("lme4")
#-----------------------------------------------------------------------
source("fe_func.R")
#-----------------------------------------------------------------------
# Load CIDAR
fe <- read.csv("../preproc/fe_original.csv")

# Clean data, remove negative EXAMTYPES
fe <- subset(fe, EXAMTYPE>=0)

# Thought disturbance calculation
# This is according to Sarpal 2016
fe$TotalThinkingDisturb_BPRS <-
  fe$UnusualThoughtContent_BPRS +
  fe$HallucinatoryBehav_BPRS +
  fe$ConceptualDisorg_BPRS +
  fe$Grandiosity_BPRS



# Load medication data
fed <- read.csv("../preproc/fe_cidar_meds_fin.csv")
fed$start <- as.Date(fed$start)
fed$stop <- as.Date(fed$stop)
#fedm <- ddply(fed, c("grid", "drug"), plyr::summarize,
#              daysondrug=max(stop, na.rm=T)-min(start, na.rm=T))
fedm <- fed %>% group_by(grid, drug) %>%
  summarize(daysondrug=max(stop, na.rm=T)-min(start, na.rm=T))

fe <- merge(fe, subset(fedm, drug=="Study Medication1"), all=TRUE)


# Date in fe file
fe$date <- as.Date(fe$EXAMDATE, "%m/%d/%Y")

# Baseline date
bld <- subset(fe, EXAMTYPE==0, select=c(grid, EXAMDATE))
bld$baselinedate <- as.Date(bld$EXAMDATE, "%m/%d/%Y")
fef <- merge(fe, bld[,c(1, 3)], all=TRUE)

# Calculate days
fef$day[fef$EXAMDATE==0] <- 0
fef$day <- as.numeric(as.Date(fef$EXAMDATE, "%m/%d/%Y")-fef$baselinedate)

# Baseline BPRS
feob <- subset(fe, EXAMTYPE==0,
               select=c(grid, TotalThinkingDisturb_BPRS))

colnames(feob) <- c("grid", "baseline")
fef <- merge(fef, feob, all=TRUE)

# Baseline Age, Sex, Race
feob2 <- subset(fe, EXAMTYPE==0,
                select=c(grid, AGE_PDI, SEX_PDI,
                         RACE_PDI, EDUC_PDI, OCCUP_PDI,
                         PATCLASS_PDI))
colnames(feob2) <- c("grid", "age", "sex", "race", "education",
                     "occupation", "patientclass")
fef <- merge(fef, feob2, all=TRUE)

#feg <- read.csv("../preproc/fe_grid.csv")
feg <- read.csv("../preproc/grid.csv")
fed <- read.csv("../preproc/fe_cidar_diagnoses.csv")

# Date in grid file
#feg$date <- as.Date(feg$EXAMDATE, "%d-%b-%y")
#feg <- subset(feg, select=c(grid, SessionNumber, date))

fego <- read.csv("../preproc/fe_moregrid.csv")
#fego$date <- as.Date(fego$EXAMDATE, "%m/%d/%Y")

fes <- read.csv("../preproc/fe_freesurferall.csv")
#fes <- read.csv("../preproc/fe_freesurfer2.csv")

# Merge grids with freesurfer
fegos <- merge(fego, fes)

# Merge Grid, freesurfer with fe
fef <- merge(fef, fegos[, c(1:403)], all=TRUE)

# Merge with diagnoses data
fef <- merge(fef, fed, all=TRUE)

# Merge with fconn dcl data
fefc <- read.csv("../preproc/fe_cidar_fconn.csv")
fef <- merge(fef, fefc, all=TRUE)


fec <- read.csv("../preproc/fe_mccb.csv")
#fes <- read.csv("../preproc/fe_freesurfer.csv")

#fef <- merge(fef, feg, all=TRUE)
#fef$SessionNumber <- fef$SessNo
# Grid with freesurfer
#fegos <- merge(fego, fes)
#fef <- merge(fef, fegos, all=TRUE)
fef <- merge(fef, fec, all=TRUE)


# Get rid of empty data lines
#fef <- subset(fef, !is.na(grid)&!is.na(EXAMTYPE)&!is.na(med)&
#                   !(EXAMTYPE==0&is.na(TotalThinkingDisturb_BPRS)))


# Now merge with final included grids
fegf <- read.csv("../preproc/fe_finalgridsincluded.csv")
colnames(fegf) <- c("majnuid", "grid", "medicationtype")
fegf$included16w <- 1

# Calculate Simpson-Angus Scale score
 fef$sas <- fef$SmpAng01_SimpsonAngus + fef$SmpAng02_SimpsonAngus +
   fef$SmpAng03_SimpsonAngus + fef$SmpAng04_SimpsonAngus +
   fef$SmpAng05_SimpsonAngus + fef$SmpAng06_SimpsonAngus +
   fef$SmpAng07_SimpsonAngus + fef$SmpAng08_SimpsonAngus +
   fef$SmpAng09_SimpsonAngus + fef$SmpAng10_SimpsonAngus +
   fef$SmpAng12_SimpsonAngus
 fef$salivation <- fef$SmpAng10_SimpsonAngus

fef$sas <- fef$SmpAng01 + fef$SmpAng04 + fef$SmpAng07 + fef$SmpAng09 +
  fef$SmpAng10

# Baseline sas
fefbls <- subset(fef, EXAMTYPE==0, select=c(grid, sas))
colnames(fefbls) <- c("grid", "baselinesas")
fef <- merge(fef, fefbls, all=TRUE)

## library("lattice")
## tmp <- fef[order(fef$grid, fef$day),]
## xyplot(Weight_LabValues_191~ EXAMTYPE|grid, data=subset(tmp, EXAMTYPE<=16),
##        ylim=c(100, 300), type="l")

## tmpm <- ddply(tmp, c("grid"), plyr::summarize,
##               Weight_LabValues_191m=mean(Weight_LabValues_191, na.rm=TRUE),
##               Weight_LabValues_191sd=sd(Weight_LabValues_191, na.rm=TRUE),
##               Weight_LabValues_191n=sum(!is.na(Weight_LabValues_191)),
##               Weight_LabValues_191se=Weight_LabValues_191sd/sqrt(Weight_LabValues_191n))

## hist(tmpm$Weight_LabValues_191m)

## tmpm <- ddply(tmp, c("grid"), plyr::summarize,
##               sasm=mean(sas, na.rm=TRUE),
##               sassd=sd(sas, na.rm=TRUE),
##               sasn=sum(!is.na(sas)),
##               sasse=sassd/sqrt(sasn))
## pdf("../output/figures/fef_cidar_sas.pdf")
## par(mar=c(5, 5, 5, 5), oma=c(1, 1, 1, 1))
## hist(log(tmpm$sasm[tmpm$sasm>=0]), col="gray",
##      main="", xlab="Log(SAS)", cex.lab=2, font=2)
## dev.off()


## tmpm <- ddply(tmp, c("grid"), plyr::summarize,
##               salivationm=mean(salivation, na.rm=TRUE),
##               salivationsd=sd(salivation, na.rm=TRUE),
##               salivationn=sum(!is.na(salivation)),
##               salivationse=salivationsd/sqrt(salivationn))
## pdf("../output/figures/fef_cidar_salivation.pdf")
## par(mar=c(5, 5, 5, 5), oma=c(1, 1, 1, 1))
## hist(tmpm$salivationm[tmpm$salivationm>=0], col="gray",
##      main="", xlab="Log(SALIVATION)", cex.lab=2, font=2)
## dev.off()


# Remove ambigous or redundant variables 
fef$SEX_PDI <- NULL
fef$AGE_PDI <- NULL
fef$EDUC_PDI <- NULL
fef$RACE_PDI <- NULL
fef$OCCUP_PDI <- NULL
fef$Sex <- NULL
#fef$Age.T0 <- NULL
#fef$Age.T1 <- NULL
fef <- merge(fef, fegf[,c(2,4)], all=TRUE)


# Thought disturbance calculation
# This is according to Sarpal 2015
#fef$TotalThinkingDisturb_BPRS <-
#  fef$UnusualThoughtContent_BPRS +
#  fef$HallucinatoryBehav_BPRS +
#  fef$ConceptualDisorg_BPRS 

fef$TotalThinkingDisturb_BPRS[fef$TotalThinkingDisturb_BPRS < 0] <- NA

write.csv(fef, "../data/fe_cidar.csv", na="", row.names=FALSE)

# Build Omega 3
feo <- read.csv("../preproc/fe_omega3_original.csv")
feo <- feo[order(feo$grid, feo$EXAMTYPE), ]

# Thought disturbance calculation
# This is according to Sarpal 2016
feo$TotalThinkingDisturb_BPRS <-
  feo$UnusualThoughtContent +
  feo$HallucinatoryBehav +
  feo$ConceptualDisorg +
  feo$Grandiosity


# Clean data, remove negative EXAMTYPES and MCCB values
feo <- subset(feo, EXAMTYPE>=0)

# Baseline date
bld <- subset(feo, EXAMTYPE==0, select=c(grid, EXAMDATE))
bld$baselinedate <- as.Date(bld$EXAMDATE, "%m/%d/%Y")
feo <- merge(feo, bld[,c(1, 3)], all=TRUE)

# Calculate days
feo$day <- as.numeric(as.Date(feo$EXAMDATE, "%m/%d/%Y")-feo$baselinedate)

# Baselines
feobl <- read.csv("../preproc/fe_omega3_baseline.csv")
feo <- merge(feo, feobl, all=TRUE)
feob <- subset(feo, EXAMTYPE==0,
               select=c(grid, TotalThinkingDisturb_BPRS))

colnames(feob) <- c("grid", "baseline")
feo <- merge(feo, feob, all=TRUE)


feosc <- read.csv("../preproc/fe_omega3_scans.csv")

# Restric scan data to baseline and grop variable EXAMTYPE
feosc <- feosc %>% filter(EXAMTYPE==0) %>% dplyr::select(-EXAMTYPE)

fec <- read.csv("../preproc/fe_omega3_cognition.csv")
fec <- subset(fec, OverallCompositeScore>=0)

# Restric cognition data to baseline and grop variable EXAMTYPE
fec <- fec %>% filter(EXAMTYPE==0) %>% dplyr::select(-EXAMTYPE)

feov <- read.csv("../preproc/fe_omega3_freesurfer.csv")
feogf <- read.csv("../preproc/fe_omega3_finalgridsincluded.csv")
feogf$included16w <- 1

feos <- merge(feov, feosc, all=TRUE)
feo <- merge(feo, feos, all=TRUE)
feo <- merge(feo, fec, all=TRUE)
feo <- merge(feo, feogf, all=TRUE)

# Calculate SAS
feo$sas <- feo$SmpAng01 + feo$SmpAng02 + feo$SmpAng03 + feo$SmpAng04 +
  feo$SmpAng05 + feo$SmpAng06 + feo$SmpAng07 + feo$SmpAng08 +
  feo$SmpAng09 + feo$SmpAng10

# Baseline SAS
feobls <- subset(feo, EXAMTYPE==0, select=c(grid, sas))
colnames(feobls) <- c("grid", "baselinesas")
feo <- merge(feo, feobls, all=TRUE)

# Thought disturbance calculation
# This is according to Sarpal 2016
feo$TotalThinkingDisturb_BPRS <-
  feo$UnusualThoughtContent +
  feo$HallucinatoryBehav +
  feo$ConceptualDisorg +
  feo$Grandiosity

# Thought disturbance calculation
# This is according to Sarpal 2015
#feo$TotalThinkingDisturb_BPRS <-
#  feo$UnusualThoughtContent +
#  feo$HallucinatoryBehav +
#  feo$ConceptualDisorg 


feo$TotalThinkingDisturb_BPRS[feo$TotalThinkingDisturb_BPRS < 0] <- NA


# Build chronic ZHII
#zhi <- read.csv("../preproc/fe_chronic.csv")
zhi <- read.csv("../preproc/fe_zhii_freesurfer_baseline.csv")
zhic <- read.csv("../preproc/fe_zhii_clinical.csv")
zhig <- read.csv("../preproc/fe_zhii_gender.csv")
zhi <- merge(zhi, zhic, all=TRUE)
zhi <- merge(zhi, zhig, all=TRUE)
zhii <- subset(zhi, !is.na((TotalThinkingDisturb_BPRS)&Time==0)&
                   !is.na(grid))

zhii$Date.T0 <- as.Date(zhii$start.of.inpatient.stay, "%m/%d/%Y")
zhii$Date.T1 <- as.Date(zhii$end.of.stay, "%m/%d/%Y")
zhii$LengthOfStay <- as.numeric(zhii$Date.T1 - zhii$Date.T0)
#zhii$day.T0 <- 0
#zhii$day.T1 <- zhii$LengthOfStay 
zhii$day[zhii$Time==0] <- 0
zhii$day[zhii$Time==1] <- zhii$LengthOfStay[zhii$Time==1]


# Write final files to disk
write.csv(fef, "../preproc/fe_cidarold.csv", na="", row.names=FALSE)
write.csv(feo, "../preproc/fe_omega3old.csv", na="", row.names=FALSE)
write.csv(zhii, "../preproc/fe_zhiiold.csv", na="", row.names=FALSE)
#-----------------------------------------------------------------------
# Additional cleaning
#-----------------------------------------------------------------------
# Harmonize the data sets with respect to variable naming
# In addition, all 68 DK-Freesurfer ROIs
rois <- read.csv("../output/tables/fe_freesurfer_roinames.csv")
rois$name <- as.character(rois$name)
#rois <- rois %>%
#    bind_rows(rois) %>%
#    mutate(name=c(paste("l", name[1:34], sep="_"),
#                  paste("r", name[35:68], sep="_")))
#-----------------------------------------------------------------------
# CIDAR 
#-----------------------------------------------------------------------
# Also load education
educ <- read.csv("../preproc/fe_demographics_education.csv")

fefnew <- fef %>% dplyr::select(
                             grid,
                             SessionNumber,
                             Age.T0,
                             sex,
                             race,
                             education,
                             StudySite,
                             med,
                             #starts_with("SANS", ignore.case=FALSE),
                             Total_1_18_BPRS,
                             TotalThinkingDisturb_BPRS,
                             baseline,
                             sas,
                             MCCB.GenCog.Tscore.T0,
                             MCCB.WM.T0,
                             MCCB.Speed.T0,
                             MCCB.Reasoning.T0,
                             MCCB.VerbalMemory.T0,
                             MCCB.VisualMemory.T0,
                             MCCB.SocialCog.T0,
                             MCCB.AttnVig.T0,
                             BMI_LabValues_191,
                             Weight_LabValues_191,
                             PATCLASS_PDI,
                             PSYMP1ST_PDI,
                             ALCOHOL_PDI,
                             DRUGABUS_PDI,
                             WRAT3ReadingSS.T0,
                             WKPSYMP_PDI,
                             included16w,
                             EXAMTYPE,
                             IntraCranialVol,
                             matches(paste(rois$name, collapse="|")),
                             day) %>%
    dplyr::select(-matches("white"), -matches("SD"),
                  -contains("grayVolume")) %>%
    rename(grid=grid,
           age=Age.T0,
           site=StudySite,
           med=med,
           sex=sex,
           race=race,
           bprs=Total_1_18_BPRS,
           bprs_td=TotalThinkingDisturb_BPRS,
           bprs_tdbl=baseline,
           mccb=MCCB.GenCog.Tscore.T0,
           mccb_wm=MCCB.WM.T0,
           mccb_speed=MCCB.Speed.T0,
           mccb_reasoning=MCCB.Reasoning.T0,
           mccb_verbalmem=MCCB.VerbalMemory.T0,
           mccb_visualmem=MCCB.VisualMemory.T0,
           mccb_socialcog=MCCB.SocialCog.T0,
           mccb_attnvig=MCCB.AttnVig.T0,
           bmi=BMI_LabValues_191,
           weight=Weight_LabValues_191,
           class=PATCLASS_PDI,
           onset=PSYMP1ST_PDI,
           alcohol=ALCOHOL_PDI,
           drugs=DRUGABUS_PDI,
           iq=WRAT3ReadingSS.T0,
           dup=WKPSYMP_PDI,
           included16w=included16w,
           visit=EXAMTYPE,
           icv=IntraCranialVol,
           day=day) %>%
    mutate(site=as.character(site)) %>%
    dplyr::select(-education) %>%  # drop initial education value
    left_join(educ)

fefrois <- fefnew %>%
    dplyr::select(matches(paste(rois$name, collapse="|"))) %>% 
    magrittr::set_names(c(paste("l", rois$name, sep="_"),
                          paste("r", rois$name, sep="_"))) 
fefnew <- fefnew %>% dplyr::select(-contains("gray")) %>%
    bind_cols(fefrois) %>%
    mutate(hasmri=ifelse(!is.na(SessionNumber), 1, 0),
           study="CIDAR") %>%
    filter(included16w==1)

fefnew <- fefnew %>%
    mutate(hasmri=ifelse(!is.na(SessionNumber), 1, 0)) %>%
    dplyr::select(grid, hasmri, day) %>%
    filter(day==0) %>%
    dplyr::select(-day) %>%
    right_join(fefnew) %>%
    dplyr::rename(sessnum=SessionNumber)

# Replace negative values by NA
fefnew <- fefnew %>%  mutate_each(funs(replace(., .<0, NA)))

# Correct Risperdone to Risperidone
fefnew$med <- as.character(fefnew$med)
fefnew$med[fefnew$med=="Risperdone"] <- "Risperidone"
fefnew$med <- factor(fefnew$med)

# Write to disk
write.csv(fefnew, "../data/fe_cidar.csv", row.names=FALSE, na="")
#-----------------------------------------------------------------------
# Omega 3
#-----------------------------------------------------------------------
feo$StudySite <- "ZHH"
feo$med <- "SGA"
feonew <- feo %>% dplyr::select("grid", "SomaticConcern", "Anxiety",
                                   "EmotionalWithdrawal",
                                   "ConceptualDisorg", "GuiltFeelings",
                                   "Tension", "MannerismsPosturing",
                                   "Grandiosity", "DepressiveMood",
                                   "Hostility", "Suspiciousness",
                                   "HallucinatoryBehav", "MotorRetard",
                                   "Uncooperativeness",
                                   "UnusualThoughtContent",
                                   "BluntedAffect", "Excitement",
                                   "Disorientation") %>%
    mutate_all(funs(replace(., .<0, NA))) %>%
    mutate(Total_1_18_BPRS=rowSums(.[2:19], na.rm=TRUE)) %>%
    dplyr::select(Total_1_18_BPRS)

#feonew <- feo %>% bind_cols(feonew) %>%

feonew <- cbind(feo, feonew)  
    

feonew <- feonew %>% dplyr::select(
                             grid,
                             SessionNumber,
                             AGE,
                             SEX,
                             RACE,
                             EDUC,
                             StudySite,
                             #starts_with("SANS", ignore.case=FALSE),
                             Total_1_18_BPRS,
                             TotalThinkingDisturb_BPRS,
                             baseline,
                             sas,
                             OverallCompositeScore,
                             WorkingMemory,
                             ReasoningProblemSolv,
                             SpeedProcessing,
                             VisualLearning,
                             VerbalLearning,
                             AttentionVigilance,
                             SocialCognition,
                             BMI,
                             Weight_BodyAsses,
                             PATCLASS,
                             PSYMP1ST,
                             ALCOHOL,
                             DRUGABUS,
                             WRAT3SS,
                             WKPSYMP,
                             included16w,
                             EXAMTYPE,
                             IntraCranialVol,
                             matches(paste(rois$name, collapse="|")),
                             day) %>%
    dplyr::select(-matches("white"), -matches("SD"),
                  -contains("grayVolume")) %>%
    rename(grid=grid,
           age=AGE,
           sex=SEX,
           race=RACE,
           education=EDUC,
           site=StudySite,
           bprs=Total_1_18_BPRS,
           bprs_td=TotalThinkingDisturb_BPRS,
           bprs_tdbl=baseline,
           mccb=OverallCompositeScore,
           mccb_wm=WorkingMemory,
           mccb_reasoning=ReasoningProblemSolv,
           mccb_speed=SpeedProcessing,
           mccb_visualmem=VisualLearning,
           mccb_verbalmem=VerbalLearning,
           mccb_attnvig=AttentionVigilance,
           mccb_socialcog=SocialCognition,
           bmi=BMI,
           weight=Weight_BodyAsses,
           class=PATCLASS,
           onset=PSYMP1ST,
           alcohol=ALCOHOL,
           drugs=DRUGABUS,
           iq=WRAT3SS,
           dup=WKPSYMP,
           included16w=included16w,
           visit=EXAMTYPE,
           icv=IntraCranialVol,
           day=day) %>%
    dplyr::select(-education) %>%  # drop initial education value
    left_join(educ)

feorois <- feonew %>%
    dplyr::select(matches(paste(rois$name, collapse="|"))) %>% 
    magrittr::set_names(c(paste("l", rois$name, sep="_"),
                          paste("r", rois$name, sep="_"))) 
feonew <- feonew %>% dplyr::select(-contains("gray")) %>%
    bind_cols(feorois) %>%
    mutate(hasmri=ifelse(!is.na(SessionNumber), 1, 0),
           study="OMEGA3") %>%
    filter(included16w==1)

feon <- feonew %>%
    rename(sessnum=SessionNumber) %>%
    mutate(hasmri=ifelse((is.na(sessnum)|
                         is.na(icv)), 0, 1)) %>%
    dplyr::select(grid, hasmri, day) %>%
    filter(day==0) %>%
    dplyr::select(-day) %>%
    right_join(feonew) 

feonew <- feon %>%
    rename(sessnum=SessionNumber) %>%
    mutate(hasmri=ifelse((is.na(sessnum)|
                          is.na(icv)), 0, 1),
           med="Risperidone")

# Replace negative values by NA
feonew <- feonew %>%  mutate_each(funs(replace(., .<0, NA)))


# Write to disk
write.csv(feonew, "../data/fe_omega3.csv", row.names=FALSE, na="")
#-----------------------------------------------------------------------
# HC data from Miklos 
#-----------------------------------------------------------------------
hcm <- read.csv("../preproc/fe_freesurfer_hcfrommiklos.csv")
hcm_cog <- read.csv("../preproc/fe_hcfrommiklos_cognition.csv")
hcm_dem <- read.csv("../preproc/fe_hcm_demographics.csv")
summary(hcm)

hcm_cog <- hcm_cog %>% dplyr::select(GRID, SessNo, age, sex,
                                     NP_race, EXAMTYPE,
                                     EXAMDATE, WRAT3SS) %>%
    rename(grid=GRID,
           sessno=SessNo,
           race=NP_race,
           visit=EXAMTYPE,
           date=EXAMDATE,
           iq=WRAT3SS)

hcm_cogm <- hcm_cog %>% dplyr::select(grid, date, sessno, race,
                                 visit, iq) %>%
    mutate_all(funs(replace(., .<0, NA))) %>%
    group_by(grid) %>%
    dplyr::summarize(
               iq=mean(iq, na.rm=TRUE)) %>%
    left_join(hcm_cog %>% dplyr::select(grid, race) %>%
              distinct() %>% na.omit())

hcmrois <- hcm %>%
    dplyr::select(matches(paste(rois$name, collapse="|"))) %>% 
    magrittr::set_names(c(paste("l", rois$name, sep="_"),
                          paste("r", rois$name, sep="_"))) 

hcmnew <- hcm %>%
    dplyr::rename(icv=IntraCranialVol) %>%
    dplyr::select(-matches("Left|Right|Brain|Ventricle|
                            Gray|WM|CC|Cortical|Cortex|
                            CSF|Vol|Edn|EXAM|Session|
                            between|Optic")) %>%
    dplyr::select(grid, sessnum, age, sex, icv) %>%
    bind_cols(hcmrois) %>%
    mutate(day=0) %>%
    left_join(hcm_cogm) %>%
    left_join(hcm_dem %>% dplyr::select(grid, education)) %>%
    #left_join(hcm_cog %>% dplyr::select(grid, race))
    filter(age<50, age>=15)

summary(hcmnew)
write.csv(hcmnew, "../data/fe_hcm.csv", row.names=FALSE, na="")
