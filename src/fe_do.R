#-----------------------------------------------------------------------
# Structural similarity network analysis
#
# Philipp Homan, 13/3/18
# <phoman@northwell.edu>
#-----------------------------------------------------------------------
source("../src/fe_load.R")
# Create file marker
system("touch ../output/R/fe_do.Rout")
#-----------------------------------------------------------------------
# Declare variables
#-----------------------------------------------------------------------
groups    <- list("cidar", "omega3")
# Random seed for the following analyses
set.seed(20180214)
niter <- 10^3
npermut <- 10^3

# params dictionary for graph metrics
graph_params_dict <- c(
              `assortativity` = "Assortativity",
              `diameter` = "Diameter",
              `E.global` = "Global efficiency",
              `E.local` = "Local efficiency",
              `Lp` = "Path length",
              `num.hubs` = "Number of hubs",
              `transitivity` = "Clustering coef.",
              `vulnerability` = "Vulnerability",
              `smallworld` = "Small worldness",
              `richclub` = "Rich club coef."
             )

# sparsity threshold
# This gave rounding problems (!)
# sparsity  <- seq(0.1, 0.7, 0.05)
sparsity <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
              0.6, 0.65, 0.7)

# Lists of data sets and slopes
all       <- list(fef, feo)
allslopes <- list(slopes1, slopes2)
#-----------------------------------------------------------------------
# Demonstrate similarity of schizophrenia cohorts
#-----------------------------------------------------------------------
exweeks <- c(0, 1, 2, 3, 4, 6, 8, 10, 12)
fef1 <- fef %>% dplyr::select(grid, bprs_td,
                              visit, included16w)  %>%
    rename(y=bprs_td) %>%
  mutate(study="CIDAR")

feo1 <- feo %>% dplyr::select(grid, bprs_td,
                              visit, included16w)  %>%
    rename(y=bprs_td) %>%
  mutate(study="OMEGA3")

p <- (plot_average(rbind(fef1, feo1) %>% filter(visit %in% exweeks)))
ggsave("../output/figures/fe_averagetx.pdf", plot=p, width=16, height=9, scale=0.5)
#-----------------------------------------------------------------------
# Extract grids for each group
#-----------------------------------------------------------------------
grid <- lapply(all,
               function(x) {
                 filter(x,
                        hasmri==1 &
                        day==0) %>%
                   dplyr::select(grid)
               })
#-----------------------------------------------------------------------
# Extract slope data
#-----------------------------------------------------------------------
slopelist <- lapply(allslopes,
               function(x) {
                 x <- filter(x, hasmri==1) %>%
                   dplyr::select(grid, slope, icv, age, sex,
                                 mccb, iq, bprs_tdbl, dup, race)
               })
#-----------------------------------------------------------------------
# Estimate shrinkage for each group
#-----------------------------------------------------------------------
dflist <- list(feo %>% filter(included16w==1),
               fef %>% filter(included16w==1))

# Anonymize grids to get lower grid numbers for illustration
dflist <- lapply(dflist, function(x) anon_grid(x))
dflist <- lapply(dflist,
              function(x) {
                  dplyr::rename(x,
                                y=bprs_td,
                                x=day,
                                visit=visit) %>%
                    mutate(grid=as.character(grid_anon),
                           x=log(x+1)) %>%
                    filter(visit <= 12,
                           included16w==1)
                })

glist <- lapply(dflist, function(x) plot_shrinkage(x, fm=NULL))

p1 <- plot_grid(plotlist=glist[[2]][2:3],
                scale=0.95, labels=c("A", "B"), label_size=40)

p2 <- plot_grid(plotlist=glist[[1]][2:3],
                scale=0.95, labels=c("C", "D"), label_size=40)
p3 <- plot_grid(p1, p2, labels=NULL, nrow=2)

p4 <- plot_grid(glist[[2]][[4]], glist[[1]][[4]], labels="AUTO",
                label_size=40, scale=0.9, nrow=2)

ggsave(filename="../output/figures/fe_freesurfer_shrinkage.pdf",
       plot=p3, width=30, height=30)

ggsave(filename="../output/figures/fe_freesurfer_shrinkagered.pdf",
       plot=p4, width=20, height=25)
#-----------------------------------------------------------------------
# Estimate graphs for each group
#-----------------------------------------------------------------------
#atlasstr = "_aparc.a2009s.annot"
#nrois = 148
#i <- 0
#graphlist <- lapply(groups,
#                 function(x) {
#                   i <<- i + 1 
#                   make_graphs(x, sparsity, nsp,
#                               grid[[i]]$grid, atlasstr, nrois)
#                 })
#saveRDS(graphlist, "../data/fe_graphlist_destrieux.rds")
fn <- paste("../data/fe_graphlist_", atlas, ".rds", sep="")
graphlist <- readRDS(fn)
rm(fn)

# Use only the brainGraph output, the first element for each group
grlist <- list(graphlist[[1]][[1]], graphlist[[2]][[1]])
#-----------------------------------------------------------------------
# Extract degree for each hub in each subject and each sparsity
#-----------------------------------------------------------------------
degreelist <- lapply(grlist,
                     function(x) {
                       lapply(c(1:length(sparsity)),
                              function (y) {
                                calc_degree(x, ncol=nrois, spi=y)
                              })
                     })
#saveRDS(degreelist, "../data/fe_degreelist_destrieux.rds")
#degreelist <- readRDS("../data/fe_degreelist_destrieux.rds")
#-----------------------------------------------------------------------
# Extract graph params for each subject and each sparsity
#-----------------------------------------------------------------------
#paramslist <- lapply(grlist,
#                     function(x) {
#                       lapply(c(1:length(sparsity)),
#                              function (y) {
#                                get_graph_params(x, y)
#                              })
#                     })
#saveRDS(paramslist, "../data/fe_paramslist_destrieux.rds")
fn <- paste("../data/fe_paramslist_", atlas, ".rds", sep="")
paramslist <- readRDS(fn)

# Reorganize grids
allgrids <- do.call("rbind", grid) %>%
  mutate(study=rep(unlist(groups),
                   times=c(length(slopes1$slope[slopes1$hasmri==1]),
                           length(slopes2$slope[slopes2$hasmri==1]))))

allslopes <- do.call("rbind", slopelist) %>%
  mutate(study=rep(unlist(groups), 
                   times=c(length(slopes1$slope[slopes1$hasmri==1]),
                           length(slopes2$slope[slopes2$hasmri==1]))))
#-----------------------------------------------------------------------
# Create data frame from the list of degrees
#-----------------------------------------------------------------------
Ddf <- do.call("rbind", do.call("rbind", degreelist)) %>%
  mutate(study=rep(rep(unlist(groups),
                       times=c(length(slopes1$slope[slopes1$hasmri==1]),
                               length(slopes2$slope[slopes2$hasmri==1]))),
                   length(sparsity)),
         grid=rep(allgrids$grid, length(sparsity)),
         density=rep(sparsity, each=nrow(allgrids)))

# Merge with slope
dms <- merge(Ddf, allslopes %>%
                  dplyr::select(-study), by="grid")
#-----------------------------------------------------------------------
# Residualize data.
#
# Start with slope. Although not needed, apply it to all levels of
# density so that the merging is identical to the next residualize
# operation
#-----------------------------------------------------------------------
fm            <- "slope ~ age * sex + icv + study + bprs_tdbl"
dmslist       <- lapply(sparsity,
                        function(x) {
                          residualize_df(subset(dms, density==x),
                                         fm, (nrois+4))
                        })
dms$slope     <- do.call("rbind", dmslist)$slope

fm            <- "y ~ age * sex + icv + study + bprs_tdbl"
#dms           <- residualize_df(dms, fm, c(2:69))
dmslist       <- lapply(sparsity,
                        function(x) {
                          residualize_df(subset(dms, density==x),
                                         fm, c(2:(nrois-1)))
                        })
dms           <- do.call("rbind", dmslist)
  
# Remove variables that are not needed any longer
dms           <- dplyr::select(dms,
                     -grid, -study, -age, -sex, -iq, -mccb, 
                     -icv, -bprs_tdbl, -dup, -race) 

colnames(dms) <- c(atlasdf$name, "density", "slope")

#-----------------------------------------------------------------------
# PLS regression for all levels of sparsity
#-----------------------------------------------------------------------
plsfitlist <- lapply(sparsity,
                     function(x) {
                       plsr(slope ~ ., data=subset(dms, density==x) %>%
                            dplyr::select(-density))
                     })
                     
predlist <- lapply(c(1:length(sparsity)),
                   function(x) {
                     predict(plsfitlist[[x]], ncomp=2,
                             newdata=subset(dms, density==sparsity[x]))
                   })
#-----------------------------------------------------------------------
# Summarize the the results 
#-----------------------------------------------------------------------
plssum <- do.call("rbind", lapply(c(1:length(sparsity)),
                           function(x) {
                             defaultSummary(
                               data.frame(
                                 obs=dms$slope[dms$density==sparsity[x]],
                                 pred=as.numeric(predlist[[x]])))
                           })) %>%
  tibble::as_data_frame() %>%
  mutate(density=sparsity)
#-----------------------------------------------------------------------
# Cross-validation k-fold
#-----------------------------------------------------------------------
cvkfold_k <- 5
cvlist <- cross_validate(subset(dms, sparsity==0.1),
                         k=cvkfold_k, niter=niter)
cvrsq <- (cvlist[[2]])
print(paste("Cross-validated R2: ", round(mean(cvrsq), 2), " (",
round(sd(cvrsq), 2), ")", sep=""))
fn <- paste("../output/tables/fe_msn_crossvalidation_", atlas,
            ".csv", sep="")
write.csv(data.frame(mean=mean(cvrsq), sd=sd(cvrsq), k=cvkfold_k,
                       density=0.1, niter=niter),
          fn,
          row.names=FALSE)
rm(fn)
#-----------------------------------------------------------------------
# Transform paramslist to data frame
#-----------------------------------------------------------------------
Pdf       <- do.call("rbind", do.call("rbind", paramslist)) %>%
  mutate(study=rep(rep(unlist(groups), 
                       times=c(length(slopes1$slope[slopes1$hasmri==1]),
                               length(slopes2$slope[slopes2$hasmri==1]))),
                   length(sparsity)),
         grid=rep(allgrids$grid, length(sparsity)))

pms       <- merge(Pdf, allslopes %>% dplyr::select(-study), by="grid")
pmsl      <- pms %>% gather(key=metric, value=estimate,
                       Cp, Lp, E.global, E.local, diameter,
                       transitivity, num.hubs, vulnerability, 
                       assortativity, smallworld, richclub) %>%
               arrange(grid, metric, density)

# Individual linear models to obtain coefficients for each metric across
# all levels of density
extr_coeff_ci <- function(df, fm=NULL, ind=c(2, 8)) {
  #
  # Extract coef and ci for given formula of linear model
  #
  lmfit <- lm(fm, data=df)
  s <- summary(lmfit)
  return(list(coef(s)[ind[1]], coef(s)[ind[2]]*1.96))
}
fm <- formula("scale(slope) ~ scale(estimate) + scale(age) +
                 scale(sex) + scale(icv) + scale(bprs_tdbl) +
                   scale(as.numeric(as.factor(study)))")

#fm <- formula("scale(slope) ~ scale(estimate) + scale(age) +
#                 scale(sex) + scale(icv) + scale(bprs_tdbl)") 

cl <- by(pmsl, list(pmsl$metric, pmsl$density), extr_coeff_ci, fm,
         c(2, 9))

clm <- as.data.frame(do.call("rbind", cl)) %>%
  mutate(metric=rep(unique(pmsl$metric), 13),
         density=rep(unique(pmsl$density), each=11)) %>%
  rename(coeff=V1, ci=V2) %>%
  mutate(coeff=unlist(coeff),
         ci=unlist(ci)) %>%
  mutate(lower=coeff-ci,
         upper=coeff+ci) %>%
  filter(metric!="Cp")

# Produce a plot for coefficients and conf intervals showing how
# coefficients differ from zero when predicting treatment response
# across levels of density
g2 <- ggplot(clm, aes(x=density, y=coeff, fill="blue")) +
  geom_line() +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, fill="blue") +
  #facet_wrap(~metric, scale="free", labeller=as_labeller(graph_params)) +
  facet_wrap(~metric, labeller=as_labeller(graph_params_dict)) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_gray() +
  theme(
    #panel.border = element_blank(),
    #panel.grid.major.y = element_line(colour = "lightgray"),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor = element_blank(),
    plot.title=element_text(size=30, face="bold", hjust=0),
    #axis.line = element_line(colour = "black"),
    axis.line = element_blank(),
    axis.text.x=element_text(size=20, angle=NULL),
    axis.text.y=element_text(size=20),
    legend.text=element_text(size=20),
    legend.title=element_text(size=18, face="bold"),
    legend.position="right",
    #legend.title=element_blank(),
    strip.text.x=element_text(size=20, face="bold"),
    #strip.background = element_blank(),
    #legend.direction="horizontal" ,
    axis.title=element_text(size=30, face="bold")) +
  xlab(label="\nDensity") +
  ylab(label="Coefficient \n") +
  #coord_cartesian(ylim=c(0, 80)) + 
  ggtitle("") 

fn <- paste("../output/figures/fe_freesurfer_graphsumresp_",
            atlas, ".pdf", sep="")
ggsave(plot=g2,
       filename=fn,
       device=pdf, width=16, height=10)
rm(fn)
#-----------------------------------------------------------------------
# Boot strap for R2
#-----------------------------------------------------------------------
# Bootstrap 95% CI for R-Squared
# function to obtain R-Squared from the data 
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- plsr(formula, data=d, ncomp=2, validation="CV")
  pred <- predict(fit, d)
  pls_eval <- data.frame(obs=d$slope, pred=as.numeric(pred))
  ds       <- defaultSummary(pls_eval)
  #print(ds)
  return(ds[2])
 # 
  #return(summary(fit)$r.square)
} 
# bootstrapping with 1000 replications 
results <- boot::boot(data=subset(dms, density==0.1), statistic=rsq, 
  	R=1000, formula=slope ~ .)

# view results
results 
plot(results)

# get 95% confidence interval 
boot.ci(results, type="bca")

# Calculate P-value
brsq_pval <- mean(abs(mean(results$t)-results$t) > abs(results$t0))
fn <- paste("../output/tables/fe_msn_bootstrap_rsq_",
            atlas, ".csv", sep="")
write.csv(data.frame(brsq_pval=brsq_pval, density=0.1, niter=niter),
          fn,
          row.names=FALSE)
rm(fn)
#-----------------------------------------------------------------------
# Permutation test to assess significance of PLS model fit
#-----------------------------------------------------------------------
permut_rsq <- function(df, ncomp) {
  # Sample indices without replacement
  idx <- sample(nrow(df), replace=FALSE)
  df$y <- df$y[idx]
  plsfit <- plsr(y ~., data=df)
  rsq <- R2(plsfit)$val[2]
}

tmpdf <- dms %>% rename(y=slope) %>% filter(density==0.1) %>%
  dplyr::select(-density)
prsq <- sapply(1:npermut, function(x) permut_rsq(tmpdf))

# Calculate P-value
rsqorig <- R2(plsfitlist[[1]], ncomp=2)$val[2]
prsq_pval <- sum(prsq >= rsqorig)/npermut

fn <- paste("../output/tables/fe_msn_permut_rsq_",
            atlas, ".csv", sep="") 
write.csv(data.frame(rsq_orig=rsqorig, prsq=prsq, prsq_pval=prsq_pval,
                     density=0.1, niter=npermut),
          fn, row.names=FALSE)
rm(fn)

# Repeat this for all levels of density
for (i in sparsity) {
  tmpdf <- dms %>% rename(y=slope) %>% filter(density==i) %>%
    dplyr::select(-density)
  prsqlist <- lapply(1:npermut, function(x) permut_rsq(tmpdf))
  print(paste("P-val: ", sum(prsqlist[[i]] >= rsqorig)/npermut, sep=""))
}
#-----------------------------------------------------------------------
# Do the whole procedure for HC also
#-----------------------------------------------------------------------
groupshc <- list("hcm")

# Lists of data sets and slopes
all       <- list(hcm)

#-----------------------------------------------------------------------
# Extract grids for each group
#-----------------------------------------------------------------------
gridhc <- lapply(all,
               function(x) {
                 filter(x,
                        !is.na(sessnum) &
                        day==0) %>%
                   dplyr::select(grid)
               })
#-----------------------------------------------------------------------
# Estimate graphs 
#-----------------------------------------------------------------------
#atlasstr <- '_aparc.a2009s.annot'
#atlasstr <- ""
#nrois <- 148
#i <- 0
#graphlisthc <- lapply(groupshc,
#                 function(x) {
#                   i <<- i + 1 
#                   make_graphs(x, sparsity, nsp, gridhc[[i]]$grid,
#                               atlasstr, nrois)
#                 })
fn <- paste("../data/fe_graphlisthcm_", atlas, ".rds", sep="")
#saveRDS(graphlisthc, fn)
graphlisthc <- readRDS(fn)
rm(fn)


# Use only the brainGraph output, the first element for each group
grlisthc <- list(graphlisthc[[1]][[1]])

#-----------------------------------------------------------------------
# Extract degree for each hub in each subject and each sparsity
#-----------------------------------------------------------------------
degreelisthc <- lapply(grlisthc,
                     function(x) {
                       lapply(c(1:length(sparsity)),
                              function (y) {
                                calc_degree(x, spi=y)
                              })
                     })
#-----------------------------------------------------------------------
# Extract graph params for each subject and each sparsity
#-----------------------------------------------------------------------
#paramslisthc <- lapply(grlisthc,
#                     function(x) {
#                       lapply(c(1:length(sparsity)),
#                              function (y) {
#                                get_graph_params(x, y)
#                              })
#                     })
#saveRDS(paramslisthc, "../data/fe_paramslisthcm_destrieux.rds")
#paramslisthc <- readRDS("../data/fe_paramslisthcm_destrieux.rds")
fn <- paste("../data/fe_paramslisthcm_", atlas, ".rds", sep="")
#saveRDS(paramslisthc, fn)
paramslisthc <- readRDS(fn)
rm(fn)

#-----------------------------------------------------------------------
# Create data frame from the list of degrees
#-----------------------------------------------------------------------
Ddfhc <- do.call("rbind", do.call("rbind", degreelisthc)) %>%
  mutate(study=rep(rep(unlist(groupshc), each=nrow(hcm)), length(sparsity)),
         grid=rep(gridhc[[1]]$grid, length(sparsity)),
         density=rep(sparsity, each=nrow(hcm)))

# Merge with covariates
dmshc <- merge(Ddfhc, hcm %>%
                  dplyr::select(age, sex, grid), by="grid")

#-----------------------------------------------------------------------
# Transform paramslist to data frame
#-----------------------------------------------------------------------
Pdfhc       <- do.call("rbind", do.call("rbind", paramslisthc)) %>%
  mutate(study=rep(rep(unlist(groupshc), each=nrow(hcm)), length(sparsity)),
         grid=rep(gridhc[[1]]$grid, length(sparsity)))

pmshc       <- merge(Pdfhc, hcm %>% dplyr::select(grid, age, sex, iq,
                                                  icv), by="grid") 

pmslhc      <- pmshc %>% gather(key=metric, value=estimate,
                       Cp, Lp, E.global, E.local, diameter,
                       transitivity, num.hubs, vulnerability,
                       assortativity, smallworld, richclub)
#-----------------------------------------------------------------------
# Plot average curves for graph metrics
#-----------------------------------------------------------------------
pm1 <- pmsl %>% dplyr::select(-slope, -dup, -race,
                              -mccb, -bprs_tdbl)
pm2 <- pmslhc  
pmslall <- rbind(pm1, pm2) %>% rename(group=study)
pmslall$group[pmslall$group %in% c("cidar", "omega3")] <- "SZ"
pmslall$group[pmslall$group %in% c("hcm")] <- "HC"
pmslmall <- pmslall %>% group_by(density, metric, group) %>%
  dplyr::summarize(mean=mean(estimate, na.rm=TRUE),
                   sd=sd(estimate, na.rm=TRUE),
                   n=sum(!is.na(estimate)),
                   se=sd/sqrt(n),
                   ci=1.96*se) %>% filter(metric!="Cp")
                         
g1 <- ggplot(pmslmall, aes(x=density, y=mean, fill=group)) +
  geom_line(aes(col=group)) +
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.3) +
  #geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~metric, scale="free",
             labeller=as_labeller(graph_params_dict)) +
  theme_grey() +
  #theme_classic() +
  theme(
    #panel.border = element_blank(),
    #panel.grid.major.y = element_line(colour = "lightgray"),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor = element_blank(),
    plot.title=element_text(size=30, face="bold", hjust=0),
    #axis.line = element_line(colour = "black"),
    #axis.line = element_line(colour = "white"),
    axis.text.x=element_text(size=20, angle=NULL),
    axis.text.y=element_text(size=20),
    legend.text=element_text(size=20),
    #legend.title=element_text(size=18, face="bold"),
    legend.title=element_blank(),
    legend.position="right",
    #legend.title=element_blank(),
    strip.text.x=element_text(size=20, face="bold"),
    #strip.background = element_blank(),
    #legend.direction="horizontal" ,
    axis.title=element_text(size=30, face="bold")) +
  xlab(label="\nDensity") +
  ylab(label="Estimate\n") +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  #coord_cartesian(ylim=c(0, 80)) + 
  ggtitle("") 

fn <- paste("../output/figures/fe_freesurfer_graphsum_", atlas, ".pdf",
            sep="")
ggsave(filename=fn, plot=g1,
       width=16, height=10)
rm(fn)

graphsumlist <- by(pmslall, list(pmslall$metric, pmslall$density),
                   lm, formula=formula("estimate ~ group"))
graphsum_tvec <- sapply(graphsumlist, function(x) coef(summary(x))[6])
graphsum_pvec <- sapply(graphsumlist, function(x) coef(summary(x))[8])
#graphsum_dfs <- sapply(graphsumlist, function(x) coef(summary(x))[6])

# Get FDR corrected p value
graphsum_qindex <- FDR(graphsum_pvec)
graphsum_qvec <- rep(1, length(graphsum_pvec))
graphsum_qvec[graphsum_pvec <= graphsum_pvec[graphsum_qindex]] <- 0.05

graphsum_tdf <- data.frame(metric=rep(unique(pmslall$metric), 13),
                           density=rep(sparsity, each=11),
                           tval=graphsum_tvec,
                           pval=graphsum_pvec,
                           qval=graphsum_qvec) %>%
  filter(tval < -2 | tval > 2) %>%
  filter(metric!="Cp") %>%
  arrange(metric, density)
graphsum_tdf[, 3:4] <- round(graphsum_tdf[, 3:4], 3)
fn <- paste("../output/tables/fe_freesurfer_graphsum_", atlas,
            ".csv", sep="")
write.table(graphsum_tdf, fn,
            sep=",", col.names=TRUE, quote=TRUE, row.names=FALSE)

#-----------------------------------------------------------------------
# individual pages
#-----------------------------------------------------------------------
fn <- paste("../output/figures/fe_freesurfer_graphsum_individual_",
            atlas, ".pdf", sep="")
pdf(fn)
for (i in 1:length(unique(pmslmall$metric))) {
  tmpdf <- pmslmall %>% filter(metric==unique(pmslmall$metric)[i])
  print(ggplot(tmpdf, aes(x=density, y=mean, fill=group)) +
    geom_line(aes(col=group)) +
    geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.3) +
    #geom_hline(yintercept=0, linetype="dashed") +
    facet_wrap(~metric, scale="free",
               labeller=as_labeller(graph_params_dict)) +
    theme_grey() +
    #theme_classic() +
    theme(
      #panel.border = element_blank(),
      #panel.grid.major.y = element_line(colour = "lightgray"),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor = element_blank(),
      plot.title=element_text(size=30, face="bold", hjust=0),
      #axis.line = element_line(colour = "black"),
      #axis.line = element_line(colour = "white"),
      axis.text.x=element_text(size=20, angle=NULL),
      axis.text.y=element_text(size=20),
      legend.text=element_text(size=20),
      #legend.title=element_text(size=18, face="bold"),
      legend.title=element_blank(),
      legend.position="right",
      #legend.title=element_blank(),
      strip.text.x=element_text(size=20, face="bold"),
      #strip.background = element_blank(),
      #legend.direction="horizontal" ,
      axis.title=element_text(size=30, face="bold")) +
    xlab(label="\nDensity") +
    ylab(label="Estimate\n") +
    scale_color_manual(values = c("#00BFC4", "#F8766D")) +
    scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
    #coord_cartesian(ylim=c(0, 80)) + 
    ggtitle("")) 
}
dev.off()

#-----------------------------------------------------------------------
# Test interaction of group and IQ on graph metrics 
#-----------------------------------------------------------------------
pmslist <- by(pmslall, list(pmslall$metric, pmslall$density),
                 lm, formula=formula("scale(estimate) ~
                     scale(age) + scale(sex) + scale(icv) +
                     scale(iq) * scale(as.numeric(as.factor(group)))"))


#pmslhclist <- by(pmsl, list(pmsl$metric, pmsl$density),
#                 lm, formula=formula("scale(estimate) ~
#                                        scale(mccb_reasoning)"))

#pmslhclist <- by(pmslall, list(pmslall$metric, pmslall$density),
#                 lm, formula=formula("scale(estimate) ~ 
#                     scale(age) + scale(sex) +
#                     scale(iq) * scale(as.numeric(as.factor(group)))"))

# Extract t-value
cvec_iq <- sapply(pmslist, function(x) coef(summary(x))[7])
svec_iq <- sapply(pmslist, function(x) coef(summary(x))[14])
tvec_iq <- sapply(pmslist, function(x) coef(summary(x))[21])
pvec_iq <- sapply(pmslist, function(x) coef(summary(x))[28])

# Build data frame for results, sort by metric and density
co_iq <- data.frame(metric=rep(sort(unique(pmslall$metric)), 13),
                     density=rep(sparsity, each=11), t=tvec_iq,
                     coef=cvec_iq, ci=1.96 * svec_iq,
                     psig=starsfromp(pvec_iq)) %>%
  arrange(metric, density)

p <- plot_summaries(co_iq %>% rename(y=coef, x=density,
                                      wrap=metric, se=ci) %>%
                    mutate(group="blue"), graph_params_dict, xlab="Density",
                    ylab="Coefficient", scale=NULL)

p + geom_hline(yintercept=0, linetype="dashed")
#-----------------------------------------------------------------------
# Univariate cortical thickness group comparison
#-----------------------------------------------------------------------
hcmrois <- hcm %>% dplyr::select(matches(paste(rois$name,
                                               collapse="|"))) %>%
    bind_cols(hcm %>% dplyr::select(grid, day, age, sex, icv)) %>%
    mutate(study="HC")  %>%
    mutate(group="HC") 

fefrois <- fef %>% dplyr::select(matches(paste(rois$name,
                                               collapse="|"))) %>% 
    bind_cols(fef %>% dplyr::select(grid, day, age, sex, icv, study)) %>%
    mutate(group="SZ") 

feorois <- feo %>% dplyr::select(matches(paste(rois$name,
                                               collapse="|"))) %>% 
    bind_cols(feo %>% dplyr::select(grid, day, age, sex, icv, study)) %>%
    mutate(group="SZ") 

allrois <- rbind(fefrois, feorois, hcmrois) %>% filter(day==0) %>% na.omit()

# Residualize for age, sex, and ICV
allroislist <- by(allrois, allrois$study, residualize_df,
                     formula("y ~ age * sex"),
                     1:68)

allroisr <- do.call("rbind", allroislist)

# remove row names
rownames(allroisr) <- NULL

allroisl <- allrois %>% gather(key=roi,
                               value=thickness,
                               paste("l", rois$name, sep="_"),
                               paste("r", rois$name, sep="_"))

# Run linear model for each region
lmroislist <- by(allroisl, allroisl$roi, lm,
                 formula=formula("thickness ~ age * sex +
                                     icv + group")) 

# Extract t-values for the effect of group
# NB: negative t-value means lower thickness in SZ compared to HC
tvec <- sapply(lmroislist, function(x) coef(summary(x))[17])

# Bar plot for all t-values reaching uncorrected significance P < 0.05
barplot(tvec[tvec < -2], las=1, cex.name=0.8)

# Data frame from t values
tval_df <- tvec %>% as.data.frame() %>%
    rename(t=".") %>%
    mutate(roi=rownames(as.data.frame(tvec)))

tval_df$Pval <- 2*pt(abs(tval_df$t), nrow(tval_df), lower=FALSE)
qindex <- FDR(tval_df$Pval)
sigindex <- which(tval_df$Pval <= tval_df$Pval[qindex])
tval_df$q <- 1
tval_df$q[sigindex] <- 0.05
if (!is.null(qindex)) {
  tval_df$qth <- tval_df$t[qindex]
}
tval_df$tsig <- 0
tval_df$tsig[sigindex] <- tval_df$t[sigindex]


# Sort according to freesurfer
tval_df$roi <- factor(tval_df$roi, levels=c(paste("l", rois$name, sep="_"),
                                            paste("r", rois$name, sep="_")))

tval_df <- tval_df[order(tval_df$roi),]


# for manuscript
fn <- paste("../output/tables/fe_freesurfer_roitvals_",
            atlas, ".csv", sep="")
write.table(tval_df,
            file=fn,
            sep=",", row.names=FALSE)
rm(fn)

# for python input
fn <- paste("../output/tables/fe_freesurfer_pyinput_roitvals_",
            atlas, ".csv", sep="")
write.table(tval_df$t,
            file=fn,
            sep=",", row.names=FALSE, col.names=FALSE)
rm(fn)
#-----------------------------------------------------------------------
# Group average graphs
#-----------------------------------------------------------------------
fefmatlist <- graphlist[[1]][[2]]
fefavgmat <- average_mat(fefmatlist)
fefbavgmat <- binarize_mat(fefavgmat, 0.1)
fefavg <- graph_from_mat(fefbavgmat, atlas=atlas)
feomatlist <- graphlist[[2]][[2]]
feoavgmat <- average_mat(feomatlist)
feobavgmat <- binarize_mat(feoavgmat, 0.1)
feoavg <- graph_from_mat(feobavgmat, atlas=atlas)
fefeoavgmat <- (fefavgmat + feoavgmat)/2
fefeobavgmat <- binarize_mat(fefeoavgmat, 0.1)
fefeoavg <- graph_from_mat(fefeobavgmat, atlas=atlas)
fefeoavg <- set_brainGraph_attr(fefeoavg)
degree(fefeoavg)

hcmmatlist <- graphlisthc[[1]][[2]]
hcmavgmat <- average_mat(hcmmatlist)
hcmbavgmat <- binarize_mat(hcmavgmat, 0.1)
hcmavg <- graph_from_mat(hcmbavgmat, atlas=atlas)
hcmavg <- set_brainGraph_attr(hcmavg)
sum(degree(hcmavg))
sum(degree(fefeoavg))
barplot(degree(hcmavg), las=2)
barplot(degree(fefeoavg), las=2, add=TRUE, name="", col="blue")

delta <- (degree(fefeoavg) - degree(hcmavg))

deltap <- array(dim=c(npermut, nrois))
for (i in 1:npermut) {
    deltap[i, ] <- delta[sample(length(delta))]
}

sum(abs(colMeans(deltap)) >= delta)/npermut   

pvals_delta<- sapply(1:nrois,
                     function(x) {
                       sum(abs(deltap[, x]) >= abs(delta[x]))/npermut
                     })

# Get ROI names of significant differences
#dk$name[which(pvals_delta<=0.05)]
destrieux$name[which(pvals_delta<=0.05)]
delta[which(pvals_delta<=0.05)]
qdeltaindex <- FDR(pvals_delta)

# Map delta degrees to roi list
fn <- paste("../output/tables/fe_freesurfer_pyinput_deltadegrees_",
            atlas, ".csv", sep="")
write.table(data.frame(delta=delta),
            fn,
          col.names=FALSE,
          row.names=FALSE)
rm(fn)

# Write significant values to disk
sigvals <- vector("numeric", nrois)
sigvals[which(pvals_delta<=0.05)] <- delta[which(pvals_delta<=0.05)]
fn <- paste("../output/tables/fe_freesurfer_pyinput_deltadegrees_sig_",
            atlas, ".csv", sep="")
write.table(sigvals, fn,
            col.names=FALSE, row.names=FALSE)
rm(fn)

# Write results to disk
fn <- paste("../output/tables/fe_freesurfer_deltadegrees_",
            atlas, ".csv", sep="")
write.csv(data.frame(region=dk$name[which(pvals_delta<=0.05)],
                     pval=pvals_delta[which(pvals_delta<=0.05)],
                     delta=delta[which(pvals_delta<=0.05)]),
                     fn,
                     row.names=FALSE)
rm(fn)

#-----------------------------------------------------------------------
# Boot strap for PLS
#-----------------------------------------------------------------------
# Bootstrap 95% CI for R-Squared
# function to obtain R-Squared from the data 
boot_plsfit <- function(formula, data, indices) {
  d <- data[indices, ] # allows boot to select sample 
  fit <- plsr(formula, data=d, ncomp=2, validation="CV")
  #pred <- predict(fit, d)
  #pls_eval <- data.frame(obs=d$slope, pred=as.numeric(pred))
  #ds       <- defaultSummary(pls_eval)
  #print(ds)
  bld <- loadings(fit)
  bld1 <- bld[, 1]
  bld2 <- bld[, 2]
  #print(bld2)
  #return(list(ld[, 1], ld[, 2]))
 # 
  #return(summary(fit)$r.square)
  return(c(bld1, bld2))
} 

bpl <- lapply(1:npermut,
              function(x) {
                boot_plsfit(formula=slope ~ .,
                            data=subset(dms, density==0.1),
                            indices=sample(nrow(dms), replace=TRUE))
              })

bpls <- do.call("rbind", bpl)
#-----------------------------------------------------------------------
# Use the PLS scores to determine contributions of brain regions to
# treatment response
#
# This is modeled after Seidlitz et al. 2018, PNAS. From the scripts
# used int that study, it is obvious that what is used for the brain
# visualization are the weights. 
#-----------------------------------------------------------------------
ld <- loadings(plsfitlist[[1]])
ld1 <- ld[, 1]

# Standardized weights
ld1z <- ld1/sd(bpls[, 1:69])

ld2 <- ld[, 2]

# Standardized weights
ld2z <- ld2/sd(bpls[, 17:138])

barplot(sort(ld1z, decreasing=FALSE), las=2)
barplot(sort(ld2z, decreasing=FALSE), las=2)

perc.rank <- function(x) trunc(rank(x))/length(x)
ld12p <- list(perc.rank(ld1z), perc.rank(ld2z))

fname1 <- paste("../output/tables/fe_freesurfer_pyinput_plsweights1_",
                atlas, ".csv", sep="")
fname2 <- paste("../output/tables/fe_freesurfer_pyinput_plsweights2_",
                atlas, ".csv", sep="")
fnames <- c(fname1, fname2)
sapply(1:2,
       function(x) {
         write.table(ld12p[[x]], fnames[x],
                     row.names=FALSE, col.names=FALSE, sep=",")
       })
rm(fnames)
rm(fname1)
rm(fname2)

comps <- c("PLS1", "PLS2")
complist <- lapply(comps,
                  function(x) {
                    plot_components(plsfitlist[[1]],
                      subset(dms, density==sparsity[1]),
                      comp=x)
                  })

# Extract correlations for both components
pls_scores <- scores(plsfitlist[[1]])
rpls1 <- (cor.test(pls_scores[, 1], dms$slope[dms$density==0.1]))
rpls2 <- (cor.test(pls_scores[, 2], dms$slope[dms$density==0.1]))

fn <- paste("../output/tables/fe_msn_plsscore_rs_", atlas,
            ".csv", sep="") 
write.csv(data.frame(r=c(rpls1$estimate, rpls2$estimate),
                     pval=c(rpls1$p.val, rpls2$p.val),
                     component=c("PLS1", "PLS2"),
                     density=0.1),
          fn, row.names=FALSE)
rm(fn)

# Run python and create brains
system(paste("../src/fe_brains.py -a", atlas))

fname1 <- paste("../output/figures/fe_freesurfer_pls1brain_",
                atlas, ".png", sep="")
fname2 <- paste("../output/figures/fe_freesurfer_pls2brain_",
                atlas, ".png", sep="")

imgs <- c(fname1, fname2)
rm(fname1)
rm(fname2)

imglist <- lapply(imgs, function(x) readPNG(x))

rglist <- lapply(imglist, function(x) rasterGrob(x, interpolate=TRUE))
p <- plot_grid(complist[[1]], rglist[[1]],
               complist[[2]], rglist[[2]], labels='AUTO', scale=1.0,
               label_size=35)
fn <- paste("../output/figures/fe_freesurfer_plsbrains_",
            atlas, ".pdf", sep="")
ggsave(plot=p, filename=fn, height=13, width=15)
rm(fn)
