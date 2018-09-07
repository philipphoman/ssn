#-----------------------------------------------------------------------
# FE functions
#
# PH, 2/13/18
#-----------------------------------------------------------------------
#
#-----------------------------------------------------------------------
# Libraries 
#-----------------------------------------------------------------------
#devtools::install_github("sjmgarnier/viridis")

if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(
          "R.matlab",
          "lm.beta",
          "cowplot",
          "mixtools",
          "ellipse",
          "car",
          "magick",
          "png",
          "grid",
          "flexmix",
          "DescTools",
          "boot",
          "lme4",
          "lsmeans",
          "graphics",
          "Hmisc",
          "tibble",
          "corrplot",
          "astsa",
          "brainGraph",
          "tidyr",
          "ggplot2",
          "dplyr",
          "MASS",
          "caret",
          "pls" 
        )

#-----------------------------------------------------------------------
# Functions 
#-----------------------------------------------------------------------
anon_grid <- function(df) {
  #
  # Add a anonymized grid variable
  # for subjects with MRI. The rest
  # keeps original grid
  #
  grid <- df %>% filter(hasmri==1&day==0) %>% dplyr::select(grid)
  i <- c(1:length(grid$grid))
  dftop <- merge(df %>% filter(hasmri==1),
                 data.frame(grid=grid, grid_anon=i))
  dfbottom <- df %>% filter(is.na(hasmri)) %>% mutate(grid_anon=grid+1000)
  return(rbind(dftop, dfbottom))
}

starsfromp <- function(pval, c1="~", c2="*") {
  # Parses pvalues, returns asterisks

  as <- vector(mode="character", length=length(pval))
  for (i in 1:length(pval)) {
    p <- pval[i]
    if (p < 0.1) as[i] <- c1 
    if (p < 0.05) as[i] <- paste(c2, sep="")
    if (p < 0.01) as[i] <- paste(c2, c2, sep="")
    if (p < 0.001) as[i] <- paste(c2, c2, c2, sep="")
  }
  return(as)
}

clean_paramslist <- function(paramslist, ind) {
  #
  # This function removes subjects
  # from a list of graph parameters
  pr <- paramslist
  pr[[1]] <- lapply(1:13, function(x) paramslist[[1]][[x]][-ind, ])
  return(pr)
}

clean_graphlist <- function(graphlist, ind) {
  #
  # This function removes subjects
  # from a list of graph parameters
  #l <- lapply(1:13, function(x) graphlist[[1]][[1]][[x]][-ind, ])
  gr <- graphlist
  gr[[1]][[1]] <- graphlist[[1]][[1]][-ind]
  return(gr)
}

set_mriflag <- function(df) {
  #
  # Sets custom flag for dataframe
  #
  dfnew <- df %>%
    dplyr::select(grid, sessnum, day,
                  Left.Amygdala) %>%
    dplyr::filter(!is.na(sessnum),
                  !is.na(Left.Amygdala),
                  day==0) %>%
    mutate(hasmri=1) %>%
    dplyr::select(-sessnum, -Left.Amygdala, -day) %>%
    right_join(df, by="grid") 
   return(dfnew)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Helper function to make a data-frame of ellipse points that 
# includes the level as a column
make_ellipse <- function(cov_mat, center, level) {
  ellipse::ellipse(cov_mat, centre = center, level = level) %>%
    as.data.frame() %>%
    add_column(level = level) %>% 
    as_tibble()
}

plot_average <- function(df, xlab="Week", legend.position="right",
                         ylab="Thinking disturbance BPRS") {
  #
  # Plot average treatment effect over 12 weeks
  #
  dfm <- df %>% group_by(study, visit) %>%
    dplyr::summarize(mean=mean(y, na.rm=TRUE),
                     sd=sd(y, na.rm=TRUE),
                     n=sum(!is.na(y)),
                     se=sd/sqrt(n),
                     ci=1.96 * se)

  pd <- position_dodge(0.5)
  g <- ggplot(dfm, aes(x=visit, y=mean, color=study)) +
    geom_line(size=1, position=pd) +
    geom_point(size=3, position=pd) +
    geom_errorbar(size=1, 
                  aes(ymin=mean-ci, ymax=mean+ci), width=.0,
                 position=pd) +
    theme_classic() +
    theme(
      axis.title=element_text(size=20, face="bold"),
      axis.text=element_text(size=18),
      legend.text=element_text(size=14),
      legend.title=element_blank(),
      legend.position=legend.position
    ) +
    scale_color_manual(values=c("dodgerblue4", "dodgerblue1")) +
    xlab(xlab) +
    ylab(ylab)
  return(g)
}

plot_summaries <- function(df, params, scale="free",
                           base_size=14,
                           xlab="Density", ylab=NULL,
                           legend.position="right") {
  #
  # Plot summary graphs across
  # densities
  #
  #print(df)
  p <- ggplot(df, aes(x=x, y=y, fill=group)) +
    geom_line(aes(col=group)) +
    geom_ribbon(aes(ymin=y-se, ymax=y+se), alpha=0.3) +
    facet_wrap(~wrap, scale=scale, labeller=as_labeller(params)) +
    theme_grey(base_size=base_size) +
    theme(
      plot.title=element_text(face="bold", hjust=0),
      axis.text.x=element_text(angle=NULL),
      legend.title=element_blank(),
      legend.position=legend.position,
      strip.text.x=element_text(face="bold"),
      axis.title=element_text(face="bold")) +
    xlab(label=paste("\n", xlab)) +
    ylab(label=paste(ylab, "\n")) +
    scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  ggtitle("") 
  return(p)
}


  
plot_shrinkage <- function(df, fm, xlab="Log (day)",
                           ylab="Thinking disturbance BPRS",
                           xlab2="Intercept",
                           ylab2="Slope") {
  #
  # Shrinkage plot
  #

  # Show individual measurements for each subject
  # plus a regression line.
  # Slope of the regression line could be seen as
  # indivgridual response (without shrinkage)
  g1 <- ggplot(df %>% filter(hasmri==1)) +
    aes(x=x, y=y) +
    stat_smooth(method="lm", se=FALSE) +
    geom_point() +
    facet_wrap(~grid) 
    #labs(x=xlab, y=ylab) 
    #scale_x_continuous(breaks = 0:4 * 2)
  #print(g1)

  # Calculate the corresponding linear models
  df1 <- lmList(y ~ x | grid, df) %>% 
    coef() %>% 
    as_tibble() %>%
    rownames_to_column("grid") %>% 
    dplyr::rename(Intercept=`(Intercept)`, Slope=x) %>% 
    add_column(Model="No shrinkage") 
  #df1$grid <- factor(df1$grid)

  # Now calculate a mixed model to obtain shrinkage
  # estimates
  lmer.fit <- lmer(y ~ x + (1 + x|grid), df)

  # Calculate the mixed model again with MRI-only
  # subjects
  lmer.fit.red <- lmer(formula(lmer.fit),
                       data=df %>% filter(hasmri==1))

  # Extract the matrix
  cov_mat <- VarCorr(lmer.fit)[["grid"]]

  # Strip off some details so that just the useful part is printed
  attr(cov_mat, "stddev") <- NULL
  attr(cov_mat, "correlation") <- NULL
  cov_mat
  #>             (Intercept)     Days
  #> (Intercept)   582.69656  9.89797
   #> Days            9.89797 35.03298


  center <- lme4::fixef(lmer.fit)
  levels <- c(.1, .3, .5, .7, .9)

  # Create an ellipse dataframe for each of the levels defined 
  # above and combine them
  df_ellipse <- levels %>%
    purrr::map_df(~ make_ellipse(cov_mat, center, level = .x)) %>% 
    dplyr::rename(Intercept = `(Intercept)`, Slope=x)


  # Produce a dataframe with shrinkage estimates
  df2 <- coef(lmer.fit)[["grid"]] %>% as_tibble() %>% 
    rownames_to_column("grid") %>% 
    dplyr::rename(Intercept=`(Intercept)`, Slope=x) %>% 
    add_column(Model="Shrinkage")
  #df2$grid <- factor(df2$grid)

  # Produce a dataframe with shrinkage estimates from reduced data set
  df2.red <- coef(lmer.fit.red)[["grid"]] %>% as_tibble() %>% 
    rownames_to_column("grid") %>% 
    dplyr::rename(Intercept=`(Intercept)`, Slope=x) %>% 
    add_column(Model="Shrinkage (reduced)")


  # Merge with df1
  dfs <- bind_rows(df1, df2) %>% left_join(df, by="grid")

  # Filter for MRI only
  dfsf <- dfs %>% filter(hasmri==1)

  # Restrict to available slopes
  completegrids <- subset(df1, select=grid, !is.na(Slope))

  dfsf2 <- dfs %>% filter(grid %in% completegrids$grid,
                          hasmri==1,
                          x==0) %>%
    dplyr::select(grid, Intercept, Slope, Model)
    

  # Plot
  g2 <- ggplot(dfsf %>% mutate(grid=factor(grid, levels=c(1:41)))) + 
    aes(x=x, y=y) + 
    geom_abline(aes(intercept=Intercept, slope=Slope, color=Model),
                size = .75) + 
    geom_point() +
    facet_wrap(~grid) +
    scale_color_manual(labels = c("No pooling",
                                  "Partial pooling"),
                       values = c("#F8766D", "#00BFC4")) +
   theme_bw() +
    theme(
      #panel.border = element_blank(),
      #panel.grid.major.y = element_line(colour = "lightgray"),
      panel.grid.major.y = element_line(colour = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title=element_text(size=30, face="bold", hjust=0),
      #axis.line = element_line(colour = "white"),
      axis.text.x=element_text(size=14, angle=NULL),
      axis.text.y=element_text(size=14),
      legend.text=element_text(size=34, face="bold"),
      legend.position="top",
      legend.title=element_blank(),
      strip.text.x=element_text(size=26, face="bold"),
      #strip.background = element_blank(),
      #legend.direction="horizontal" ,
      axis.title=element_text(size=40, face="bold")) +
    xlab(label=paste("\n", xlab)) +
    ylab(label=paste(ylab, "\n")) +
    #coord_cartesian(ylim=c(0, 80)) + 
    ggtitle("") 

  # Show the shrinkage effect
  df3 <- data_frame(
    Model="Average",
    Intercept=fixef(lmer.fit)[1],
    Slope=fixef(lmer.fit)[2])


  # Contour plot 
  g3 <- plot_contour(dfsf2,
                     dfsf2 %>% filter(Model=="No shrinkage"),
                     df3, df_ellipse, xlab=xlab2, ylab=ylab2) +
  scale_color_manual(labels = c("Average",
                                "No pooling",
                                "Partial pooling"),
                       #values = c("#619CFF", "#F8766D", "#00BFC4")) 
                       values = c("darkblue", "#F8766D", "#00BFC4")) 

  # For contour plot with shrinkage and reduced shrinkage
  dfs <- bind_rows(df2.red, df2) %>% left_join(df, by="grid")

  # Filter for MRI only
  dfsf <- dfs %>% filter(hasmri==1)

  # Restrict to available slopes
  completegrids <- subset(df1, select=grid, !is.na(Slope))

  dfsf2 <- dfs %>% filter(grid %in% completegrids$grid,
                          hasmri==1,
                          x==0) %>%
    dplyr::select(grid, Intercept, Slope, Model)
  
  g4 <- plot_contour(dfsf2,
                     dfsf2 %>% filter(Model=="Shrinkage (reduced)"),
                     df3, df_ellipse, xlab=xlab2, ylab=ylab2) +
  scale_color_manual(labels = c("Average",
                                "Partial pooling (full data set)",
                                "Partial pooling (reduced data set)"),
                       #values = c("#619CFF", "#00BFC4", "darkred")) 
                       values = c("darkblue", "#00BFC4", "darkred")) 
  
  return(list(g1, g2, g3, g4))


}

plot_contour <- function(df1, df2, df3, df_ellipse,
                         xlab="Intercept",
                         ylab="Slope") {
  #
  # Produce contour plot
  # for shrinkage estimated
  #

  #df <- rbind(df1, df2)
  
  g <- ggplot(df1) + 
    aes(x=Intercept, y=Slope, color=Model) + 
    geom_point(size=6) + 
    geom_point(data=df3, size=8) + 
    geom_path(data=df1, aes(group=grid, color=NULL), 
              arrow=arrow(length=unit(.02, "npc"), ends="last")) + 
    # Jitter the labels away from the points
    ggrepel::geom_text_repel(size=10,
               aes(label=grid, color=NULL), 
               data=df2) +
    #aes(x=Intercept, y=Slope, color=Model) +

    # Draw contour lines from the distribution of effects
    geom_path(aes(group=level, color=NULL), data=df_ellipse,
              linetype="dashed", color="grey40") +
    #scale_color_brewer(palette = "Dark2")  +

    theme_gray() +
    theme(#panel.border = element_blank(),
      #panel.border = element_rect(colour = "black", fill=NA, size=1),
      #panel.grid.major.y = element_line(colour = "white"),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor = element_blank(),
      plot.title=element_text(size=30, face="bold", hjust=0),
      #axis.line = element_line(colour="black", size=0.1),
      axis.text.x=element_text(size=28, angle=NULL),
      axis.text.y=element_text(size=28),
      axis.title=element_text(size=40, face="bold"), 
      legend.text=element_text(size=34, face="bold"),
      legend.position="top",
      legend.title=element_blank(),
      strip.text.x=element_text(size=22, face="bold"),
      strip.background = element_blank()) +
      #legend.direction="horizontal" ,
    xlab(label=paste("\n", xlab)) +
    ylab(label=paste(ylab, "\n")) +
    #coord_cartesian(ylim=c(0, 80)) + 
    ggtitle("") 
  return(g)
}


plot_example_density <- function(d1, x1, d2, x2, nsp=2^9) {
  #
  # Plot examples
  pdf("../lib/fe_freesurfer_kl.pdf")
  par(oma=c(5, 5, 5, 5))
  plot(x=x1, y=d1, "l", col="violet", lwd=7, ylim=c(0, 0.8),
       xlab="", ylab="", bty="n", cex=2, cex.axis=2)
  lines(x=x2, y=d2, col="darkblue", lwd=7)
  #axis(1, cex=2)
  #axis(2, cex=2)
  mtext("Cortical thickness", 1, 3, font=2, cex=2.5)
  mtext("Density", 2, 3, font=2, cex=2.5)
  legend("topright", legend=c("P", "Q"),
         cex=1.5, lwd=c(7, 7), col=c("violet", "darkblue"), bty="n")
  dev.off()
}  



estimate_density <- function(group, grid, nsp=2^9,
                             atlasstr="", nrois=68){
  # Create matrix of hemispheres
  # with estimated densities per region
  pwd <- getwd()
  setwd("~/tmp")
  print(paste("Estimating density for ", grid, " with atlas ", atlasstr,
              sep=""))
  fn1 <- paste(group, "_vertices/", grid, "_lh", atlasstr,  ".mat", sep="")
  fn2 <- paste(group, "_vertices/", grid, "_rh", atlasstr,  ".mat", sep="")
  lh <- readMat(fn1)
  rh <- readMat(fn2)

  # Note: The order of the vertices is the one from Freesurfer
  # Need to apply the same order as used in brainGraph
  l <- lh$vertices
  r <- rh$vertices
  if (atlasstr == "") {
    atlasdf <- dk
  } else {
    atlasdf <- destrieux
    ind <- atlasdf$ind
    inds <- order(ind)
    l <- l[inds[1:74]]
    r <- r[inds[1:74]]
  }
  
    # Remove list entries 43 which include empty regions in Desistrieux
  # atlas
  if (atlasstr == "_aparc.a2009s.annot") { 
    #l[[43]] <- NULL
    #r[[43]] <- NULL
    #print("ok")
  # Remove list entries 1 and 5 which include empty regions
  # otherwise
  } else {
    l[[5]] <- NULL
    l[[1]] <- NULL
    r[[5]] <- NULL
    r[[1]] <- NULL
  }

  # merge lists of hemispheres
  lr <- cbind(l, r)

  # estimate densities for each region
  # there needs to be a fallback for missing values
  # for regions without data
  lrl <- lapply(lr,
                function(x) {
                  if (length(x[[1]]) == 0) {
                    x <- NA
                  } else {
                    x <- density(x[[1]], n=nsp)
                    x$y
                  }
                })

  #l1 <- lapply(l,
  #              function(x) {
  #                x <- density(x[[1]], n=nsp)
  #                x$y
  #              })

  #r1 <- lapply(r,
  #              function(x) {
  #                x <- density(x[[1]], n=nsp)
  #                x$y
  #              })


  #ml <- matrix(unlist(l1), ncol=34, nrow=nsp)
  #rl <- matrix(unlist(r1), ncol=34, nrow=nsp)

  #mm2 <- cbind(ml, rl)
  

  # make matrix from list
  mm <- matrix(unlist(lrl), ncol=nrois, nrow=nsp)
  # print("done\n")
  setwd(pwd)
  return(mm)
}

calc_kl_div <- function(M) {
  # Calculate KL divergence
  Mkl <- PairApply(M, FUN=function(x, y){
    exp(-sum(KLdiv(cbind(x, y))))})
  return(Mkl)
}

binarize_mat <- function(M, sparsity) {
  # How many edges with that sparsity level?
  #print(sparsity)
  nEdges <- ceiling(sparsity * dim(M)[1] * (dim(M)[1]-1))
  #print(nEdges)
  if (nEdges %% 2 == 1)
    nEdges <- nEdges + 1

  # Subject specific threshold
  Mv <- as.vector(M)
  Mv <- sort(Mv, decreasing=TRUE)
  thr <- Mv[nEdges]
  #Mt <- matrix(0, ncol=dim(M)[1], nrow=dim(M)[1])
  Mbin <- M
  
  Mbin[M>thr] <- 1
  Mbin[M<=thr] <- 0

  # Set self connection to zero (cf. Wang et al.)
  diag(Mbin) <- 0
  return(Mbin)
}

density_plot <- function(d1, d2) {
  #
  # Create density plot
  # To be used in Figure
  pdf("../lib/fe_freesurfer_density.pdf")
  layout(matrix(c(1:4), nrow=2, byrow=TRUE))
  par(oma=c(2, 2, 2, 2))
  plot(d1, col="violet", lwd=3, bty="n", xlim=c(0, 5),
       xlab="Cortical thickness", font.lab=2,
       ylab="Density", cex.lab=1.5, font=2, main="")
  #polygon(c( x[x>=1250], 1250 ),  c(y[x>=1250],0 ), col="red")
  #polygon(d1$x, d1$y,  col="darkred", border="darkred")
  lines(d2, col="darkblue", lwd=3)
  #polygon(d2$x, d2$y,  col="orange", border="orange")
  op <- par(font=2)
  legend("topright", col=c("violet", "darkblue"),
         lty=c(1, 1), cex=1, 
         lwd=c(2, 2), c("left IPL", "left OG"), bty="n")
  par(op)
  dev.off()

}
  
  

corr_plot <- function(M, cl.pos=0.5, col=NULL,
                      cl.lim=c(0, 1)) {
  #
  # Wrapper for corr plot
  #
  jet.colors <- colorRampPalette(c("#00007F",
                                   "blue", "#007FFF",
                                   "cyan", "#7FFF7F",
                                   "yellow", "#FF7F00",
                                   "red", "#7F0000"))

  jc <- jet.colors(100) 
  #col <- ifelse(is.null(col), c(jc, jc), col)
  if(is.null(col)) {
      col <- c(jc, jc)
  }
  p <- par(cex=1.5)
  corrplot::corrplot(M,
                     is.corr=FALSE,
                     method="color",
                     col=col,
                     cl.lim=cl.lim,
                     tl.pos="n",
                     cl.pos=cl.pos)
  par(p)
  #filled.contour(volcano, color = jet.colors, asp = 1, nlevels=100)
}

validation_plot <- function(plsfit, val.type="RMSE", col="black",
                            main=NULL) {
  validationplot(plsfit, lwd=3,
                 val.type=val.type, bty="n", main="",
                 cex.lab=1.5, font.lab=2, col=col)
  mtext(main, 3, cex=1.5, font=2)
}
  

make_graphs <- function(group, sparsity, nsp=2^9, grid, atlasstr="",
                        nrois=68) {
  #
  # Calculate KL-matrix for
  # each subject
  #
  pwd <- getwd()
  setwd("~/tmp")
  lgs <- vector("list", length(grid))
  Ms  <- vector("list", length(grid))
  bMs  <- vector("list", length(grid))

  if (atlasstr=="") {
    atlasdf <- dk
  } else {
    atlasdf <- destrieux
  }
  

  # Create list of connectivity matrices
  Ms <- lapply(grid,
           function(x) {
             x <- estimate_density(group, x, nsp=2^9, atlasstr, nrois)
             x <- as.data.frame(x)
             #colnames(x) <- dk$name
             colnames(x) <- atlasdf$name
             x <- calc_kl_div(x)
             diag(x) <- 0
             return(x)
           })

  # Binarize the matrices for each level of sparsity
  bMs <- lapply(Ms,
                function(x) {
                  lapply(sparsity,
                         function(y){
                           binarize_mat(x, y)
                         })
                })

  # Create a list of brain graphs from binarized connectivity matrices
  switch(atlasstr,
         "_aparc.a2009s.annot" = {at <- "destrieux"},
         at <- "dk")
           
  lgs <- lapply(bMs,
                function(x) {
                  lapply(c(1:length(sparsity)),
                         function(y) {
                           x <- graph_from_adjacency_matrix(x[[y]],
                                                            diag=FALSE,
                                                            mode="undirected")
                           make_brainGraph(x, atlas=at)
                         })
                  #return(x)
                })

  # Create average matrices
  avM <- Reduce("+", Ms) / length(Ms)
  
  setwd(pwd)
  return(list(lgs, Ms, avM))
}

average_mat <- function(Matlist) {
  #
  # Create average matrix from list of matrices
  #
  avM <- Reduce("+", Matlist) / length(Matlist)
  return(avM)
}



calc_degree <- function(graphlist, ncol=68, norm=FALSE, spi=1) {
  # Assign graph degrees to matrix
  # spi, sparsity index
  l <- sapply(graphlist,
         function(x) {
           as.numeric(degree(x[[spi]], normalized=norm))
         })
  X <- matrix(l, ncol=ncol, nrow=length(graphlist))
  
  #for (i in 1:nsubs){
  #  #X[i, ] <- as.numeric(degree(network[[5]][[i]][[spi]], normalize=norm))
  #  X[i, ] <- as.numeric(degree(network[[5]][[i]][[spi]], normalize=norm))
  #}
  return (as.data.frame(X))
}

graph_from_mat <- function(M, atlas="dk") {
  #
  # Create graph from matrix
  g <- graph_from_adjacency_matrix(M, diag=FALSE, mode="undirected")
  g <- make_brainGraph(g, atlas=atlas)
  g <- set_brainGraph_attr(g)
  
  # Random graph for small world parameter
  randg <- sim.rand.graph.par(g, N=1)
  sigma <- small.world(g, randg)$sigma
  g$sigma <- sigma

  # Rich club calculation
  l <- rich_club_coeff(g)
  g$phi <- l[[1]]
  g$richclub_params <- l
  return(g)
}


get_graph_params <- function(graphlist, spi=1) {
  # Extract graph parameters of interest
  #
  # Input needs to be a list of lists of lists with
  # one element for each subject at a specific sparsity
  # threshold
  #
  # By default, we extract:
  # - Global efficiency
  # - Local efficiency
  # - Clustering coefficient
  # - Betweenness
  # - Characteristic path length
  # - Vulnerability
  df <- data.frame(ind=c(1:length(graphlist)), Cp=NA, Lp=NA,
                   E.global=NA, E.local=NA, diameter=NA,
                   transitivity=NA, num.hubs=NA, density=NA,
                   assortativity=NA, vulnerability=NA, smallworld=NA,
                   richclub=NA)


  cat("Computing graph metrics...")
  # Use apply instead of loop. Result is list of lists 
  gl <- lapply(graphlist,
               function(x){
                 set_brainGraph_attr(x[[spi]])
               })


  df$Cp            <- unlist(lapply(gl, function(x) x$Cp))
  df$Lp            <- unlist(lapply(gl, function(x) x$Lp))
  df$E.global      <- unlist(lapply(gl, function(x) x$E.global))
  df$E.local       <- unlist(lapply(gl, function(x) x$E.local))
  df$diameter      <- unlist(lapply(gl, function(x) x$diameter))
  df$num.hubs      <- unlist(lapply(gl, function(x) x$num.hubs))
  df$density       <- unlist(lapply(gl, function(x) x$density))
  df$transitivity  <- unlist(lapply(gl, function(x) x$transitivity))
  df$assortativity <- unlist(lapply(gl, function(x) x$assortativity))
  df$vulnerability <- unlist(lapply(gl, function(x) x$vulnerability))

  df$smallworld    <- unlist(lapply(gl,
                                    function(x) {
                                      x <- get_smallworld(x)
                                      x$sigma
                                    }))
  df$richclub     <- unlist(lapply(gl,
                                    function(x) {
                                      x <- get_richclub(x)
                                      x$phi
                                    }))
  df$ind           <- NULL
  cat("done\n")
  return(df)
}

get_smallworld <- function(g)  {
  randg <- sim.rand.graph.par(g, N=1)
  sigma <- small.world(g, randg)$sigma
  g$sigma <- sigma
  return(g)
}

get_richclub <- function(g) {
  l <- rich_club_coeff(g)
  g$phi <- l[[1]]
  g$richclub_params <- l
}
  

    

residualize_df <- function(df, fm, rng) {
  # Regress-out age, sex, intracranial volume
  #print(rng)
  dfr <- df
  for (i in rng) {
    df$y <- df[, i]
    #l1 <- lm(y ~ age + sex + icv, data=tmpdf)
    l <- lm(fm, data=df)
    dfr[, i] <- as.numeric(resid(l))
  }
  return(dfr)
}


plot_components <- function(plsfit, df, comp=c("PLS1", "PLS2")) {
  # Extract loadings
  ld <- scores(plsfit)

  # Plot first two components with slope 
  df <- cbind(df, ld[, c(1:2)]) %>% dplyr::rename(PLS1="Comp 1",
                                                  PLS2="Comp 2") 
  r1 <- (cor.test(df$PLS1, df$slope)$estimate)
  r2 <- (cor.test(df$PLS2, df$slope)$estimate)
  
  dfl <- df %>% gather(key=Component, value=Loading, PLS1, PLS2)

  g <- ggplot(dfl %>% filter(Component %in% comp),
              aes(x=Loading, y=slope)) +
    geom_point(size=8.2) +
    geom_smooth(method="lm", col="black") + 
    facet_wrap(~Component, scale="free") +
    #theme_bw() +
    theme_classic() +
    theme(
      #panel.border = element_blank(),
      #panel.grid.major.y = element_line(colour = "lightgray"),
      #panel.grid.major.x = element_blank(),
      #panel.grid.minor = element_blank(),
      plot.title=element_text(size=30, face="bold", hjust=0),
      axis.line = element_line(colour = "black"),
      axis.text.x=element_text(size=26, angle=NULL),
      axis.text.y=element_text(size=26),
      legend.text=element_text(size=26),
      legend.position="top",
      legend.title=element_blank(),
      strip.text.x=element_text(size=22, face="bold"),
      strip.background = element_blank(),
      #legend.direction="horizontal" ,
      axis.title=element_text(size=30, face="bold")) +
    xlab(label="Score") +
    ylab(label="Slope") +
    #coord_cartesian(ylim=c(0, 80)) + 
    ggtitle("") 
  #scale_fill_manual(values=c("darkblue", "firebrick3"))

  #p <- plot_grid(g1, g2, nrow=2, labels="AUTO", label_size=35, align="h")
  #print(p1)
  return(g)
}


cross_validate <- function(df, k=5, niter=10^3) {
  #
  # Perform k-fold cross_validateidation
  #
  predslope <- vector("numeric", nrow(df))
  obsslope  <- vector("numeric", nrow(df))
  rmse      <- vector("numeric", niter)
  r2        <- vector("numeric", niter)
  cat("Starting cross-validation...")
  for (i in 1:niter) {
    idx      <- sample(1:nrow(df), k, replace=FALSE)
    train    <- df[-idx, ]
    test     <- df[idx, ]
    plsfit   <- plsr(slope ~ ., data=train, validation="CV")
    pred     <- predict(plsfit, ncomp=2, newdata=test)
    pls.eval <- data.frame(obs=test$slope, pred=as.numeric(pred))
    ds       <- defaultSummary(pls.eval)
    #print(ds)
    rmse[i]  <- as.numeric(ds[1]) 
    r2[i]    <- as.numeric(ds[2]) 
  }
  cat("done\n")
  #print(mean(rmse))
  #print(mean(r2, na.rm=TRUE))
  #ds <- defaultSummary(data.frame(obs=train, pred=plsfit))
  return(list(rmse, r2))
}


vals_to_colormap <- function(vals, cols=c("blue", "white")) {
  #
  # Color map stuff for freesurfer
  #
  rr <- range(vals)
  f <- colorRamp(cols)
  svals <- (vals-rr[1])/diff(rr)
  colors <- col2rgb(rgb(f(svals)/255))
  return(colors)
}

z_scored_thickness <- function(df1, df2, vars) {
  #
  # Normalize case (df1) thickness by control (df2) thickness
  #
  zvals <- sapply(vars,
                  function(x) {
                    (df1[, x] - mean(df2[, x]))/sd(df2[, x])
                  })
  return(zvals)
}

compute_slopes <- function(df, fm) {
    #
    # Compute the random slopes
    # NB: Time predictor needs to be last in the formula
    #
    lmerfit <- lmer(formula=fm, data=df)
    # Retrieve index of last predictor which is the time predictor
    ind <- length(fixef(lmerfit))
    slopes <-data.frame(grid=as.numeric(rownames(coef(lmerfit)$grid)),
                        intercept=coef(lmerfit)$grid[, 1],
                        slope=coef(lmerfit)$grid[, ind])
    return(slopes)
}

#zvals <- z_scored_thickness(fefrois, hcdrois, colnames(hcdrois)[1:68])

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_example_individual <- function(df, fn,
                                    visits=c(0, 1, 2, 3, 4,
                                             6, 8, 10, 12)) {
  #
  # Plot illustrative individual
  #

  # Print grid, med, age, sex, dup, iq, bmi
  ddf <- df %>% filter(day==0) %>%
        dplyr::select(grid, age, sex, dup, iq, bmi, med)
  print(ddf)
  
  df$bprsz <- scale(df$bprs)
  df$bprs_tdz <- scale(df$bprs_td)
  dfl <- df %>% gather(key=Scale, value=Score,
                       bprsz, bprs_tdz)
  
 

  cols <- gg_color_hue(2)
  p <- ggplot(dfl %>% filter(visit %in% visits),
              aes(x=visit, y=Score, color=Scale)) +
         geom_point(aes(color=Scale), size=4) +
         geom_line(aes(color=Scale), size=2) +
         theme_gray(base_size=28) +
         theme(
           legend.title=element_blank()
         ) +
         scale_color_manual(labels=c("Positive", "Total"),
                            values=cols) +
         xlab("Week") +
         ylab("Symptoms (normalized)")
  ggsave(plot=p, filename=fn)
}

spaghetti_plot <- function(df, visits=c(0, 1, 2, 3, 4, 6, 8, 10, 12)) {
  #
  # Plot individual and average time courses
  #
  p <- ggplot(df %>% filter(visit %in% visits),
                aes(x=visit, y=bprs_td, group=grid)) +
           geom_point(aes(col=med)) +
           geom_line(aes(col=med)) +
           theme_gray(base_size=28) +
           theme(
             legend.title=element_blank()
           ) +
           xlab("Week") +
           ylab("Positive Symptoms")
  return(p)
}

empty_plot <- function() {
  # Plot an empty plot
  g <- ggplot(data.frame(x=rnorm(10),
                    y=rnorm(10)),
                    aes(x=x, y=y)) +
    geom_point(col="white") +
    theme_void() +
    xlab("") +
    ylab("")
  return(g)
}


 

 
