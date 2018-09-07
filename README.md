
# Structural similarity networks and treatment response



# Author

Philipp Homan <phoman1 at northwell dot edu>


# Getting Started

This repository contains all the data and analysis code to reproduce the
manuscript " Structural similarity networks predict clinical outcome in
early-phase schizophrenia". These instructions describe how to obtain a
copy of the project up and running on your local machine for reproducing
the analysis described in the manuscript. The repository contains a
Makefile which reflects the dependencies of the analysis; analysis,
figures and manuscript can be produced by simply typing 'make' from the
Unix command line.


## Prerequisites

All analyses were conducted with the R software 
R version 3.3.2 (2016-10-31).  Mixed models were
estimated using the lme4 library, partial least squares regression were
computed with the pls library, and brain graph metrics with the
brainGraph and igraph libraries.  
Python 2.7.13 and 
pysurfer (0.7) were used for visualizing the imaging
results. The full session info under R can be found at the end of this
file


# Installing

Clone the repository or download the zip file and run 'unzip ssn.zip'
from the command line.


# Running the analysis

Change to the ssn directory and run 'make analysis'.


# Producing the figures

Change to the ssn directory and run 'make figures'. The figures can then
be found in output/figures.


# Producing the manuscript

Change to the ssn directory and run 'make manuscript'. The manuscript
will be in src/fe<sub>ms.pdf</sub>


# Built With

Org-mode 9.1.9.


# Session info

    R version 3.3.2 (2016-10-31)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 17.04
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    attached base packages:
    [1] grid      stats     graphics  grDevices utils     datasets  methods  
    [8] base     
    
    other attached packages:
     [1] pls_2.6-0         caret_6.0-78      MASS_7.3-47       dplyr_0.7.4      
     [5] tidyr_0.8.0       brainGraph_2.0.0  cairoDevice_2.24  RGtk2_2.20.31    
     [9] igraph_1.1.2      astsa_1.8         corrplot_0.84     tibble_1.4.2     
    [13] Hmisc_4.0-3       Formula_1.2-2     survival_2.40-1   lsmeans_2.27-61  
    [17] lme4_1.1-14       Matrix_1.2-8      boot_1.3-18       DescTools_0.99.23
    [21] flexmix_2.3-14    lattice_0.20-34   png_0.1-7         magick_1.7       
    [25] car_2.1-6         ellipse_0.4.1     mixtools_1.1.0    cowplot_0.9.2    
    [29] ggplot2_2.2.1     lm.beta_1.5-1     R.matlab_3.6.1    pacman_0.4.6     
    
    loaded via a namespace (and not attached):
      [1] TH.data_1.0-8       minqa_1.2.4         colorspace_1.3-2   
      [4] RcppEigen_0.3.3.3.1 class_7.3-14        modeltools_0.2-21  
      [7] estimability_1.3    htmlTable_1.9       base64enc_0.1-3    
     [10] DRR_0.0.3           MatrixModels_0.4-1  ggrepel_0.7.0      
     [13] lubridate_1.7.1     prodlim_1.6.1       manipulate_1.0.1   
     [16] mvtnorm_1.0-6       codetools_0.2-15    splines_3.3.2      
     [19] R.methodsS3_1.7.1   mnormt_1.5-5        robustbase_0.92-8  
     [22] knitr_1.20          RcppRoll_0.2.2      ade4_1.7-10        
     [25] nloptr_1.0.4        pbkrtest_0.4-7      broom_0.4.5        
     [28] ddalpha_1.3.1.1     kernlab_0.9-25      cluster_2.0.5      
     [31] sfsmisc_1.1-1       R.oo_1.21.0         backports_1.1.1    
     [34] assertthat_0.2.0    lazyeval_0.2.1      acepack_1.4.1      
     [37] htmltools_0.3.6     quantreg_5.34       tools_3.3.2        
     [40] bindrcpp_0.2        coda_0.19-1         gtable_0.2.0       
     [43] glue_1.2.0          reshape2_1.4.2      Rcpp_0.12.13       
     [46] nlme_3.1-131        iterators_1.0.9     psych_1.7.8        
     [49] timeDate_3042.101   gower_0.1.2         stringr_1.3.1      
     [52] lpSolve_5.6.13      DEoptimR_1.0-8      zoo_1.8-0          
     [55] scales_0.5.0        ipred_0.9-6         parallel_3.3.2     
     [58] sandwich_2.4-0      expm_0.999-2        SparseM_1.77       
     [61] RColorBrewer_1.1-2  gridExtra_2.3       rpart_4.1-10       
     [64] segmented_0.5-3.0   latticeExtra_0.6-28 stringi_1.2.3      
     [67] foreach_1.4.4       RNifti_0.7.1        checkmate_1.8.5    
     [70] permute_0.9-4       lava_1.6            oro.nifti_0.9.1    
     [73] rlang_0.2.0         pkgconfig_2.0.1     bitops_1.0-6       
     [76] purrr_0.2.4         bindr_0.1           recipes_0.1.2      
     [79] htmlwidgets_0.9     tidyselect_0.2.3    CVST_0.2-1         
     [82] plyr_1.8.4          magrittr_1.5        R6_2.2.2           
     [85] dimRed_0.1.0        multcomp_1.4-8      withr_2.1.0        
     [88] pillar_1.2.2        foreign_0.8-67      mgcv_1.8-16        
     [91] abind_1.4-5         nnet_7.3-12         data.table_1.10.4-3
     [94] ModelMetrics_1.1.0  digest_0.6.12       xtable_1.8-2       
     [97] mediation_4.4.6     R.utils_2.6.0       stats4_3.3.2       
    [100] munsell_0.4.3      

