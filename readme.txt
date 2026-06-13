Readme for replication files for "Spillovers from US monetary policy: Evidence from a time-varying parameter GVAR model" forthcoming in the Journal of the Royal Statistical Society: Series A

- You need to install the threshtvp package from source. Windows and Mac versions are included in the zip file

- The script file necessary to replicate the main results is 01b_baseline_estimation_sp_levels.r . The remaining files include either the data (i.e. the ones ending with .rda) or additional helper files (ending with .r) that are necessary to estimate the TVP-GVAR model

- Before estimating the model, make sure that all files in the zip are unpacked into a single folder and this folder is set as your working directory in R (you can use the setwd("  ") command).

- You need to adjust the number of draws (lines 30-32). If you use this script on a standard personal computer, make sure that you change the number of CPU cores used by your computer (see line 71) to be equal to the number of available CPU cores (i.e. CPU = 4 if you have a quad core processor)

- After doing so, simply source the 01b_baseline_estimation_sp_levels.r file. This will start the main estimation routine. Notice that, due to the complexity of the model and the computations involved, this could take a few hours. 

- The first part estimates the N+1 country models. The second part of the code then proceeds by computing impulse responses identified through standard sign restrictions

- All results are then stored in a separate folder Results. These include impulse responses as well as plots that show the normalized log-determinant of the process innovation variances over time

- In case of questions, contact Florian.Huber@sbg.ac.at. Website: https://sites.google.com/site/fhuber7/home 


