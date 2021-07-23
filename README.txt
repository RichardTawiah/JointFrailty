######## Supplementary Material ##############
#Title: Joint frailty modelling of time-to-event data to elicit the evolution pathway of events: A
        generalised linear mixed model approach

The compressed file contains a readme file "README.txt", an example data "data_GitHub" with two covariates 
and an R-command file "JointFrailty.R".

1. data_GitHub  # example data 

2. Description of variables in the data
   hosp       # Hospital
   id         # Patient's identification number 
   time       # Gap time for recurrent events or death (death time is the last available time for each patient)
   delta1     # Censoring indicator for recurrent events 
   delta2     # Censoring indicator for death
   age        # Covariate
   male       # Covariate (female=0, male=1)

3. JointFrailty.R  # R code for the proposed joint frailty model
   
   # Notes:
   # R for Windows Version 4.1.0 was used in the illustration 
   # Install "dplyr" package through install.packages("dplyr")
   # Library "dplyr" package 
   # Set work directory: setwd("C:/Users/...")
   # Load the data into R: see the preamble of the code
   # theta01, theta02, rho0 are initial values for thetau^2, thetav^2, and rho, respectively. 
   # The main function is "joint.frailty"; it calls the output of the model 
   # (with initial value set to 0.2 for all variance components and maximum iteration = 300)
   # The output shows the progress of the estimation procedure till the final results.  
   




   

 
