# Bupropion vs. Mirtazapine
This repository contains code needed to replicate the experimental results in the paper entitled ``Consistent Differential Effects of Bupropion and Mirtazapine in Major Depression.'' All code was tested in R version 4.3.1.

The STAR*D and CO-MED datasets can be downloaded from the [NIMH Data Archive](https://nda.nih.gov/) with a limited access data use certificate. 

# Installation
First install the BiocManager and qvalue packages. Then:

> library(devtools)

> install_github("ericstrobl/BUPvsMIRT")

> library(BUPvsMIRT)

# Analysis Pipeline

We illustrate the pipeline with synthetic data, since the STAR*D and CO-MED datasets require prior approval:

> data = generate_synth() # generate synthetic data with 1000 samples, five treatments and three latent factors

> Y1 = lm.fit(cbind(data$X,1),data$Y)$residuals # partial out nuisance variables

> nc=num_components(data$Tx,Y1,ncs=2:5) # determine number of components

> mod=SV(Tx=data$Tx,Y=Y1,nc=nc) # run the SV algorithm with nc components

> res=permutation_test_SelInf(mod,data$Tx,Y1,nc,ncs=2:5) # perform permutation testing with post-model selection inference

> print(res)

> res_rep=permutation_test_SelInf(mod,data$Tx,Y1,nc,ncs=2:5,Tx_focus=c(1,2)) # perform permutation testing with post-model selection inference but comparing two particular treatments in replication testing (e.g., with CO-MED)

# Permutation Testing Outputs

A list with:

`omnibus` = omnibus absolute sum statistic and p-value

`factors` = absolute sum statistic, uncorrected p-value, and corrected p-value for each factor

`tx_pairs` = difference statistic, uncorrected p-value, and corrected p-value for each treatment pair





