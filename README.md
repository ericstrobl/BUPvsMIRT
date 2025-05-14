# Bupropion vs. Mirtazapine
This repository contains code needed to replicate the experimental results in the paper entitled ``Consistent Differential Effects of Bupropion and Mirtazapine in Major Depression.'' All code was tested in R version 4.3.1.

The STAR*D and CO-MED datasets can be downloaded from the [NIMH Data Archive](https://nda.nih.gov/) with a limited access data use certificate. 

# Installation
First install the BiocManager and qvalue packages. Then:

> library(devtools)

> install_github("ericstrobl/BUPvsMIRT")

> library(BUPvsMIRT)

# Analysis Pipeline

Since the STAR*D and CO-MED datasets requires prior approval, we illustrate the pipeline with synthetic data.

> data = generate_synth(nsamps=1000, nF=3) # generate synthetic data with 1000 samples, five treatments and three latent factors

> Y1 = lm.fit(cbind(data$X,1),data$Y)$residuals # partial out nuisance variables

> nc=num_components(data$Tx,Y1,ncs=2:5) # determine number of components

> mod=SV(Tx=data$Tx,Y=Y1,nc=nc) # run the SV algorithm with nc components

> res=SInf_permutation_test(mod,data$Tx,Y1,nc,ncs=2:5) # perform permutation testing with post-model selection inference

> print(res)

# Outputs

A list with:

`omnibus` = omnibus statistic and p-values

`factors` = uncorrected and corrected p-values for each factor

`tx_pairs` = uncorrected and corrected p-values for each treatment pair





