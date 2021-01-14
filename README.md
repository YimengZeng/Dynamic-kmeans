# Dynamic-kmeans
## Checklist for Data Analysis Conducted in the manuscript
  Lead author(s): Yimeng Zeng
  
  Created by Yimeng Zeng on September, 2020
  
  Updated by Yimeng Zeng on January 14th, 2021
## Summary and Requirements:
  This repository contains Matlab, Python scriptes needed to produce key results in the manuscript including figures, statstic test, fmri data preprocess and machine learning procedure.
## Test Environment(s):
MATLAB version: R2015b,R2020a  Python Version: 3.6  Operating System: windows 10, linux centos6_x86_64
## SCL raw data
  * preprocsss procedure for SCLs (conducted only on MATLAB R2020a)
  
    script: scr_preprocess.m

## Fmri data analysis

### Seed-based whole-brain functional connectivity from resting-state fMRI （Figure 1）
  * whole-brain fc analysis
  
    script: fconnect_wm_csf_nogs_final.m,fconnect_wm_csf_nogs_final_config.m
  
    
### Dynamic functional connectivity analysis with sliding window and k-means clustering (Figure 1)
  * sliding window analysis
  * k-means clustering
  * cluster result evalution
  
    script: sliding_window_and_k_means_clustering.m (line 1---line 69)
  
    Figures in the MS: Figure.1A, plotted by using MRIcron and MRIcroGL；Figure.1B, plotted by workbench; Figure.1C, plotted by python. Others plotted by MATLAB
### State based analysis of their temporal and spatial properties(Figure 2)
  * correlation analysis 
  * time-lagged cross correlation analysis
  
    script: sliding_window_and_k_means_clustering.m (line 70---line 209)
  
    Figures in the MS: Figure.2A plotted by python; Figure.2B,2C,2D,2E plotted by MATLAB
### Elastic-net regression for predicting skin conductance level (Figure 3,Figure 4)
  * elastic-net regression analysis
  
    script: state_based_analysis.m
    
    Figures in the MS: Figure.3A,3B plotted by MATLAB; Figure.3C plotted by Python
    
    Figure.4A,4C plotted by MATLAB and Excel; Figure.4B,4D plotted by BrainNet viewer (Xia et al., 2013)
