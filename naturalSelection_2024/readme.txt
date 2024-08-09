Goal: Run natrual selection analysis on genotype data 
Date: 2023
Main: selscan

This directory contains the following:
1. Bash scripts to run sequenially 
step1: run cmd1_runTemplate to run Eagle (phasing, imputation), and selscan 
step2: run cmd2_norm_ihs to generate normalized iHS scroes
2. bash template file to be used by program in step1.
   One for run Eagle, two for run selscan
cmdTemplate_runEagle
cmdTemplate_runSelscan_gmap, using genetic map
cmdTemplate_runSelscan: not used 

Ref:
1. https://academic.oup.com/mbe/article/31/10/2824/1012603?login=true
Zachary A. Szpiech and Ryan D. Hernandez
selscan: An Efficient Multithreaded Program to Perform EHH-Based Scans for Positive Selection Molecular Biology and Evolution 2014, 31 (10):2824–2827 doi:10.1093/molbev/msu211
