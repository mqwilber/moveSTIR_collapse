# moveSTIR_collapse

This repository contains code, simulations and documents relating to a generalization of the moveSTIR framework for linking spatial movement data with epidemiological transmission dynamics. The original framework was developed by Mark Wilber and published in Ecology Letters in 2022. The goal is to make the framework more general and to provide equations to analyze different types of data, or with different assumptions relating to specific case study scenarios. The framework allows to estimate a force of infection (FOI) based on movement and spatial use data. 
The framework is called Probabilistic MoveSTIR (PMoveSTIR) and the manuscript will be published soon. The main steps in the workflow are 1) the calculation of a utilization distribution based on detailed movement data, 2) the calculation of the local correlation in space use, and 3) the combination of both of these into a single calculation, scaled by epidemiological processes like parasite shedding decay rates. We included a simulation study where we explore how some factors like the contact distance, decay rate, and attraction among hosts affects the expected FOI. We also included an empirical analysis where we estimate the FOI for five deer using the PMoveSTIR framework.

The data folder contains a csv file with GPS tracking data for multiple individuals of white-tailed deer. These data were used for the empirical example in the manuscript
The docs folder contains LaTeX files with the main manuscript and supplementary materials. It also contains a figures subfolder with figures included in the manuscript.
The code folder contains all the files used to generate the results for the simulation and empirical components of the paper. The main code files are:
- `functions.R` and `deer_functions.R`: These are R scripts that define the functions used to analyze the simulated movement data and empirical data, respectively, to estimate the FOI
- `deer_analysis.R` is the main code for the empirical analysis. It is a script that performs all the calculations and some visualizations for the deer data
- `run_sims_cluster.R` is the R script used to run a wide variety of different scenarios for the simulation study in a HPC cluster.
- `single_sim_example.R` is a basic script that generates movement data, estimates the utilization distribution, the correlation, and FOI.
