# Homeostatic_circadian_light_model-factors_driving_sleep_phenotypes

This repository contains data and code used to produce the results in the paper:
"Method to determine whether sleep phenotypes are driven by endogenous circadian
rhythms or environmental light by combining longitudinal data and personalised mathematical models"
Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.

Please direct all queries to: a.skeldon@surrey.ac.uk, University of Surrey, 2023.

There are four sub-folders:

data_raw:
data collected from participants and used in the analyses. These include
(i) daily sleep diary data for each participant
(ii) light data for the study period (one minute resolution)
(iii) mean slow wave activity data for the laboratory night for each participant

data_generated:
data generated by fitting and simulation of the HCL model. These include
(i) data generated from simulating the model for different parameter values - used in
the sensitivity analysis
(ii) data generated to construct the Arnold tongue entrainment regions
(iii) fitted parameter values
(iv) light metrics
(v) simulated timeseries for each participant
(vi) value of the objective function for a grid of points

utilities:
matlab functions used to load in data, calculate light metrics, specify
the equations, integrate the equations, evaluate objective functions, carry out
the fitting, some plotting functionality

main_code:
master routines to run the different analyses / fitting procedures. Calls
functions in utilities and uses data in data_raw.
~                                                                   
