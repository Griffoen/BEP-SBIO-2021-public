# BEP-SBIO-2021
In this project, a kinetic model of a generalized extracellular molecule sensor (GEMS) platform that activates the Jak-STAT3 pathway is developed and analyzed. The goal is to explore receptor activation and signal transduction in the GEMS system. This gives the opportunity to explore optimal design parameters which can be implemented experimentally. 


## Table of Contents
* [General Info](#general-information)
* [Technologies Used](#technologies-used)
* [Features](#features)
* [Usage](#usage)
* [Project Status](#project-status)
* [Acknowledgements](#acknowledgements)


## General Information
BEP-SBIO-2021 contains a model that describes a synthetic biological system of GEMS system and the Jak-STAT3 pathway. The model consist of 22 ordinary differental equations (ODEs), that are solved with the ode15s solver.
The purpose of this project is to investigate the behavior of the system (the receptor activation and the signal transduction) and to explore the effect of changing some constants.


## Technologies Used
- Matlab - version 2018b


## Features
The following three scripts simulate the model:
- ODEs.m                    - Contains the constants and the ODEs.
- Defaults.m                - The default settings.
- OdeSolver.m               - Solves the ODEs and plot the results.

The following scripts can be used to investigate the system: 
- TEST_ODE.m			          - To execute some simple tests to check whether the model satisfy predefined boundary conditions.
- Plot_BellCurve_LvsR.m		  - To create a bell-shaped dose response curve.
- Plot_LigandvsReceptor.m	  - To plot the effect of variations of the ligand-receptor ratio on the signal transduction.
- Plot_SEAP_Variations.m	  - To plot the effect of variations of some constants on the SEAP expression.
- CreatePlots.m             - To execute the above described functions
- Vary_Dimerization_KD.m	  - The effect of variations in the Kd value on the required time for STAT3npd to stabalize


## Setup
To solve the model, the following files are required: ODEs.m, Defaults.m, and OdeSolver.m. These codes are to solve the system.
There are several other files that contain one or more functions to analyse the system. In CreatePlots, all those functions can be executed. 


## Usage
To simulate the model and to plot the results, run the file OdeSolver.m.
To execute the simple tests or to variate parameters, the file CreatePlots.m should be used. When running this code, you can choose which test or function you want to run by entering the number of that test.


## Project Status
Project is complete.


## Acknowledgements
- This project was inspired by the work of prof. Tom de Greef.
- Many thanks to Tom de Greef, Glenn, Alex, and Anna.

