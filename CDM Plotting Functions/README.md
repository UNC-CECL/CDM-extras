## Description
This directory contains a series of R functions (`CDM_plotting_functions`) for generating exploratory data analysis figures and animations of Coastal Dune Model (CDM) simulations using ggplot2 and plotly (3D plot and 3D animation).

R library dependencies are: plotly, ggplot2, plyr, gridExtra, and viridis.

#### Example
An example batch CDM simulation series (`Example_Batch_CDM_Simulation`) is also included in which the vegetation line (Lveg) is parameterized at 20, 40, and 60m (`model_iter1`-`model_iter3`). 

Example plots of a CDM simulations: dune height over time for each Lveg parameterization (`h.timeseries.png`), cross-shore profile of dune height at simulation year 26 for each Lveg parameterization (`h.final.png`), and an animation of the entire CDM timeseries (`CDM_animation.html`; Lveg = 40m).

R code to generate the example plots and animations from the included dataset in `Example_CDM_Plots.R`
