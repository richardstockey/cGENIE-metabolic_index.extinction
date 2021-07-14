################################################################
### readme.txt #################################################
################################################################

For: 'Decreasing Phanerozoic extinction intensity as a consequence of Earth surface oxygenation and metazoan ecophysiology'

Richard G. Stockey1*, Alexandre Pohl2,3, Andy Ridgwell2, Seth Finnegan4, and Erik A. Sperling1

1 Department of Geological Sciences, Stanford University, Stanford, CA 94305, USA
2 Department of Earth and Planetary Sciences, University of California, Riverside, CA, USA
3 Biogéosciences, UMR 6282, UBFC/CNRS, Université Bourgogne Franche-Comté, 6 boulevard Gabriel, F-21000 Dijon, France
4 Department of Integrative Biology, University of California, Berkeley, Berkeley, CA, USA
* Corresponding author: rstockey@stanford.edu

################################################################

This repository contains the required R code to run the ecophysiological model analyses presented in this study, and generate the associated figures.

For cGENIE model code - see https://github.com/derpycode/cgenie.muffin/genie-userconfigs/MS/stockeyetal.PNAS.2021

This repository contains two directories:
> Viability extractions [contains R code for all ecophysiological model analyses presented in this study]
> Figures [contains R code to generate all figures]

To replicate the analyses and plots presented here, the following R packages are required:
AquaEnv
deeptime
dplyr
ggplot2
ncdf4
pals
reshape2
RColorBrewer
RNetCDF

The deeptime package is currently only available on GitHub, and will need to be installed from there to reproduce all multipanel figures and plots (due to their dependency on the function ggarrange2). This can be achieved by running the following commands in your R console (ignore the first line if you already have devtools installed).
install.packages("devtools")  
devtools::install_github("willgearty/deeptime") 

All other packages can be installed from CRAN. These scripts have been tested using R version 3.6.2 - Copyright (C) 2019 The R Foundation for Statistical Computing.

Note that the R script "Hofmannn functions.R" is also required for ecophysiological model analyses. If not using the directory structure described above this will need to be manually sourced. 

These R scripts assume that the cGENIE ensembles presented in this study are stored in separate directories as follows:
> fm0450.CO2.O2 [directory for both modern and shallow remineralisation depth Ordovician experiments]
> p0251.CO2.O2.ramp
> p0055.CO2.O2.ramp 

These directory names generally match those in https://github.com/derpycode/cgenie.muffin/genie-userconfigs/MS/stockeyetal.PNAS.2021 - however note that the outputs for user config directories stockeyetal.PNAS.2021/fm0450.CO2.O2 and stockeyetal.PNAS.2021/fm0450.CO2.O2.shallow.remin are here stored in the same output directory.
These scripts further assume that the directories described above are stored in a parent directory named "GENIE.output". For full replication of these methods - the full recommended directory structure for storing cGENIE outputs and ecophysiological modeling results is as follows:

> GENIE.output
	> fm0450.CO2.O2
	> p0251.CO2.O2.ramp
	> p0055.CO2.O2.ramp 
> Habitat.viability.summaries
> Habitat.viability.netcdfs
	> fm0450.CO2.O2
	> p0251.CO2.O2.ramp
	> p0055.CO2.O2.ramp 

For the scripts in this repository, it is assumed that the whole "cGENIE-metabolic_index.extinction" directory has been downloaded and this directory structure was added within it (such that five sub-directories are stored within the primary directory).

Viability extractions that should be run before each figure can be generated:
Figure 4 	- 	Viability_extraction_fm0450.CO2.O2_1EFD_shelves.no.poles.R
Figure 5 	- 	Viability_extraction_p0055.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles.R
			Viability_extraction_p0251.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles.R
			Viability_extraction_fm0450.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles.R
Figure S4 	- 	Viability_extraction_p0055.CO2.O2_1EFD_shelves.no.poles.R
			Viability_extraction_p0251.CO2.O2_1EFD_shelves.no.poles.R
			Viability_extraction_fm0450.CO2.O2_0.34EFD_shelves.no.poles.R
Figure S5 	- 	Viability_extraction_fm0450.CO2.O2_1EFD_global.no.poles.R
Figure S6 	- 	Viability_extraction_p0055.CO2.O2_1EFD_global.no.poles.R
			Viability_extraction_p0251.CO2.O2_1EFD_global.no.poles.R
			Viability_extraction_fm0450.CO2.O2_0.34EFD_global.no.poles.R	
Figure S7 	- 	Viability_extraction_fm0450.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles.R	
			Viability_extraction_fm0450.CO2.O2_0.34EFD_shelves.ecotype.sensitivity.no.poles.R	
Figs S8-10, S12 - 	Viability_extraction_p0055.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles.R
			Viability_extraction_p0251.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles.R
			Viability_extraction_fm0450.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles.R
Figure S11	- 	Viability_extraction_p0055.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_lowA0_lowEo.R
			Viability_extraction_p0055.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_highA0_lowEo.R
			Viability_extraction_p0055.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_lowA0_highEo.R
			Viability_extraction_p0055.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_highA0_highEo.R
			Viability_extraction_p0251.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_lowA0_lowEo.R
			Viability_extraction_p0251.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_highA0_lowEo.R
			Viability_extraction_p0251.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_lowA0_highEo.R
			Viability_extraction_p0251.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_highA0_highEo.R
			Viability_extraction_fm0450.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_lowA0_lowEo.R
			Viability_extraction_fm0450.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_highA0_lowEo.R
			Viability_extraction_fm0450.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_lowA0_highEo.R
			Viability_extraction_fm0450.CO2.O2_1EFD_shelves.ecotype.sensitivity.no.poles_highA0_highEo.R
			
Quantitative paleobiological analyses, along with code to produce and save associated figures and tables, can all be found in the script "Extinction_rates_ectotherms.R".

################################################################
################################################################
################################################################
