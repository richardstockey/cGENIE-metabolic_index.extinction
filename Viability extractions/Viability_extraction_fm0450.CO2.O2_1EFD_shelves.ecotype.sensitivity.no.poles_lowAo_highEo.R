library(AquaEnv)
library(dplyr)
library(ncdf4)
library(reshape2)
library(RNetCDF)

################################
## Aim: to read in cGENIE output files for the ensembles of Earth system modeling experiments used in Stockey et al. 
## "Decreasing Phanerozoic extinction intensity as a consequence of Earth surface oxygenation and metazoan ecophysiology"
## Generate datasets (± NetCDFs) describing how metabolic habitat viability responds to this...
################################

################################
## OPTIONS!
################################

# Set working directory that cGENIE outputs are stored in
genie.directory <- "/cGENIE-metabolic_index.extinction/GENIE.output/"

# Set save directory for RData file summaries
rdata.directory <- "/cGENIE-metabolic_index.extinction/Habitat.viability.summaries/"

# Set save directory for NetCDFs (ignore this if not creating NetCDFs, as prescribed by options below)
netcdf.directory <- "/cGENIE-metabolic_index.extinction/Habitat.viability.netcdfs/"

# Source functions from Hofmann et al. 2011 Deep Sea Research (set directory if necessary)
source("Hofmann functions.R")

# Prescribe vectors of O2 forcings of interest (CO2 forcings prescribed within tectonic config!)
O2.vec <- c("0.05", "0.1", "0.2", "0.4", "0.6", "0.8", "1")

# Set e-folding depth
EFD <- 1

# Set tectonic configuration (this changes directory and folder names) - if not one of those used in oxygen-extinction paper then define path (see examples below) as tectonics parameter
tectonics <- "Ordovician"

# Did you use a ramp for this configuration? Yes = ramp <- 1; No = ramp <- 0
ramp <- 0

if(ramp == 1){
  ramp.text <- ".ramp"
}else{
  ramp.text <- ""
}


# Generate a NetCDF of Metabolic Habitat Viability for each model?
# Use habitat.viability.netcdf == 0 as default for speed/storage benefits. habitat.viability.netcdf == 1 will generate NetCDFs in your chosen directory.
habitat.viability.netcdf <- 0

# Surface layer only?  (only activate shelves OR surface - by setting to == 1 - at any given time)
surface.only <- 0

# Shelves only? (only activate shelves OR surface - by setting to == 1 - at any given time)
shelves.only <- 1

# Set file naming conventions based on *options*
if(surface.only == 1){
  ocean.mask <- "surface"
}else if(shelves.only == 1){
  ocean.mask <- "shelves"
} else{
  ocean.mask <- "global" # global is not a read-in variable below, just for record keeping
}

if(tectonics == "Ordovician"){
  tectonic.config <- "fm0450.CO2.O2/muffin.CB.fm0450dc.BASES."
  config.short <- "fm0450"
  CO2.vec <- c(4,8,16,32,64,128) 
}else if(tectonics == "Permian"){
  tectonic.config <- "p0251.CO2.O2.ramp/cgenie.eb_go_gs_ac_bg.p0251b.BASES."
  config.short <- "p0251"
  CO2.vec <- c(2,4,8,16,32,64) 
}else if(tectonics == "Paleocene"){
  tectonic.config <- "p0055.CO2.O2.ramp/cgenie.eb_go_gs_ac_bg.p0055c.BASES."
  config.short <- "p0055"
  CO2.vec <- c(0.5,1,2,4,8,16) 
}else{
  tectonic.config <- tectonics
  config.short  <- tectonics # if using, you might want to update this for filing purposes
  CO2.vec <- c(4,8,16,32,64,128)  # Change as appropriate if using this script for a configuration not used in this study
}

################################
## Loop through 100 different ecotype populations!
################################

## data frame initiation moved up!
# Initiate vectors and dataframes
ecotypes.summary <-  data.frame(ecotype=double(), cells = double(), O2.pal=double(), CO2.pal=double(), iteration=double())
ecotypes.summary.short <-  data.frame(ecotypes.number = double(), O2.pal=double(), CO2.pal=double(), eq.temp=double(), iteration=double())

# get a vector of set.seed numbers
seed.vec <- runif(100, min=1, max=5000)

# start loop
for(it in 1:100){
  
################################
## Define parameters for metabolic index model
################################

# Set number of ecophysiotypes in model - abbreviated to 'ecotype' 
# throughout this script for concision
init.ecotypes <- 1000

# set seed for random number generation
set.seed(seed.vec[it])

# input Penn et al. 2018 parameter values
u_A <- 3.01 # log(atm ^-1)
sig_A <- 0.49 # log(atm ^-1)
u_A_low <- log(qlnorm(0.25, meanlog = u_A, sdlog = sig_A)) # find 25th percentile and then retransform to log space
u_A_high <- log(qlnorm(0.75, meanlog = u_A, sdlog = sig_A)) # find 75th percentile and then retransform to log space
sig_A_half <- log(sqrt((exp(sig_A)^2)/2)) # here we half variance, recompute sd and then retransform to log space

u_C <- 1.10 # unitless
sig_C <- 0.42 # unitless

u_E <- 0.41 # eV
sig_E <- 0.29 # eV
u_E_low <- qnorm(0.25, mean = u_E, sd = sig_E) # find 25th percentile
u_E_high <- qnorm(0.75, mean = u_E, sd = sig_E) # find 75th percentile
sig_E_half <- sqrt((sig_E^2)/2) # here we half variance and recompute sd

# sample Ao values from PDF (log-normal distribution)
A0.xxx <- rlnorm(init.ecotypes, meanlog = u_A_low, sdlog = sig_A_half)

# sample Eo values from PDF (normal distribution)
E0.xxx <- rnorm(init.ecotypes, mean = u_E_high, sd = sig_E_half)

# sample phi crit values from PDF (log-normal distribution)
phi_crit.xxx <- rlnorm(init.ecotypes, meanlog = u_C, sdlog = sig_C)

# Set Penn et al. 2018 constants
kB <- 8.61733326e-5
Tref <- 15+273.15

print(system.time(
for (O2.factor in O2.vec){
    for (CO2.factor in CO2.vec){

    Nc <- open.nc(paste0(genie.directory, tectonic.config, CO2.factor, "CO2.", O2.factor, "O2.", EFD, "EFD", ramp.text, ".config/biogem/fields_biogem_3d", ".nc"))

    ## =============================================================
    ## Extract annual means of key variables
    ## Only need to extract lat-lon from one file because all NetCDF files use same 1 degree bins
    ## Optional surface plots below value extraction
    ## =============================================================
    lat <- var.get.nc(Nc, "lat") # units: degrees north
    lat.edges <- var.get.nc(Nc, "lat_edges")
    lon <- var.get.nc(Nc, "lon") # units: degrees east
    lon.edges <- var.get.nc(Nc, "lon_edges")
    depth <- var.get.nc(Nc, "zt") # units: metres
    depth.edges <- var.get.nc(Nc, "zt_edges") # units: metres
    time <- var.get.nc(Nc, "time") # units: year mid-point
    oxy <- var.get.nc(Nc, "ocn_O2") # units: mol kg-1
    sal <- var.get.nc(Nc, "ocn_sal") # units: PSU
    temp <- var.get.nc(Nc, "ocn_temp") # units: degrees C
    topo <- var.get.nc(Nc, "grid_topo") # units: meters
    
    # =============================================================
    # No density output from cGENIE - so here we use seadensity(S,t) function from Hofmann et al 2011
    # Seadensity takes S in PSU and t in degrees C so no change from cGENIE required
    # =============================================================

    ## =============================================================
    ## set desired time bin (last of model simulation)
    ## =============================================================

    time <- as.numeric(length(time))

    ## =============================================================
    ## set desired dimensions of array
    ## =============================================================

    xdim <- as.numeric(length(lon))
    ydim <- as.numeric(length(lat))
    if(surface.only == 1){
    # If surface only - set maximum z dimension to 1  
    zdim <- 1
    }else if(shelves.only == 1){
    # If shelves only - set maximum z dimension to 3
    zdim <- 3

    ## =============================================================
    ##  save mean equatorial sea surface temperature for plotting 
    ## =============================================================
    eq.temp <- mean(temp[,18:19,1,time], na.rm=T)
    
    ## =============================================================
    ## if only evaluating shelf settings, apply mask to cut out shelves only
    ## these masks are defined as manual 'cut-outs', setting all ocean cells not adjacent to
    ## the continents as NAs. Only surface 3 ocean layers are used for shelves.
    ## =============================================================

    #### Ordovician shelves cutout
    if(grepl("fm0450dc", tectonic.config, fixed = TRUE) == TRUE ){
      value <- NA

      temp[11:23,1,1,time] <- value
      temp[12:22,2,1,time] <- value
      temp[c(13:23,26:27),3,1,time] <- value
      temp[12:29,4,1,time] <- value
      temp[c(11:16,18:30),5,1,time] <- value
      temp[c(10:14,19:31),6,1,time] <- value
      temp[c(9:13,21:31),7,1,time] <- value
      temp[c(7:14,16,22:29,31),8,1,time] <- value
      temp[c(7:16,22:28,32),9,1,time] <- value
      temp[c(3,7:13,16,22:28,32),10,1,time] <- value
      temp[c(3,7:12, 22:28),11,1,time] <- value
      temp[c(3,7:12, 21:29),12,1,time] <- value
      temp[c(6:10,17,21:29,33),13,1,time] <- value
      temp[c(3,5:9,17:18,20:29,33),14,1,time] <- value
      temp[c(3:9,18:29),15,1,time] <- value
      temp[c(3:10,12:13,18:29,33),16,1,time] <- value
      temp[c(2:12,18:29,33),17,1,time] <- value
      temp[c(2:11,18:29),18,1,time] <- value
      temp[c(2:11,18:30,32),19,1,time] <- value
      temp[c(1:12,18:19,21:32),20,1,time] <- value
      temp[c(1:15,17:18,22:33,36),21,1,time] <- value
      temp[c(1:17,22:32,36),22,1,time] <- value
      temp[c(1:17,22:32,36),23,1,time] <- value
      temp[c(1:17,22:31,36),24,1,time] <- value
      temp[c(1:15,17,23:30,35:36),25,1,time] <- value
      temp[c(1:14,23:31,33:36),26,1,time] <- value
      temp[c(1:15,17,23:36),27,1,time] <- value
      temp[c(1:18,23:36),28,1,time] <- value
      temp[c(1:18,23:36),29,1,time] <- value
      temp[c(1:18,22:36),30,1,time] <- value
      temp[c(1:19,21:36),31,1,time] <- value
      temp[c(1:36),32,1,time] <- value
      temp[c(1:36),33,1,time] <- value
      temp[c(1:36),34,1,time] <- value
      temp[c(1:36),35,1,time] <- value
      temp[c(1:36),36,1,time] <- value

      temp[11:23,1,2,time] <- value
      temp[12:22,2,2,time] <- value
      temp[c(13:23,26:27),3,2,time] <- value
      temp[12:29,4,2,time] <- value
      temp[c(11:16,18:30),5,2,time] <- value
      temp[c(10:14,19:31),6,2,time] <- value
      temp[c(9:13,21:31),7,2,time] <- value
      temp[c(7:14,16,22:29,31),8,2,time] <- value
      temp[c(7:16,22:28,32),9,2,time] <- value
      temp[c(3,7:13,16,22:28,32),10,2,time] <- value
      temp[c(3,7:12, 22:28),11,2,time] <- value
      temp[c(3,7:12, 21:29),12,2,time] <- value
      temp[c(6:10,17,21:29,33),13,2,time] <- value
      temp[c(3,5:9,17:18,20:29,33),14,2,time] <- value
      temp[c(3:9,18:29),15,2,time] <- value
      temp[c(3:10,12:13,18:29,33),16,2,time] <- value
      temp[c(2:12,18:29,33),17,2,time] <- value
      temp[c(2:11,18:29),18,2,time] <- value
      temp[c(2:11,18:30,32),19,2,time] <- value
      temp[c(1:12,18:19,21:32),20,2,time] <- value
      temp[c(1:15,17:18,22:33,36),21,2,time] <- value
      temp[c(1:17,22:32,36),22,2,time] <- value
      temp[c(1:17,22:32,36),23,2,time] <- value
      temp[c(1:17,22:31,36),24,2,time] <- value
      temp[c(1:15,17,23:30,35:36),25,2,time] <- value
      temp[c(1:14,23:31,33:36),26,2,time] <- value
      temp[c(1:15,17,23:36),27,2,time] <- value
      temp[c(1:18,23:36),28,2,time] <- value
      temp[c(1:18,23:36),29,2,time] <- value
      temp[c(1:18,22:36),30,2,time] <- value
      temp[c(1:19,21:36),31,2,time] <- value
      temp[c(1:36),32,2,time] <- value
      temp[c(1:36),33,2,time] <- value
      temp[c(1:36),34,2,time] <- value
      temp[c(1:36),35,2,time] <- value
      temp[c(1:36),36,2,time] <- value

      temp[11:23,1,3,time] <- value
      temp[12:22,2,3,time] <- value
      temp[c(13:23,26:27),3,3,time] <- value
      temp[12:29,4,3,time] <- value
      temp[c(11:16,18:30),5,3,time] <- value
      temp[c(10:14,19:31),6,3,time] <- value
      temp[c(9:13,21:31),7,3,time] <- value
      temp[c(7:14,16,22:29,31),8,3,time] <- value
      temp[c(7:16,22:28,32),9,3,time] <- value
      temp[c(3,7:13,16,22:28,32),10,3,time] <- value
      temp[c(3,7:12, 22:28),11,3,time] <- value
      temp[c(3,7:10, 21:29),12,3,time] <- value
      temp[c(6:9,17,21:29,33),13,3,time] <- value
      temp[c(3,5:9,17:18,20:29,33),14,3,time] <- value
      temp[c(3:9,18:29),15,3,time] <- value
      temp[c(3:10,12:13,18:29,33),16,3,time] <- value
      temp[c(2:12,18:29,33),17,3,time] <- value
      temp[c(2:11,18:29),18,3,time] <- value
      temp[c(2:11,18:30,32),19,3,time] <- value
      temp[c(1:12,18:19,21:32),20,3,time] <- value
      temp[c(1:13,17:18,22:33,36),21,3,time] <- value
      temp[c(1:17,22:32,36),22,3,time] <- value
      temp[c(1:17,22:31,36),23,3,time] <- value
      temp[c(1:17,22:30,36),24,3,time] <- value
      temp[c(1:15,17,23:30,35:36),25,3,time] <- value
      temp[c(1:14,23:31,33:36),26,3,time] <- value
      temp[c(1:15,17,23:36),27,3,time] <- value
      temp[c(1:18,23:36),28,3,time] <- value
      temp[c(1:18,23:36),29,3,time] <- value
      temp[c(1:18,22:36),30,3,time] <- value
      temp[c(1:19,21:36),31,3,time] <- value
      temp[c(1:36),32,3,time] <- value
      temp[c(1:36),33,3,time] <- value
      temp[c(1:36),34,3,time] <- value
      temp[c(1:36),35,3,time] <- value
      temp[c(1:36),36,3,time] <- value
    }
    #### Permian shelves cutout
    if(grepl("p0251b", tectonic.config, fixed = TRUE) == TRUE ){
      value <- NA

      temp[c(1:14,34:36),1,1,time] <- value
      temp[c(1:13,34:36),2,1,time] <- value
      temp[c(1:13,33:36),3,1,time] <- value
      temp[c(1:12,33:36),4,1,time] <- value
      temp[c(1:13,33:36),5,1,time] <- value
      temp[c(1:13,33:36),6,1,time] <- value
      temp[c(1:13,25:26,32:36),7,1,time] <- value
      temp[c(1:13,25:27,32:36),8,1,time] <- value
      temp[c(1:13,25:28,32:36),9,1,time] <- value
      temp[c(1:12,24:29,31:36),10,1,time] <- value
      temp[c(1:12,24:29,31:36),11,1,time] <- value
      temp[c(1:12,24:28,32:36),12,1,time] <- value
      temp[c(1:12,24:25,28,32:36),13,1,time] <- value
      temp[c(1:13,24,32:36),14,1,time] <- value
      temp[c(1:13,23:25,28:29,31:36),15,1,time] <- value
      temp[c(1:13,22:28,30:36),16,1,time] <- value
      temp[c(1:13,22:27,31:36),17,1,time] <- value
      temp[c(1:13,21:26,31:36),18,1,time] <- value
      temp[c(1:14,21:26,31:36),19,1,time] <- value
      temp[c(1:14,22:27,31:36),20,1,time] <- value
      temp[c(1:13,23:28,30:36),21,1,time] <- value
      temp[c(1:13,23:36),22,1,time] <- value
      temp[c(1:14,23:29,31:36),23,1,time] <- value
      temp[c(1:14,24:28,32:36),24,1,time] <- value
      temp[c(1:14,26:27,32:36),25,1,time] <- value
      temp[c(1:14,27,32:36),26,1,time] <- value
      temp[c(1:14,27,33:36),27,1,time] <- value
      temp[c(1:14,33:36),28,1,time] <- value
      temp[c(1:14,32:36),29,1,time] <- value
      temp[c(1:14,22,31:36),30,1,time] <- value
      temp[c(1:15,21,30:36),31,1,time] <- value
      temp[c(1:15,20,29:36),32,1,time] <- value
      temp[c(1:16,20,28:36),33,1,time] <- value
      temp[c(1:17,19,28:36),34,1,time] <- value
      temp[c(1:20,28:36),35,1,time] <- value
      temp[c(1:21,28:36),36,1,time] <- value

      temp[c(1:14,34:36),1,2,time] <- value
      temp[c(1:13,34:36),2,2,time] <- value
      temp[c(1:13,33:36),3,2,time] <- value
      temp[c(1:12,33:36),4,2,time] <- value
      temp[c(1:13,33:36),5,2,time] <- value
      temp[c(1:13,33:36),6,2,time] <- value
      temp[c(1:13,25:26,32:36),7,2,time] <- value
      temp[c(1:13,25:27,32:36),8,2,time] <- value
      temp[c(1:13,25:28,32:36),9,2,time] <- value
      temp[c(1:12,24:29,31:36),10,2,time] <- value
      temp[c(1:12,24:29,31:36),11,2,time] <- value
      temp[c(1:12,24:28,32:36),12,2,time] <- value
      temp[c(1:12,24:25,28,32:36),13,2,time] <- value
      temp[c(1:13,24,32:36),14,2,time] <- value
      temp[c(1:13,23:25,28:29,31:36),15,2,time] <- value
      temp[c(1:13,22:28,30:36),16,2,time] <- value
      temp[c(1:13,22:27,31:36),17,2,time] <- value
      temp[c(1:13,21:26,31:36),18,2,time] <- value
      temp[c(1:14,21:26,31:36),19,2,time] <- value
      temp[c(1:14,22:27,31:36),20,2,time] <- value
      temp[c(1:13,23:28,30:36),21,2,time] <- value
      temp[c(1:13,23:36),22,2,time] <- value
      temp[c(1:14,23:29,31:36),23,2,time] <- value
      temp[c(1:14,24:28,32:36),24,2,time] <- value
      temp[c(1:14,26:27,32:36),25,2,time] <- value
      temp[c(1:14,27,32:36),26,2,time] <- value
      temp[c(1:14,27,33:36),27,2,time] <- value
      temp[c(1:14,33:36),28,2,time] <- value
      temp[c(1:14,32:36),29,2,time] <- value
      temp[c(1:14,22,31:36),30,2,time] <- value
      temp[c(1:15,21,30:36),31,2,time] <- value
      temp[c(1:15,20,29:36),32,2,time] <- value
      temp[c(1:16,20,28:36),33,2,time] <- value
      temp[c(1:17,19,28:36),34,2,time] <- value
      temp[c(1:20,28:36),35,2,time] <- value
      temp[c(1:21,28:36),36,2,time] <- value

      temp[c(1:14,34:36),1,3,time] <- value
      temp[c(1:13,34:36),2,3,time] <- value
      temp[c(1:13,33:36),3,3,time] <- value
      temp[c(1:12,33:36),4,3,time] <- value
      temp[c(1:13,33:36),5,3,time] <- value
      temp[c(1:13,33:36),6,3,time] <- value
      temp[c(1:13,32:36),7,3,time] <- value
      temp[c(1:13,32:36),8,3,time] <- value
      temp[c(1:13,27,32:36),9,3,time] <- value
      temp[c(1:12,25:26,31:36),10,3,time] <- value
      temp[c(1:12,25,31:36),11,3,time] <- value
      temp[c(1:12,25:26,32:36),12,3,time] <- value
      temp[c(1:12,25,32:36),13,3,time] <- value
      temp[c(1:13,32:36),14,3,time] <- value
      temp[c(1:13,25,29,31:36),15,3,time] <- value
      temp[c(1:13,25:26,30:36),16,3,time] <- value
      temp[c(1:13,24:26,31:36),17,3,time] <- value
      temp[c(1:13,24:26,31:36),18,3,time] <- value
      temp[c(1:13,23:26,31:36),19,3,time] <- value
      temp[c(1:13,23:27,31:36),20,3,time] <- value
      temp[c(1:13,23:28,30:36),21,3,time] <- value
      temp[c(1:13,23:36),22,3,time] <- value
      temp[c(1:13,23:29,31:36),23,3,time] <- value
      temp[c(1:14,24:28,32:36),24,3,time] <- value
      temp[c(1:14,26:27,32:36),25,3,time] <- value
      temp[c(1:14,27,32:36),26,3,time] <- value
      temp[c(1:14,27,33:36),27,3,time] <- value
      temp[c(1:14,33:36),28,3,time] <- value
      temp[c(1:14,32:36),29,3,time] <- value
      temp[c(1:14,31:36),30,3,time] <- value
      temp[c(1:15,30:36),31,3,time] <- value
      temp[c(1:15,29:36),32,3,time] <- value
      temp[c(1:15,20,28:36),33,3,time] <- value
      temp[c(1:16,19,28:36),34,3,time] <- value
      temp[c(1:16,19:20,28:36),35,3,time] <- value
      temp[c(1:15,20:21,28:36),36,3,time] <- value
    }
    #### Paleocene shelves cutout
    if(grepl("p0055c", tectonic.config, fixed = TRUE) == TRUE ){
      value <- NA

      # No cells to set to NA at latitude 1
      # No cells to set to NA at latitude 2
      temp[c(1:10,14:20,34:36),3,1,time] <- value
      temp[c(1:10,15:26,34:36),4,1,time] <- value
      temp[c(1:10,14:26,34:36),5,1,time] <- value
      temp[c(1:10,14:18,20:26,34:36),6,1,time] <- value
      temp[c(1:10,15:17,21:27,34:36),7,1,time] <- value
      temp[c(1:10,15:17,22:28,33:36),8,1,time] <- value
      temp[c(1:10,16:17,22:29,33:36),9,1,time] <- value
      temp[c(1:10,16:17,22:29,33:36),10,1,time] <- value
      temp[c(1:10,16:17,22:30,32:36),11,1,time] <- value
      temp[c(1:10,17,23:36),12,1,time] <- value
      temp[c(1:10,17,23:36),13,1,time] <- value
      temp[c(1:9,17,22:36),14,1,time] <- value
      temp[c(1:9,17,22:23,25:36),15,1,time] <- value
      temp[c(1:9,26:36),16,1,time] <- value
      temp[c(1:9,16,26:36),17,1,time] <- value
      temp[c(1:9,16,27:28,30:36),18,1,time] <- value
      temp[c(1:9,15,26:27,31:36),19,1,time] <- value
      temp[c(1:10,14:15,26:27,31:36),20,1,time] <- value
      temp[c(1:11,13:14,23,25:26,31:36),21,1,time] <- value
      temp[c(1:15,23:26,30:36),22,1,time] <- value
      temp[c(1:8,10:15,23:24,30:36),23,1,time] <- value
      temp[c(1:7,11:14,23,30:36),24,1,time] <- value
      temp[c(1:7,11:14,30:36),25,1,time] <- value
      temp[c(1:7,11:15,20,22,31:36),26,1,time] <- value
      temp[c(1:7,11:16,20:22,31:36),27,1,time] <- value
      temp[c(1:6,12:16,19:21,31:36),28,1,time] <- value
      temp[c(1:6,13:15,20,32:36),29,1,time] <- value
      temp[c(1:6,15:16,21,23,32:36),30,1,time] <- value
      temp[c(1:6,21:23,32:36),31,1,time] <- value
      temp[c(1:6,15,33:36),32,1,time] <- value
      temp[c(1:5,33:36),33,1,time] <- value
      temp[c(1:2,35:36),34,1,time] <- value
      # No cells to set to NA at latitude 35
      temp[c(23),36,1,time] <- value

      # No cells to set to NA at latitude 1
      # No cells to set to NA at latitude 2
      temp[c(1:10,14:20,34:36),3,2,time] <- value
      temp[c(1:10,15:26,34:36),4,2,time] <- value
      temp[c(1:10,14:26,34:36),5,2,time] <- value
      temp[c(1:10,14:18,20:26,34:36),6,2,time] <- value
      temp[c(1:10,15:17,21:27,34:36),7,2,time] <- value
      temp[c(1:10,15:17,22:28,33:36),8,2,time] <- value
      temp[c(1:10,16:17,22:29,33:36),9,2,time] <- value
      temp[c(1:10,16:17,22:29,33:36),10,2,time] <- value
      temp[c(1:10,16:17,22:30,32:36),11,2,time] <- value
      temp[c(1:10,17,23:36),12,2,time] <- value
      temp[c(1:10,17,23:36),13,2,time] <- value
      temp[c(1:9,17,22:36),14,2,time] <- value
      temp[c(1:9,17,22:23,25:36),15,2,time] <- value
      temp[c(1:9,26:36),16,2,time] <- value
      temp[c(1:9,16,26:36),17,2,time] <- value
      temp[c(1:9,16,27:28,30:36),18,2,time] <- value
      temp[c(1:9,15,26:27,31:36),19,2,time] <- value
      temp[c(1:10,14:15,26:27,31:36),20,2,time] <- value
      temp[c(1:11,13:14,23,25:26,31:36),21,2,time] <- value
      temp[c(1:15,23:26,30:36),22,2,time] <- value
      temp[c(1:8,10:15,23:24,30:36),23,2,time] <- value
      temp[c(1:7,11:14,23,30:36),24,2,time] <- value
      temp[c(1:7,11:14,30:36),25,2,time] <- value
      temp[c(1:7,11:15,20,22,31:36),26,2,time] <- value
      temp[c(1:7,11:16,20:22,31:36),27,2,time] <- value
      temp[c(1:6,12:16,19:21,31:36),28,2,time] <- value
      temp[c(1:6,13:15,20,32:36),29,2,time] <- value
      temp[c(1:6,15:16,21,32:36),30,2,time] <- value
      temp[c(1:6,21:22,32:36),31,2,time] <- value
      temp[c(1:6,15,33:36),32,2,time] <- value
      temp[c(1:5,33:36),33,2,time] <- value
      temp[c(1:2,35:36),34,2,time] <- value
      # No cells to set to NA at latitude 35
      temp[c(23),36,2,time] <- value

      # No cells to set to NA at latitude 1
      # No cells to set to NA at latitude 2
      temp[c(1:10,14:20,34:36),3,3,time] <- value
      temp[c(1:10,15:26,34:36),4,3,time] <- value
      temp[c(1:10,14:26,34:36),5,3,time] <- value
      temp[c(1:10,14:18,20:26,34:36),6,3,time] <- value
      temp[c(1:10,15:17,21:27,34:36),7,3,time] <- value
      temp[c(1:10,15:17,22:28,33:36),8,3,time] <- value
      temp[c(1:10,16:17,22:29,33:36),9,3,time] <- value
      temp[c(1:10,16:17,22:29,33:36),10,3,time] <- value
      temp[c(1:10,16:17,22:30,32:36),11,3,time] <- value
      temp[c(1:10,17,23:36),12,3,time] <- value
      temp[c(1:10,17,23:36),13,3,time] <- value
      temp[c(1:9,17,22:36),14,3,time] <- value
      temp[c(1:9,17,22:23,25:36),15,3,time] <- value
      temp[c(1:9,26:36),16,3,time] <- value
      temp[c(1:9,16,26:36),17,3,time] <- value
      temp[c(1:9,16,27:28,30:36),18,3,time] <- value
      temp[c(1:9,15,27,31:36),19,3,time] <- value
      temp[c(1:10,14:15,27,31:36),20,3,time] <- value
      temp[c(1:11,13:14,31:36),21,3,time] <- value
      temp[c(1:15,26,30:36),22,3,time] <- value
      temp[c(1:8,10:15,23:24,30:36),23,3,time] <- value
      temp[c(1:7,11:14,23,30:36),24,3,time] <- value
      temp[c(1:7,11:14,30:36),25,3,time] <- value
      temp[c(1:7,11:15,22,31:36),26,3,time] <- value
      temp[c(1:7,11:16,20:22,31:36),27,3,time] <- value
      temp[c(1:6,12:16,19:21,31:36),28,3,time] <- value
      temp[c(1:6,13:15,20,32:36),29,3,time] <- value
      temp[c(1:6,15:16,21,32:36),30,3,time] <- value
      temp[c(1:6,32:36),31,3,time] <- value
      temp[c(1:6,15,33:36),32,3,time] <- value
      temp[c(1:5,33:36),33,3,time] <- value
      temp[c(1:2,35:36),34,3,time] <- value
      # No cells to set to NA at latitude 35
      # No cells to set to NA at latitude 36
    }
    # Shelf cutouts over.
    }else{
    # If not cutting out surface or shelves - set maximum z dimension to number of depth levels
    zdim <- as.numeric(length(depth))
    }


    ## =============================================================
    ## make array of pO2 - units are atm
    ## =============================================================
    oxygen_umol_kg <- oxy*10^6  # convert to umol/kg from mol/kg

    depth_array <- array(rep(depth, each=xdim*ydim), dim=c(xdim,ydim,zdim))
    lat_array <- array(rep(lat, each = xdim, times = zdim), dim=c(xdim,ydim,zdim))
    ## pass oxygen, salinity, temp and depth to Hoffman et al. 2011's pO2 function
    pO2_array <- pO2(O2=oxygen_umol_kg[,,,time], S=sal[,,,time], t=temp[,,,time], d=depth_array, lat=lat_array)/1000


    ## =============================================================
    ## Convert temp to Kelvin
    ## =============================================================

    temp_K <- temp[,,,time]+273.15

    ## =============================================================
    ## Habitat viability model - calculate viability per cGENIE model cell
    ## =============================================================

    # initiate array for habitat viability (only saved if NetCDF saving activated)
    habitat_viability <- array(dim=c(xdim,ydim,zdim))

    # initiate data frame for viable ecotype summary (one per cGENIE simultion)
    phi.crit.xxx.summary <-  data.frame(A0.xxx = double(), phi_crit.xxx=double(), temp_K=double(), pO2=double(), phi.xxx=double(), ecotype=double(), depth=double())


    # loop through by cell (we are evaluating xxx physiological ecotypes per cell)
    for(x in 1:xdim){
      for(y in 2:(ydim-1)){
        for(z in 1:zdim){
          pO2.cell <- pO2_array[x,y,z]
          temp.cell <- temp_K[x,y,z]

          if(is.na(pO2.cell) == TRUE | is.na(temp.cell) == TRUE){
            # Do nothing! If these are true - then this is land/crust!

          }else{
            # Calculate phi values for all xxx ecotypes (eqn. from Penn et al. 2018)
            phi.xxx <- A0.xxx * (pO2.cell/(exp((-E0.xxx/kB)*((1/temp.cell)-(1/Tref)))))

            # Combine relevant data into a dataframe
            phi.cell <- cbind(A0.xxx, phi_crit.xxx, rep(temp.cell, init.ecotypes), rep(pO2.cell, init.ecotypes), phi.xxx, seq(1:init.ecotypes), depth[z])

            # Name variables in cell ecotypes dataframe so that we can refer to them easily
            names(phi.cell) <- c("A0.xxx", "phi_crit.xxx", "temp_K", "pO2", "phi.xxx", "ecotype", "depth")

            # Filter the cell ecotypes dataframe to only include viable ecotypes (phi > phi_crit)
            phi.crit.xxx <- phi.cell %>%
                            as.data.frame()  %>%
                            filter(phi.xxx >= phi_crit.xxx)

            # Add habitat viability (percentage of total ecotypes that are viable in a given cell) to summary array (for building NetCDFs if so inclined)
            habitat_viability[x,y,z] <- nrow(phi.crit.xxx)/init.ecotypes*100

            # Add viable ecotypes and associated details to summary dataframe for specific cGENIE simulation
            # This has been vectorized (at the loss of some readability) for speed
            if(nrow(phi.crit.xxx)>0){
            phi.crit.xxx.summary[(nrow(phi.crit.xxx.summary)+1):(nrow(phi.crit.xxx.summary)+nrow(phi.crit.xxx)) , ] <- phi.crit.xxx
            }
          }
        }}}

    # Tally the number of cells that each ecotype can live in
    phi.crit.xxx.summary.condensed <- phi.crit.xxx.summary %>%
                                      group_by(ecotype) %>%
                                      tally

    # Add this tally to a summary frame - again, this has been vectorized (at the loss of some readability) for speed
    if(nrow(phi.crit.xxx.summary.condensed)>0){
    i <- nrow(ecotypes.summary)+1
    j <- nrow(ecotypes.summary)+nrow(phi.crit.xxx.summary.condensed)
    ecotypes.summary[i:j , 1:2] <- phi.crit.xxx.summary.condensed
    ecotypes.summary[i:j , 3] <- O2.factor
    ecotypes.summary[i:j , 4] <- CO2.factor
    ecotypes.summary[i:j , 5] <- it
    }

    # Summarise this data further for habitat viability plots (1 point per cGENIE simulation) - here add equatorial temperature for plotting.
    row <- nrow(ecotypes.summary.short)+1
    ecotypes.summary.short[row , 1] <- nrow(phi.crit.xxx.summary.condensed)
    ecotypes.summary.short[row , 2] <- O2.factor
    ecotypes.summary.short[row , 3] <- CO2.factor
    ecotypes.summary.short[row , 4] <- eq.temp
    ecotypes.summary.short[row , 5] <- it
    
    ## =============================================================
    ## Make NetCDF (if option selected towards beginning of script)
    ## =============================================================

    if(habitat.viability.netcdf == 1){

    # Define the name of netcdf you will make below
    filename <- paste0(netcdf.directory, tectonic.config, CO2.factor, "CO2.", O2.factor, "O2.", EFD, "EFD.", ocean.mask, ".habitat.viability.nc")

    # Define dimensional variables ready to add them to your array
    xvals <- lon
    yvals <- lat
    nx <- length(xvals)
    ny <- length(yvals)
    lon1 <- ncdim_def("longitude", "degrees_east", xvals)
    lat2 <- ncdim_def("latitude", "degrees_north", yvals)
    depth3 <- ncdim_def("Depth", "meters", depth[1:zdim])
    mv <- NA #missing value to use (for land/crust areas)

    # Create a Habitat Viability variable (with dimensions defined above)
    var_via <- ncvar_def("Habitat Viability", paste(expression("(%)")), list(lon1, lat2, depth3), longname="Habitat Viability", mv)

    # Create a NetCDF file with the file name and Habitat Viability variable defined above
    ncnew <- nc_create(filename, list(var_via))

    # Add your Habitat Viability array (made in loop) above to the corresponding variable in your new NetCDF file
    ncvar_put(ncnew, var_via, habitat_viability)

    # Close the new NetCDF file - it is ready to open in your chosen NetCDF reader in the directory defined above!
    nc_close(ncnew)
    }

    # Print O2 and CO2 levels of last finished cGENIE model simulation
    # (adds a negligible amount of additional runtime, but great for checking progress!)
    print(O2.factor)
    print(CO2.factor)
    print(paste("iteration number", it))
  }}
))}

# Save the ecotype summary files generated in the extractions and manipulations performed above
save(ecotypes.summary, file=paste0(rdata.directory,config.short,".CO2.O2.", EFD, "EFD ecotypes summary ", ocean.mask, ".ecotype.sensitivity.no.poles_lowAo_highEo.Rdata"))
save(ecotypes.summary.short, file=paste0(rdata.directory,config.short,".CO2.O2.", EFD, "EFD ecotypes summary ", ocean.mask, " short.ecotype.sensitivity.no.poles_lowAo_highEo.Rdata"))

