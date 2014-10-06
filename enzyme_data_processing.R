# First, read in the files that you need to calculate enzyme calculations.  This script is based on the enzyme protocol last modified by SB on 07/2014
# The following are needed:
# 1. The assay plate (soil homogenate + substrates and std. curves)
# 2. The blank plate (no soil, just substrate blanks and std. curves...there should only be one a day)

assay_plate<-read.csv(file.choose())
blank_plate<-read.csv(file.choose(),row.names=1,check.names=FALSE)


# Check to make sure that the files you read in look fine.

head(assay_plate)
head(blank_plate)

# The steps below will walk you through how to perform the calculations for an individual sample.  Immediately following will be steps to perform the calculations for a series of samples.

# We need the actual mass of MUB and MUC weighed out on the microbalance, change this per batch of standards made
act_mub<-.88
act_muc<-.88

# The first step is to calculate concentrations of substrates
MUB_blank_slope(blank_plate, act_mub)
MUC_blank_slope(blank_plate, act_muc)

# If the plots look good, assign the slope to an object
mub_blank_slope<-MUB_blank_slope(blank_plate, act_mub)[[2]]
muc_blank_slope<-MUC_blank_slope(blank_plate, act_muc)[[2]]

# Next we calculate the MUB and MUC slope on the assay plate, note that I specified the rows for the first sample here
# These functions will also not produce plots, though this can be changed if wanted

mub_assay_slope<-MUB_assay_slope(assay_plate[1:8,],act_mub)
muc_assay_slope<-MUC_assay_slope(assay_plate[1:8,],act_muc)

# Next we calculate the quench and emission coefficients 
mub_quench<-mub_assay_slope/mub_blank_slope
muc_quench<-muc_assay_slope/muc_blank_slope

mub_emission<-mub_blank_slope/.25
muc_emission<-muc_blank_slope/.25

# Now we calculate the activity from the plate, here we will just do BG (column 1 of the plate)  the function below will do all of the appropriate columns based on the plate layout from the protocol

activity(assay_plate[1:8,],1,mub_quench, mub_emission,mean(assay_plate[,3]),mean(assay_plate[,5]))

# Now we want to be able to do it all in one loop...first all parameters need to be calculated for the blank plate as seen above

enzyme_loop(assay_plate)

# You can put that in an object (it returns a data.frame) and run your stats!  This would be a good time to run a test dataset that has already been through other calculations
