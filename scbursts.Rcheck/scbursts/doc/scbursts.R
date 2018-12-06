## ---- fig.show='hold', eval=FALSE, include = TRUE------------------------
#  # Load the library
#  library(scbursts)
#  
#  infile <- "data/100uM.evt"
#  
#  # Import the evt as a table
#  transitions <- evt.read(infile)
#  
#  # Turn the transition times into dwells
#  dwells <- evt.to_dwells(transitions)
#  
#  # With dwells defined, we can start doing an actual analysis

## ---- fig.show='hold', eval=FALSE, include = TRUE------------------------
#  # Write the corrected transition times to disk.
#  evt.write(dwells_corrected, file="100uMc.evt")

## ---- fig.show='hold', eval=FALSE, include = TRUE------------------------
#  dwells <- dwt.read("60uM.dwt")
#  
#  # ...
#  #
#  # Correct the dwells or do an analysis
#  #
#  # ...
#  
#  dwt.write(corrected_dwells, file="60uMc.dwt")

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# Load the library
library(scbursts)

# Import a pre-packaged file (stored inside the folder extdata)
dwt_example <- system.file("extdata", "example.dwt", package = "scbursts")

# Import the evt as a table
dwells <- dwt.read(dwt_example)

# Just transition times and states
head(dwells[[1]])

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# Load the library
library(scbursts)

# Import a pre-packaged file (stored inside the folder extdata)
infile <- system.file("extdata", "example.evt", package = "scbursts")


# Import the evt as a table
tables <- evt.read(infile)

# Turn the transition times into dwells
records <- evt.to_dwells(tables)

# Correct the risetime (default time in seconds)
records_c <- risetime.correct_gaussian(Tr=35.0052278,records, unit="us")

length(records_c)

# Define critical time (tcrit=0.1s)
bursts <- bursts.defined_by_tcrit(records_c , 0.1, units="ms")

length(bursts)

## Now you can carry out analysis of the bursts

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# Load the library
library(scbursts)

# Import a pre-packaged file (stored inside the folder extdata)
infile <- system.file("extdata", "example.evt", package = "scbursts")


# Import the evt as a table
tables <- evt.read(infile)
records <- evt.to_dwells(tables)

# Correct the risetime (default time in seconds)
records_c <- risetime.correct_gaussian(Tr=35.0052278,records, unit="us")

# Define critical time (tcrit=0.1s)
bursts <- bursts.defined_by_tcrit(records_c , 0.1, units="ms")

high_popen <- function (seg) {
    segment.popen(seg) > 0.7
}

high_bursts <- bursts.select(bursts, high_popen)

## ---- fig.show='hold', eval=FALSE, include = TRUE------------------------
#  high_bursts <- bursts.select(bursts, high_popen, one_file=TRUE)

## ---- fig.show='hold', eval=FALSE, include = TRUE------------------------
#  library(scbursts)
#  infile <- "faulty_dwells.dwt"
#  
#  # This will raise a warning message to alert you that the data has problems
#  dwells <- dwt.read(infile)
#  
#  dwells_c <- risetime.correct_gaussian(Tr=0.0000350052278,dwells)
#  
#  # This will also raise a warning message to alert you that specific bursts have problems
#  bursts <- bursts.defined_by_tcrit(dwells_c,0.1)
#  
#  # This will remove the problems. It will leave only the good bursts.
#  bursts <- bursts.select(bursts, segment.verify)

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# If you have multiple records, you can recombine them with
# This is now just one list of spaced out segments.
records <- bursts.space_out(records, sep_factor=1000)
record <- bursts.recombine(records)

## ---- eval=TRUE, include = TRUE------------------------------------------
# Create a list of bursts, sorted by your chosen function
sorted <- bursts.sort(bursts, segment.popen, reverse=TRUE)

# In some cases, it might be that multiple bursts share the same value
# and so the "order" is a bit arbitrary in those cases.
sorted[[1]]

## ---- eval=TRUE, include = TRUE------------------------------------------
mean(bursts.popens(bursts))

## ---- eval=FALSE, include = TRUE-----------------------------------------
#  bursts.popens <- function (bursts) { sapply(bursts, segment.popen) }

## ---- eval=TRUE, include = TRUE------------------------------------------
mean(sapply(bursts, segment.duration))

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# Correct the risetime

corrected_records <- list()
for (i in 1:length(records)) {
    corrected_records[[i]] <- risetime.correct_gaussian(Tr=35.0052278, records[[i]], units="us")
}

# Write the corrected record to a .dwt file
dwt.write(corrected_records, file="60uMc.dwt")

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# Correct the risetime
records_c <- risetime.correct_gaussian(Tr=35.0052278, records, units="us")

# Write the corrected record to a .dwt file
dwt.write(records_c, file="60uMc.dwt")

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# Load the library
library(scbursts)

# Import a pre-packaged file (stored inside the folder extdata)
infile <- system.file("extdata", "example.evt", package = "scbursts")


# Import the evt as a table
tables <- evt.read(infile)

# Turn the transition times into dwells
records <- evt.to_dwells(tables)

# Correct the risetime (default time in seconds)
records_c <- risetime.correct_gaussian(Tr=35.0052278,records, unit="us")

evt.write(records_c, file="example_corrected.evt")

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
open_dwells <- segment.open_dwells(record) / 1000
hist(log10(open_dwells), axes=FALSE, breaks=30)
cplot.log_root_axes(open_dwells)

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
closed_dwells <- segment.closed_dwells(record) / 1000
hist(log10(closed_dwells), axes=FALSE, breaks=30)
cplot.log_root_axes(closed_dwells)

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
popens <- bursts.popens(bursts)
hist(popens, breaks=20)

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
pcloseds <-bursts.pcloseds(bursts)
hist(pcloseds, breaks=20)

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# To make this more visible, you can also export it as a large `.png` file
cplot.popen_ts(bursts)

## ---- fig.show='hold', eval=TRUE, include = TRUE-------------------------
# Or, so look at a subregion 
cplot.popen_ts(bursts, xlim=c(0,0.1))

