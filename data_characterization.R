#
# data_characterization.R
# ADAPTED FROM: gsi.R

#######################################################################
# GSI
# 
# by Rita Ribeiro, LIACC (rita@liacc.up.pt)
#    with a small collaboration of Luis Torgo (ltorgo@liacc.up.pt) and
#    Carlos Soares (csoares@liacc.up.pt)
#
# ----------------------------------------------------------------------
#
# This file contains R functions to support characterization of datasets 
# in METAL format (data and names files)
#
# ----------------------------------------------------------------------

# BY: Carlos Soares

# COMMENTS:
# -assumes two files to describe the data set
# -minimal check strategy: assumes data is in correct format, problem
#  is of correct type, etc...
# -to implement extensibility, several calculations are repeated
# -not applicable values are printed as "?"

# CHANGES:
# (14/05/04 carlos) writeFlat.default can handle both matrices and arrays.
# (11/06/02 carlos) writeFlat.table is not specific for attribute/class
# measures
# (22/05/02 carlos) v1.1: METAL regression file format supported by 
# expecting  option target.declared in the structure describing the 
# measures; added support to print out matrices; temporary implementation 
# of binarization of symbolic attributes
# (19/03/02 carlos) data set files can be in directory different than R sources
# (14/03/02 carlos) v1.0: changed structure to simplify extensibility: 
# reading dataset and characterizing it are now separate (actually 
# in different files) and the latter uses the former; data characterization
# is done using a simple API defined in this file, which includes timing 
# the calculation of data characteristics; calculating characteristics and 
# printing them is also separate; provides functions to identify symbolic 
# and continuous attributes; provide support for version information (of data
# characterization measures as well)

# TODO:
# -binarization: move to read_data.R: user assesses attributes using lists of attributes included in the data description (cont.attrs, symb.attrs, original.attrs, binarized.attrs), replacing functions in ContAttrs and SymbAttrs; update original names as well
# -implement parameters to the data characterization as structured object
# -improve support for different sets of measures and different type of problems: choosing regression/classification with target.declared option is a good solution?
# -add more parameters of ReadDataSet to CharacterizeDataSet
# -print.out is not really necessary
# -extensibility: how to include other parameters in the calls to functions? how to name arguments of data characterization functions to avoid positional arguments: should be possible with do.call
# -extensibility: should we to take time of already calculated measures into account? if so, how?
# -writeFlat: should tables be printed class-first?

# HOW TO USE:
# The main function is CharacterizeDataSet. It essentially receives a name of 
# a data set and a set of measures, and calculates all of the latter on 
# the former. It essentially returns a list containing the values computed 
# and the times required to obtain them.
# The functions that calculate measures are named after them and have a 
# specific signature. They may return a simple value, a structured value 
# or NA, if the measure is not applicable. To perform calculations, the 
# object containing the data and its description, will often be necessary. 
# A few functions are also provided, e.g. to obtain previously calculated 
# measures (using GetMeasure) and lists of attributes (ContAttrs and SymbAttrs).
# Several functions are provided to display or save a data characterization.
# Function writeMETAL, in particular, uses a format which is suitable for 
# J. Petrak's parse_results script.

# MORE INFORMATION:
# -gsi.R for examples
# -read_data.R for the data set object structure
# -SUPPORT FOR CALCULATING MEASURES section, below
# -DISPLAYING AND SAVING DATA CHARACTERIZATIONS section, below
# -the function headers contain useful information, in general
# -if everything else fails, contact csoares@liacc.up.pt

# ######################################################################

# LIBRARIES
source("read_data.R") # To enable working in a directory different than the one containing the scripts

kDCTVersion <- "1.1"

# ======================================================================
CharacterizeDataSet <- 		# was BuildGSI
	function(inFilestem,
		dc.measures,
		anonymize = TRUE) 

# Calculates a set of data measures on the given data set

# IN:
# -inFilestem: filestem of data files
# -dc.measures: a list with components:
#  -name: name of the type of characterization
#  -version
#  -target.declared
#  -binarize
#  -measures: list of measure names
# -anonymize: anonymize class values, attribute names and values

# RETURNS:
# -list with components:
#  -value: a list of simple or structure measure values
#  -time: a vector of times to calculate the measures (note that it is 
#   a vector not a list!)
#  -version: a list with names representing the software and values
#   representing versions

# COMMENTS:
# -regarding the measure calculation functions:
#   -named after the corresponding measure
#   -arguments: the data set object and the current data characteristics
#    object,
#   -return: a value, a more complex structure (vector, list, etc) or NA
#    if the measure is not applicable
#   -to use the value of another measure in a calculation, get it with 
#    GetMeasure

# ***
{
# Setup data structures
wkDCSet <- list(value = list(), time = list(), version = list(dct = kDCTVersion))
wkDCSet$version[[dc.measures$name]] <- dc.measures$version

# Read data set
wkDCSet$time$read.data <- system.time(wkDataSet <- ReadDataSet(inFilestem, target.declared = dc.measures$target.declared, anonymize = anonymize))[1]
 
if (! is.null(dc.measures$binarize))
	{
	# Binarize symbolic attributes, individually
	wkBinAttrNames <- c(sapply(SymbAttrs(wkDataSet), function(x, df) {paste(x,levels(df[, x])[-1],sep=".")}, wkDataSet$frame, USE.NAMES=FALSE), recursive=TRUE)
	wkBinAttrs <- matrix(0, ncol=length(wkBinAttrNames), nrow=nrow(wkDataSet$frame), dimnames=list(NULL, wkBinAttrNames))
	for (wkName in wkBinAttrNames)
		{
		wkDotIndex <- regexpr("[.]", wkName)[1]
		wkBinAttrs[wkDataSet$frame[[substr(wkName, start=1, stop=wkDotIndex - 1)]] == substr(wkName, start=wkDotIndex + 1, stop=nchar(wkName)), wkName] <- 1
		}

	# Replace symbolic attributes with binarized ones
	wkBinAttrsDF <- data.frame(wkBinAttrs)
	wkContAttrNames <- ContAttrs(wkDataSet)
	wkDataSet$frame <- data.frame(wkDataSet$frame[, wkDataSet$attributes$attr.name[! wkDataSet$attributes$attr.name %in% SymbAttrs(wkDataSet)]], wkBinAttrsDF)

	# Correct attribute descriptions
	wkNewNames <- unlist(sapply(wkDataSet$attributes$attr.name, function(x,type) {if (type[[x]] != "continuous") paste(x, type[[x]][-1], sep=".") else x}, wkDataSet$attributes$attr.type))
	colnames(wkDataSet$frame) <- c(wkDataSet$attributes$attr.name[! wkDataSet$attributes$attr.name %in% SymbAttrs(wkDataSet)], wkBinAttrNames)
	wkDataSet$frame <- wkDataSet$frame[, wkNewNames]

	if (anonymize)
		{
		# Re-anonymize attribute names
		wkTargetType <- wkDataSet$attributes$attr.type[wkDataSet$attributes$target.attr]
		wkDataSet$attributes$target.attr <- as.character(match(wkDataSet$attributes$target.attr, colnames(wkDataSet$frame)))
		colnames(wkDataSet$frame) <- as.character(1:length(colnames(wkDataSet$frame)))
		wkDataSet$attributes$attr.name <- colnames(wkDataSet$frame) 
		wkNewTypes <- rep(list("continuous"), times=length(colnames(wkDataSet$frame)))
		names(wkNewTypes) <- colnames(wkDataSet$frame)
		wkNewTypes[wkDataSet$attributes$target.attr] <- wkTargetType
		wkDataSet$attributes$attr.type <- wkNewTypes
		}
	}

# Calculate all measures in the list
for (wkMeasure in dc.measures$measures)
	{
	if (is.null(wkDCSet$value[[wkMeasure]]))
		{
		wkTime <- system.time(wkValue <- do.call(wkMeasure, list(wkDataSet, wkDCSet)))
		wkDCSet$value[[wkMeasure]] <- wkValue
		wkDCSet$time[[wkMeasure]] <- wkTime[1]
		}
	}

wkDCSet
}

# ****************** SUPPORT FOR CALCULATING MEASURES ******************

# ======================================================================
GetMeasure <- 
	function(inDCName,
		inDCSet,
		component.name = "value")

# Retrieve the value of a previously computed measure

# IN:
# -inDCName: name of data characteristics
# -inDCSet: set of data characteristics already computed
# -component.name: name of component (e.g. time or value) to retrieve; if NULL
#  retrieve all

# OUT:
# -simple or structured value

# COMMENTS:
# -if measure is not available, stop execution with error

# ***
{
if (is.null(inDCSet$value[[inDCName]]))
	stop(message = "WARNING:requires uncomputed measure (", inDCName, ")")

if (is.null(component.name))
	inDCSet[[inDCName]]
else
	inDCSet[[component.name]][[inDCName]]
}

# ======================================================================
SymbAttrs <-
	function(dataset)

# Retrieve names of symbolic attributes (not including the target)

# IN:
# -dataset: structure describing the data set (see read_data.R)

# OUT:
# -list of strings

# ***

{
dataset$attributes$attr.name[(dataset$attributes$attr.type != "continuous") & (dataset$attributes$attr.name != dataset$attributes$target.attr)]
}

# ======================================================================
ContAttrs <-
	function(dataset)

# Retrieve names of continuous attributes (not including the target)

# IN:
# -dataset: structure describing the data set (see read_data.R)

# OUT:
# -list of strings

# ***
{
dataset$attributes$attr.name[(dataset$attributes$attr.type == "continuous") & (dataset$attributes$attr.name != dataset$attributes$target.attr)]
}

# ************ DISPLAYING AND SAVING DATA CHARACTERIZATIONS ************

# ======================================================================
writeMETAL <- 
	function(inDC,	
		fd = stdout())

# Write a set of data charateristics in a format suitable for J. Petrak's
# parse_results script

# IN:
# -inDC: data characteristization as returned by CharacterizeDataSet
# -fd: file descriptor to send the output to

# ***
{
writeVersions(inDC$version, fd = fd)
writeFlat(inDC$value, fd = fd)
print.out(fd, paste("read_time", inDC$time$read.data, sep =": "))
print.out(fd, paste("total_time", sum(unlist(inDC$time)), sep =": ")) # I wonder if unlist is intended for this...
}

# ======================================================================
writeFlat <- 
	function(inList,
		inCurrentMeasureName = NULL,
		fd = stdout())

# Generic function to unstructures a set of data characteristics and
# writes it in the format NAME: VALUE, suitable for J. Petrak's 
# parse_results script format, if data is anonymized

# IN:
# -inList: list containing simple, structured (lists or tables)
#  and NA values
# -inCurrentMeasureName: concatenation of the names of the previously
#  unstructed parts of the current measure, which will be used to 
#  identify it
# -fd: file descriptor to send the output to

# COMMENTS:
# -currently implemented:
#  -default function (handles vectors, matrices and lists)
#  -table

# ***
{
UseMethod("writeFlat")
}

# ======================================================================
writeFlat.default <- 
	function(inList,
		inCurrentMeasureName = NULL,
		fd = stdout())

# COMMENTS: see generic function above

# ***
{
wkComponentNames <- names(inList)

if(is.null(wkComponentNames))
	{
	if (is.na(inList))
		print.out(fd, paste(gsub("\\.", "_", as.character(inCurrentMeasureName)), "?", sep=": "))
	else if (is.array(inList))
          sapply(1:dim(inList)[1], function(i, m, fd)
                 {
                   writeFlat(do.call("[", c(list(m), list(i), lapply(dim(m)[-1], function(l) {1:l}))), paste(c(inCurrentMeasureName, i), sep = "_"), fd)
                 }, inList, fd)
	else
		print.out(fd, paste(gsub("\\.", "_", as.character(inCurrentMeasureName)), inList, sep=": "))
	} 
else 
	{
	for(wkElement in wkComponentNames) 
		{
		writeFlat(inList[[wkElement]], paste(c(inCurrentMeasureName, wkElement), sep = "", collapse = "_"), fd)
		}
	}
}

# ======================================================================
writeFlat.table <- 
	function(inTable,
		inCurrentMeasureName = NULL,
		fd = stdout())

# COMMENTS:
# -assumes matrix with single value in each position
# -prints column-first
# -see generic function above

# ***
{
# TODO:(28/07/23 carlos) bug corrected?
if (length(dim(inTable)) == 1) # table with single row
{
	col <- 1
	for (row in names(inTable))
		print.out(fd, paste(gsub("\\.", "_", paste(as.character(inCurrentMeasureName), col, row, sep = "_")), inTable[row], sep=": "))	
}
else
{
	for (col in colnames(inTable))
		for (row in rownames(inTable))
			print.out(fd, paste(gsub("\\.", "_", paste(as.character(inCurrentMeasureName), col, row, sep = "_")), inTable[row, col], sep=": "))
}
}

# ======================================================================
writeVersions <-
	function(versions,
		 fd = stdout())

# Write DCT and data characteristics versions used

# IN:
# -versions: structure with version information (see CharacterizeDataSet)
# -fd: file descriptor to send the output to
 
# ***
{
for (wkComponent in names(versions))
	print.out(fd, paste(wkComponent, ": version ", versions[[wkComponent]], sep =""))
}

# ======================================================================
print.out <- 
	function(fd,...) 
# ***
{
	writeLines(paste(...,sep=""),con=fd)
}
