#######################################################################
# GSI
# 
# by Rita Ribeiro, LIACC (rita@liacc.up.pt)
#    wit a small collaboration of Luis Torgo (ltorgo@liacc.up.pt) and
#    Carlos Soares (csoares@liacc.up.pt)
#
# ----------------------------------------------------------------------
#
# This file a set of general, statistical and information-theoretic
# measures for data characterization, generated using the R data
# characterization API (see data_characterization.R)
#
# ----------------------------------------------------------------------

# COMMENTS:
# -to use in a directory different than the R sources, use 
#  'source("PATH_TO_SOURCES/gsi.R", chdir = TRUE)'

# CHANGES:
# (05/01/2024 carlos) corrected symb.target.fstat.pval for the case when there is a single value and NAs;
# (28/07/2023 carlos) corrected name of avg.attr.correlation to 
# avg.abs.attr.correlation; completed minimal set of classification 
# measures; created the statlog-type and complete sets of GSI measures
# for regression; corrected bug in attr.correlation when dealing with
# constant attributes; target.sparsity renamed target.cv.sparsity;
# added gaussian kernel measures; completed set of statlog measures
# for classification; renamed class.freq to class.rel.freq;
# reimplemented measures concerning canonical correlation analysis;
# corrected class.covar to be able to handle empty matrices; additional corrections which should be further tested 
# (22/05/02 carlos) implemented regression measures (making use of 
# adaptations to data_characterization.R structure to support METAL 
# regression data format), general and specific for MLP parameter tuning;
# implemented correlation between numeric variables; reorganized
# measures
# (19/03/02 carlos) data set files can be in directory different than R sources
# (15/03/02 carlos) v1.0 changed structure to simplify extensibility: reading 
# data set and characterizing it are now separate (actually in different 
# files) and the latter uses the former; the latter is implemented as a 
# simple API, which is used here to implement GSI data characterization
# (30/11/01 carlos) corrected bug in output of class frequencies

# REGRESSION TODO:

# CLASSIFICATION TODO:
# -complete 7 basic measures
# -complete classifications measures

# KERNEL TODO:
# -could use outer product/lower triangular to calculate ideal matrix

# GENERAL TODO:
# -there should be a simple way to obtain the independent variables.
# -improve handling of symbolic values in the names file that do not appear in the data
# -improve handling of attributes with 0 standard-deviation (correlation stuff and dispersion gain)
# -does changing the names copy the data? (see r.squared and dispersion.gain)
# -a function to define the formula could be useful and another one to get the names
# of the attributes (in data_characterization.R"?)
# -handling of missing values should be carefully considered in attr.correlation and also the entropy based stuff
# -check for file handling errors

# EXAMPLE
# in R, do:
# > METALCharacterization("abalone", "abalone.gsi", kCompleteClassificationGSI)
#
# This small example assumes the existence of files "abalone.data" and
# "abalone.names", that should conform METAL syntax. It will save in file
# "abalone.gsi" a list of the basic set of data characteristics given in 
# kMinimalClassificationGSI. The measures are anonymized and written in a
# format suitable for J. Petrak's parse_results script.
# ######################################################################

# ======================================================================
# 
# This file implements the following measures:

### GENERAL ###

# N.EXAMPLES
# 	number of examples.
# LOG.N.EXAMPLES
#       log(10) of the number of examples.
# N.ATTRS
# 	number of attributes.     
# N.EXAMPLES.REL.N.ATTRS
# 	ratio of number of examples by number of attributes.
# LOG.N.EXAMPLES.REL.N.ATTRS
#       log of the ratio of number of examples by number of attributes.
# N.CONTINUOUS.ATTRS
#	number of continuous attributes.
# N.SYMBOLIC.ATTRS
#	number of symbolic attributes.     
# N.BINARY.ATTRS
#       number of binary attributes.
# PROP.CONTINUOUS.ATTRS
#	proportion of continuous attributes.
# PROP.SYMBOLIC.ATTRS
#	proportion of symbolic attributes.
# PROP.BINARY.ATTRS
#	proportion of binary attributes.
# N.MISSING.VALUES
#	number of missing values.     
# PROP.MISSING.VALUES
#	proportion of missing values.
# SYMB.VALUE.FREQ
#	for each k value of each j symbolic attribute, the value frequency.
#	It's stored in list indexed by each symbolic attribute and that each 
#	element contains itself a list indexed by each value of a specified 
#	symbolic attribute. 
# SYMB.PAIR.VALUE.FREQ
#	Frequency of pairs of values of symbolic variables.
# SYMB.PAIR.MUTUAL.INFORMATION
#	Mutual information between pairs of symbolic attributes.
# AVG.SYMB.PAIR.MUTUAL.INFORMATION
#	Average mutual information between pairs of symbolic attributes.
# RES.SYMB.PAIR.MUTUAL.INFORMATION
#	Rescaled mutual information between pairs of symbolic attributes.
# AVG.RES.SYMB.PAIR.MUTUAL.INFORMATION
#	Average rescaled mutual information between pairs of symbolic attributes.
# MIN.RES.SYMB.PAIR.MUTUAL.INFORMATION
#	Minimum rescaled mutual information between pairs of symbolic attributes.
# MAX.RES.SYMB.PAIR.MUTUAL.INFORMATION
#	Maximum rescaled mutual information between pairs of symbolic attributes.
# H.OUTLIER
#	for each continuous attribute the ratio between the standard deviation
#	 and the standard deviation of alpha trimmed mean. If the standard 
# 	deviation is 0, then the ratio is set to 1. It's stored in list 
#	indexed by each continuous attribute.   
# N.H.OUTLIER 
#	number of continuous attributes with outliers. 
# PROP.H.OUTLIER
#	proportion of continuous attributes with outliers.    
# ATTR.ENTROPY
#	for each symbolic attribute, the attribute entropy. It's stored in a 
#	list indexed by each symbolic attribute.    
# AVG.ATTR.ENTROPY
#	attribute entropy average.
# RES.ATTR.ENTROPY
#	for each symbolic attribute, the rescaled attribute entropy. It's
#       stored in a list indexed by each symbolic attribute.    
# AVG.RES.ATTR.ENTROPY
#	rescaled attribute entropy average.
# ATTR.CORRELATION
#	Correlation between continuous attributes.
# AVG.ABS.ATTR.CORRELATION
#	Average absolute correlation between continuous attributes.
# MIN.ABS.ATTR.CORRELATION
#	Minimum absolute correlation between continuous attributes.
# MAX.ABS.ATTR.CORRELATION
#	Maximum absolute correlation between continuous attributes.
# CLASS.ATTR.CORRELATION
#	Correlation between continuous attributes, per class.
# AVG.ABS.CLASS.ATTR.CORRELATION
#       Average absolute correlation between continuous attributes, per class.
# SKEWNESS
#       Skewness of each attribute.
# AVG.SKEWNESS
#       Mean skewness of attributes.
# AVG.ABS.SKEWNESS
#       Mean absolute skewness of attributes.
# CLASS.SKEWNESS
#       Skewness of each attribute, per class.
# AVG.CLASS.SKEWNESS
#       Mean class skewness of attributes.
# AVG.ABS.CLASS.SKEWNESS
#       Mean absolute class skewness of attributes.
# KURTOSIS
#       Kurtosis of each attribute.
# AVG.KURTOSIS
#       Mean kurtosis of attributes.
# CLASS.KURTOSIS
#       Kurtosis of each attribute, per class.
# AVG.CLASS.KURTOSIS
#       Mean class kurtosis of attributes.
# CLASS.COVAR
#       Class covariance matrices. Stored in a array (class is first dimension). 
#       Empty matrices will be assigned covariance of 0.
# POOLED.COVAR
#       Pooled covariance matrix.
# M.STAT
#       Box's M Statistic

### CLASSIFICATION ###

# N.CLASSES
#       number of class values.
# CLASS.ABS.FREQ
#	absolute class frequency (vector).
# CLASS.REL.FREQ
#	relative class frequency. Stored in a list indexed by each class value. 
# CLASS.ENTROPY
#	class entropy.
# CLASS.SYMB.VALUE.FREQ
#	for each k value of each j symbolic attribute and each i class, the 
# 	proportion of cases that have the k value in the j attribute and 
#	belong to the i class. It's stored in list indexed by each symbolic 
#	attribute and containing flat contingency tables that combine the 
#	values of the symbolic attribute with the class values.
# MUTUAL.INFORMATION
#	for each symbolic attribute, the mutual information between the 
#	attribute and the class. It's stored in a list indexed by each 
#	symbolic attributes.     
# AVG.MUTUAL.INFORMATION
#	mutual information average.    
# NOISE.SIGNAL.RATIO
#	noise signal ratio.
# TOTAL.SUM.SQUARES.CROSS.PRODUCTS
#	matrix of total sums of squares and cross products.
# WITHIN.GROUPS.SUM.SQUARES.CROSS.PRODUCTS
# 	matrix of within-groups sums of squares and cross products.
# BETWEEN.GROUPS.SUM.SQUARES.CROSS.PRODUCTS
#	difference between total and within-groups sums of squares and
#	cross products
# EIGEN.VALUES.LINEAR.DISCRIMINANT.FUNCTIONS
# 	Eigen values of linear discriminant functions.
# CANONICAL.CORRELATION.BEST.LINEAR.COMBINATION
#	Canonical correlation of the best linear combination of attributes
#	to distinguish between classes.
# RELATIVE.PROP.BEST.LINEAR.COMBINATION
#       Proportion of the total discrimination power explained by the
#       best linear combination.

### REGRESSION ###

# TARGET.COEFFICIENT.VARIATION
# 	ratio of the standard-deviation and the mean of the target 
#	attribute
# ABS.TARGET.COEFFICIENT.VARIATION
# 	absolute value of the ratio of the standard-deviation and
#       the mean of the target attribute
# TARGET.CV.SPARSITY
# 	0 if not sparse, 1 if sparse and 2 if extremely sparse, based 
#	on the coefficient of variation
# TARGET.ABSCV.SPARSITY
# 	0 if not sparse, 1 if sparse and 2 if extremely sparse, based 
#	on the absolute coefficient of variation
# TARGET.H.OUTLIER
# 	h.outlier value, as calculated for the continuous attributes
# TARGET.HAS.OUTLIERS
# 	1 if it has, 0 otherwise, based on the notion of outliers used 
#	for continuous attributes
# TARGET.STATIONARITY
# 	check is standard deviation is larger than mean
# R.SQUARED
# 	R^2 coeficient of multiple linear regression (only numerical 
#	attributes are used)
# R.SQUARED.BIN.SYMB
# 	R^2 coeficient of multiple linear regression (symbolic attributes 
#	are binarized)
# TARGET.CORRELATION
# 	correlation matrix between attributes and target
# AVG.ABS.TARGET.CORRELATION
# 	average correlation continuous attribute/target
# MIN.ABS.TARGET.CORRELATION
# 	minimum correlation continuous attribute/target
# MAX.ABS.TARGET.CORRELATION
# 	maximum correlation continuous attribute/target
# DISPERSION.GAIN
# 	gain in dispersion obtained with the best split for each attribute
# AVG.DISPERSION.GAIN
# 	average gain in dispersion for all attributes
# MIN.DISPERSION.GAIN
# 	minimum gain in dispersion for all attributes
# MAX.DISPERSION.GAIN
# 	maximum gain in dispersion for all attributes
# RES.SYMB.VALUE.TARGET.RANGE
#	for each k value of each j symbolic attribute, the range
#       of target values divided by the total target range.
# AVG.RES.SYMB.VALUE.TARGET.RANGE
#	average value of rescaled target range for each j symbolic attribute.
# MIN.AVG.RES.SYMB.VALUE.TARGET.RANGE
#	minimum of the average value of rescaled target range for each j
#       symbolic attribute.
# TARGET.HIST.SPARSITY
#       standard-deviation of the proportions of a histogram with 100 bins
#       of target values.
# SYMB.TARGET.FSTAT.PVAL
#       p value of the F statistic for each symbolic attribute and the target
# AVG.SYMB.TARGET.FSTAT.PVAL
#       average p value of the F statistic for each symbolic attribute
#       and the target
# MIN.SYMB.TARGET.FSTAT.PVAL
#       minimum p value of the F statistic for each symbolic attribute
#       and the target
# MAX.SYMB.TARGET.FSTAT.PVAL
#       maximum p value of the F statistic for each symbolic attribute
#       and the target
# MEAN.RES.DIST.ADJACENT.TARGET
#       mean distance between each target value and its two neighbors (sorted by
#       value.
# AVG.MEAN.RES.DIST.ADJACENT.TARGET
#       average mean distance between each target value and its two neighbors 
#       (sorted by value).
# MAX.MEAN.RES.DIST.ADJACENT.TARGET
#       maximum mean distance between each target value and its two neighbors 
#       (sorted by value=.


### OTHER ###
# MIN.DIST.CASE
# JAAKKOLA

### GAUSSIAN KERNEL ###
# GAUSSIAN.KERNEL.MEASURES

# ======================================================================

#source("/home/csoares/Investigacao/GSI/develop/gsi.R", chdir=TRUE)
#source("../gsi.R", chdir=TRUE)
#METALCharacterization("abalone", measures=kStatLogSelectedClass)
#ds <- ReadDataSet("adult_100", target.declared=F); dc <- list(); dc$value <- list()
#dc$value$n.examples <- n.examples(ds, NULL)


# LIBRARIES
source("data_characterization.R") # To enable working in a directory different than the one containing the scripts

library(rpart)
library(abind)

# PRE-DEFINED SETS OF MEASURES

kStatLogSelectedClass <- list(name = "gsi_statlog_selected_classification", version = "1.0", target.declared = FALSE, measures = c("n.examples", "n.attrs", "n.classes", "n.continuous.attrs", "n.symbolic.attrs", "prop.symbolic.attrs", "n.missing.values", "prop.missing.values", "symb.value.freq", "symb.pair.value.freq", "symb.pair.mutual.information", "avg.symb.pair.mutual.information", "attr.correlation", "avg.abs.attr.correlation", "class.rel.freq", "class.entropy", "class.symb.value.freq", "mutual.information", "avg.mutual.information", "attr.entropy", "avg.attr.entropy", "equivalent.n.attrs", "class.abs.freq", "class.covar", "total.sum.squares.cross.products", "within.groups.sum.squares.cross.products", "between.groups.sum.squares.cross.products", "eigen.values.linear.discriminant.functions", "canonical.correlation.best.linear.combination", "h.outlier", "n.h.outlier", "prop.h.outlier", "noise.signal.ratio"))

kStatLogClass <- list(name = "gsi_statlog_classification", version = "1.0", target.declared = FALSE, measures = c("n.examples", "n.attrs", "n.classes", "n.binary.attrs", "n.continuous.attrs", "class.abs.freq", "class.covar", "pooled.covar", "m.stat", "sd.ratio", "class.attr.correlation", "avg.abs.class.attr.correlation", "total.sum.squares.cross.products", "within.groups.sum.squares.cross.products", "between.groups.sum.squares.cross.products", "eigen.values.linear.discriminant.functions", "canonical.correlation.best.linear.combination", "relative.prop.best.linear.combination", "class.skewness", "avg.abs.class.skewness", "class.kurtosis", "avg.class.kurtosis", "class.rel.freq", "class.entropy", "symb.value.freq", "class.symb.value.freq", "mutual.information", "avg.mutual.information", "attr.entropy", "avg.attr.entropy", "equivalent.n.attrs", "noise.signal.ratio"))

kCompleteClassificationGSI <- list(name = "gsi_classification", version = "1.1", target.declared = FALSE, measures = c("n.examples", "n.attrs", "n.continuous.attrs", "n.symbolic.attrs", "prop.symbolic.attrs", "n.missing.values", "prop.missing.values", "class.rel.freq", "class.entropy", "symb.value.freq", "class.symb.value.freq", "mutual.information", "avg.mutual.information", "h.outlier", "n.h.outlier", "prop.h.outlier", "attr.entropy", "avg.attr.entropy", "noise.signal.ratio", "attr.correlation", "avg.abs.attr.correlation", "symb.pair.value.freq", "symb.pair.mutual.information", "avg.symb.pair.mutual.information", "total.sum.squares.cross.products", "class.covar", "class.abs.freq", "within.groups.sum.squares.cross.products", "between.groups.sum.squares.cross.products", "eigen.values.linear.discriminant.functions", "canonical.correlation.best.linear.combination"))

kContAttrsClassificationGSI <- list(name = "gsi_classification", version = "1.1", target.declared = FALSE, measures = c("n.examples", "n.attrs", "n.continuous.attrs", "n.missing.values", "prop.missing.values", "class.rel.freq", "class.entropy", "h.outlier", "n.h.outlier", "prop.h.outlier", "attr.correlation", "avg.abs.attr.correlation", "total.sum.squares.cross.products", "class.covar", "class.abs.freq", "within.groups.sum.squares.cross.products", "between.groups.sum.squares.cross.products", "eigen.values.linear.discriminant.functions", "canonical.correlation.best.linear.combination"))

kMinimalClassificationGSI <- kCompleteClassificationGSI

kStatLogRegressionGSI <- list(name="gsi_statlog_regression", version="1.0", target.declared=TRUE, measures=c("n.examples", "n.attrs", "n.binary.attrs", "target.coefficient.variation", "abs.target.coefficient.variation", "target.abscv.sparsity", "attr.correlation", "avg.abs.attr.correlation", "target.correlation", "max.abs.target.correlation", "skewness", "avg.abs.skewness", "kurtosis", "avg.kurtosis", "mean.res.dist.adjacent.target", "avg.mean.res.dist.adjacent.target", "symb.value.freq", "attr.entropy", "avg.attr.entropy", "res.attr.entropy", "avg.res.attr.entropy", "symb.target.fstat.pval", "avg.symb.target.fstat.pval"))

kSelectedRegressionGSI <- list(name="gsi_selected_regression", version="1.0", target.declared=TRUE, measures=c("n.examples", "log.n.examples", "n.attrs", "log.n.examples.rel.n.attrs", "n.continuous.attrs", "n.symbolic.attrs", "n.binary.attrs", "prop.continuous.attrs", "prop.symbolic.attrs", "prop.binary.attrs", "h.outlier", "n.h.outlier", "prop.h.outlier", "symb.value.freq", "res.attr.entropy", "avg.res.attr.entropy", "target.h.outlier", "mean.res.dist.adjacent.target", "avg.mean.res.dist.adjacent.target", "max.mean.res.dist.adjacent.target", "attr.correlation", "avg.abs.attr.correlation", "symb.pair.value.freq", "symb.pair.mutual.information", "attr.entropy", "res.symb.pair.mutual.information", "avg.res.symb.pair.mutual.information", "target.correlation", "max.abs.target.correlation", "symb.target.fstat.pval", "max.symb.target.fstat.pval", "avg.abs.target.correlation", "avg.symb.target.fstat.pval"))

kCompleteRegressionGSI <- list(name="gsi_complete_regression", version="1.0", target.declared=TRUE, measures=c("n.examples", "log.n.examples", "n.attrs", "n.examples.rel.n.attrs", "log.n.examples.rel.n.attrs", "n.continuous.attrs", "n.symbolic.attrs", "n.binary.attrs", "prop.continuous.attrs", "prop.symbolic.attrs", "prop.binary.attrs", "attr.correlation", "avg.abs.attr.correlation", "min.abs.attr.correlation", "max.abs.attr.correlation", "h.outlier", "n.h.outlier", "prop.h.outlier", "skewness", "avg.skewness", "kurtosis", "avg.kurtosis", "symb.value.freq", "symb.pair.value.freq", "symb.pair.mutual.information", "attr.entropy", "symb.pair.mutual.information", "res.symb.pair.mutual.information", "avg.res.symb.pair.mutual.information", "min.res.symb.pair.mutual.information", "max.res.symb.pair.mutual.information", "target.correlation", "avg.abs.target.correlation", "min.abs.target.correlation", "max.abs.target.correlation", "res.symb.value.target.range", "avg.res.symb.value.target.range", "min.avg.res.symb.value.target.range", "symb.target.fstat.pval", "avg.symb.target.fstat.pval", "min.symb.target.fstat.pval", "max.symb.target.fstat.pval", "target.stationarity", "target.coefficient.variation", "abs.target.coefficient.variation", "target.cv.sparsity", "target.abscv.sparsity", "target.hist.sparsity", "target.h.outlier", "target.has.outliers", "mean.res.dist.adjacent.target", "avg.mean.res.dist.adjacent.target", "r.squared", "r.squared.bin.symb", "dispersion.gain", "avg.dispersion.gain", "min.dispersion.gain", "max.dispersion.gain"))

## RETHINK (only landmarkers?)
#kLTRegressionGSI <- list(name = "gsi_lt_regression", version = "1.0", target.declared = TRUE, measures = c("n.examples", "n.attrs", "n.continuous.attrs", "n.symbolic.attrs", "prop.symbolic.attrs", "n.missing.values", "prop.missing.values", "h.outlier", "n.h.outlier", "prop.h.outlier", "target.coefficient.variation", "target.cv.sparsity", "target.h.outlier", "target.has.outliers", "target.stationarity", "r.squared", "r.squared.bin.symb", "target.correlation", "avg.abs.target.correlation", "dispersion.gain", "avg.dispersion.gain"))

kMinimalRegressionGSI <- kCompleteRegressionGSI

kMLPRegressionParTuningGSI <- list(name = "mlp_partuning_regression", version = "1.0", target.declared = TRUE, binarize=TRUE, measures = c("n.examples", "n.attrs", "n.examples.rel.n.attrs", "target.coefficient.variation", "r.squared", "attr.correlation", "avg.abs.attr.correlation"))

kJaakkola <- list(name="jaakkola_heuristic", version="1.0", target.declared = TRUE, binarize=TRUE, measures = c("min.dist.case", "jaakkola"))

kGaussianKernel <- list(name="gaussian_kernel_measures", version="1.0", target.declared = TRUE, binarize=TRUE, measures = c("gaussian.kernel.measures"))

# ======================================================================
METALCharacterization <-
	function(filestem,
		 output.file = NULL,
		 measures)

# IN:
# -filestem: filestem of the data files
# -output.file: name of the file to save the measures in (if NULL, they
#  will be written to stdout)
# -measures: list containing information about the data set and the measures to calculate (see data_characterization.R)

# ***
{
gsi <- CharacterizeDataSet(filestem, measures)

if (is.null(output.file))
	writeMETAL(gsi)
else
	{
	fd <- file(output.file, open="w")
	writeMETAL(gsi, fd = fd)
	close(fd)
	}
}

# ****************** DATA CHARACTERIZATION FUNCTIONS ******************

### GENERAL ###

# ***
n.examples <- function(dataset, data.char)
# ***
{
nrow(dataset$frame)
}

# ***
log.n.examples <- function(dataset, data.char)
# ***
{
log10(nrow(dataset$frame))
}

# ***
n.attrs <- function(dataset, data.char)
# ***
{
ncol(dataset$frame) - 1
}

# ***
n.examples.rel.n.attrs <- function(dataset, data.char)
# ***
{
nrow(dataset$frame) / (ncol(dataset$frame) - 1)
}

# ***
log.n.examples.rel.n.attrs <- function(dataset, data.char)
# ***
{
log10(nrow(dataset$frame) / (ncol(dataset$frame) - 1))
}

# ***
n.continuous.attrs <- function(dataset, data.char)

# -changed implementation
# ***
{
length(ContAttrs(dataset))
}

# ***
n.symbolic.attrs <- function(dataset, data.char)

# -changed implementation
# ***
{
length(SymbAttrs(dataset))
}

# ***
n.binary.attrs <- function(dataset, data.char)

# -subset of the symbolic attributes
# -this should be supported in data_characterization.R

# ***
{
sa <- SymbAttrs(dataset)
nvals <- sapply(dataset$attributes$attr.type[sa], length)
length(nvals[nvals == 2])
}

# ***
prop.continuous.attrs <- function(dataset, data.char)
# ***
{
GetMeasure("n.continuous.attrs", data.char) / GetMeasure("n.attrs", data.char)
}

# ***
prop.symbolic.attrs <- function(dataset, data.char)
# ***
{
GetMeasure("n.symbolic.attrs", data.char) / GetMeasure("n.attrs", data.char)
}

# ***
prop.binary.attrs <- function(dataset, data.char)
# ***
{
GetMeasure("n.binary.attrs", data.char) / GetMeasure("n.attrs", data.char)
}

# ***
n.missing.values <- function(dataset, data.char)
# ***
{
l <- unlist(as.list(dataset$frame),use.names=FALSE)
length(subset(l,is.na(l)))
}

# ***
prop.missing.values <- function(dataset, data.char)
# ***
{
n.missing.values <- GetMeasure("n.missing.values", data.char)
n.attrs <- GetMeasure("n.attrs", data.char)
n.examples <- GetMeasure("n.examples", data.char)

n.missing.values / (n.attrs * n.examples)
}

# ***
symb.value.freq <- function(dataset, data.char)

# -changed implementation
# ***
{
n.examples <- GetMeasure("n.examples", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
	value.freq <- list()
	for(wkAttr in wkSymbAttrs)
	{
		# Calculate frequencies of each attribute
		value.freq[[wkAttr]] <- table(dataset$frame[ , wkAttr])[dataset$attributes$attr.type[[wkAttr]]] / n.examples
	}
}
else
	value.freq <- NA

value.freq
}

# ***
h.outlier <- function(dataset, data.char)

# -changed implementation
# -needs parameter (alpha)
# ***
{
wkContAttrs <- ContAttrs(dataset)
if (length(wkContAttrs) > 0)
{
	alpha <- 0.05 # should be parameter?

	h.outl <- numeric(0)
	for(wkAttr in wkContAttrs) 
	{
		v <- as.numeric(subset(dataset$frame[,wkAttr],!is.na(dataset$frame[,wkAttr])))
		num.exs.alpha <- floor(length(v) * alpha)
	      	upper.bound <- length(v) - num.exs.alpha
		lower.bound <- num.exs.alpha+1	                
		v <- sort(v)[lower.bound:upper.bound]
		x.alpha <- mean(v)
					
		sd <- sd(dataset$frame[ , wkAttr], na.rm=TRUE)

		if(sd == 0)
			h.outl[[wkAttr]] <- 1
		else 
			h.outl[[wkAttr]] <- sqrt(sum((v-x.alpha)^2) / (length(v)-1)) / sd
	}
}
else
	h.outl <- NA

h.outl
}

# ***
n.h.outlier <- function(dataset, data.char)
# ***
{
h.outl <- GetMeasure("h.outlier", data.char)

if (! any(is.na(h.outl))) # TODO:(28/07/23 carlos) bug corrected?
#if (! is.na(h.outl))
		length(h.outl[h.outl < 0.7])
else
	NA
}

# ***
prop.h.outlier <- function(dataset, data.char)
# ***
{
n.continuous.attrs <- GetMeasure("n.continuous.attrs", data.char)
n.h.outl <- GetMeasure("n.h.outlier", data.char)

if (n.continuous.attrs == 0)
	NA
else
	n.h.outl / n.continuous.attrs
}

# ***
attr.entropy <- function(dataset, data.char)

# -changed implementation
# ***
{
symb.value.freq <- GetMeasure("symb.value.freq", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
	attr.entropy <- list()	
	for(wkAttr in wkSymbAttrs) 
	{
		attr.entropy[[wkAttr]] <- - sum(subset(symb.value.freq[[wkAttr]], symb.value.freq[[wkAttr]] != 0) * Log2(subset(symb.value.freq[[wkAttr]], symb.value.freq[[wkAttr]] != 0)))
	}
}
else
	attr.entropy <- NA

attr.entropy
}

# ***
avg.attr.entropy <- function(dataset, data.char)
# ***
{
mean(as.numeric(GetMeasure("attr.entropy", data.char)))
}

# ***
res.attr.entropy <- function(dataset, data.char)
{
symb.value.freq <- GetMeasure("symb.value.freq", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
	attr.entropy <- list()	
	for(wkAttr in wkSymbAttrs) 
	{
		attr.entropy[[wkAttr]] <- ( - sum(subset(symb.value.freq[[wkAttr]], symb.value.freq[[wkAttr]] != 0) * Log2(subset(symb.value.freq[[wkAttr]], symb.value.freq[[wkAttr]] != 0))) ) / Log2(length(subset(symb.value.freq[[wkAttr]], symb.value.freq[[wkAttr]] != 0)))
	}
}
else
	attr.entropy <- NA

attr.entropy
}

# ***
avg.res.attr.entropy <- function(dataset, data.char)
# ***
{
mean(as.numeric(GetMeasure("res.attr.entropy", data.char)))
}

# ***
class.attr.correlation <- function(dataset, data.char)
# ***
{
# Determine continuous attributes
wkContAttrs <- ContAttrs(dataset)

# Calculate correlation of numeric attributes for each class
if (length(wkContAttrs) > 0)
  do.call("abind", c(lapply(split(dataset$frame[, wkContAttrs], dataset$frame[,dataset$attributes$target.attr]), cor, use="pairwise.complete.obs"), along=0, use.anon.names=F))
else
  NA
}

# ***
avg.abs.class.attr.correlation <- function(dataset, data.char)
# ***
{
ac <- GetMeasure("class.attr.correlation", data.char)
mean(abs(ac[! is.na(ac)]))
}

# ***
attr.correlation <- function(dataset, data.char)
# ***
{
wkContAttrs <- ContAttrs(dataset)
attr.correlation <- NA
if (length(wkContAttrs) > 1)  # TODO:(28/07/23 carlos) bug corrected?
#if (length(wkContAttrs) > 0)  
	{
	cm <- cor(dataset$frame[, wkContAttrs], use="pairwise.complete.obs")
	xx <- rep(1:length(wkContAttrs), times=1:length(wkContAttrs)-1)
	yy <- unlist(sapply(2:length(wkContAttrs)-1, function(x) {1:x}))
	attr.correlation <- sapply(1:length(xx), function(x,m,xx,yy) {m[xx[x],yy[x]]}, cm, xx, yy)
	names(attr.correlation) <- sapply(1:length(xx), function(x, xx, yy) {paste(xx[x], yy[x], sep="_")}, xx, yy)
	}

attr.correlation[is.na(attr.correlation)] <- 0
attr.correlation
}

# ***
avg.abs.attr.correlation <- function(dataset, data.char)
# ***
{
ac <- GetMeasure("attr.correlation", data.char)
mean(abs(ac[! is.na(ac)]))
}

# ***
min.abs.attr.correlation <- function(dataset, data.char)
# ***
{
ac <- GetMeasure("attr.correlation", data.char)
min(abs(ac[! is.na(ac)]))
}

# ***
max.abs.attr.correlation <- function(dataset, data.char)
# ***
{
ac <- GetMeasure("attr.correlation", data.char)
max(abs(ac[! is.na(ac)]))
}


# **********************************************************************
symb.pair.value.freq <- function(dataset, data.char)
#
{
n.examples <- GetMeasure("n.examples", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 1)
{
	symb.value.pair.freq <- list()
	for(i in 1:(length(wkSymbAttrs) - 1))
	{
		# Identify first attribute values
		wkSymbAttrValues1 <- dataset$attributes$attr.type[[wkSymbAttrs[i]]]

		symb.value.pair.freq[[wkSymbAttrs[i]]] <- list()
		for (j in (i+1):length(wkSymbAttrs))
			{
			# Identify second attribute values
			wkSymbAttrValues2 <- dataset$attributes$attr.type[[wkSymbAttrs[j]]]

			# Calculate frequencies of each attribute value per class
			symb.value.pair.freq[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]] <- as.table(table(dataset$frame[,c(wkSymbAttrs[i], wkSymbAttrs[j])])[wkSymbAttrValues1, wkSymbAttrValues2]) / n.examples

			# Set correct names
			dimnames(symb.value.pair.freq[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]]) <- list(attr1 = wkSymbAttrValues1, attr2 = wkSymbAttrValues2)
		}
	}
}
else
	symb.value.pair.freq <- NA

symb.value.pair.freq
}

# **********************************************************************
symb.pair.mutual.information <- function(dataset, data.char)

# -checks frequencies (useful if used with samples)
#
{
svf <- GetMeasure("symb.value.freq", data.char)
spvf <- GetMeasure("symb.pair.value.freq", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 1)
{
	spmi <- list()
	for(i in 1:(length(wkSymbAttrs) - 1))
	{
		spmi[[wkSymbAttrs[i]]] <- list()
		for (j in (i+1):length(wkSymbAttrs))
			{
			selected <- (as.vector(spvf[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]]) != 0)
			spmi[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]] <- sum(as.vector(spvf[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]])[selected] * as.vector(Log2(as.vector(spvf[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]])[selected] / (rep(svf[[wkSymbAttrs[i]]], ncol(spvf[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]]))[selected] * rep(svf[[wkSymbAttrs[j]]], each=nrow(spvf[[wkSymbAttrs[i]]][[wkSymbAttrs[j]]]))[selected]))))
			}
	}
}
else
	spmi <- NA

spmi
}

# **********************************************************************
avg.symb.pair.mutual.information <- function(dataset, data.char)
#
{
spmi <- GetMeasure("symb.pair.mutual.information", data.char)
mean(unlist(spmi))
}

# **********************************************************************
res.symb.pair.mutual.information <- function(dataset, data.char)
#
{
spmi <- GetMeasure("symb.pair.mutual.information", data.char) #TODO: this doesn't look good

rspmi <- NA
if (! any(is.na(spmi)))  # TODO:(28/07/23 carlos) bug corrected?
#if (! is.na(spmi)) 		
  {
    attr.entr <- GetMeasure("attr.entropy", data.char)
    
    rspmi <- sapply(names(spmi), function(a1, mi, ae)
           {
             as.list(sapply(names(mi[[a1]]), function(a2, a1, mi, ae)
                            {
                              mi[[a1]][[a2]] / min(ae[[a1]], ae[[a2]])
                            }, a1, mi, ae))
           }, spmi, attr.entr)
  }
rspmi
}

# **********************************************************************
avg.res.symb.pair.mutual.information <- function(dataset, data.char)
#
{
rspmi <- GetMeasure("res.symb.pair.mutual.information", data.char)
mean(unlist(rspmi))
}

# **********************************************************************
min.res.symb.pair.mutual.information <- function(dataset, data.char)
#
{
rspmi <- GetMeasure("res.symb.pair.mutual.information", data.char)
min(unlist(rspmi))
}

# **********************************************************************
max.res.symb.pair.mutual.information <- function(dataset, data.char)
#
{
rspmi <- GetMeasure("res.symb.pair.mutual.information", data.char)
max(unlist(rspmi))
}

# ***
skewness <- function(dataset, data.char)
# ***
{
	wkContAttrs <- ContAttrs(dataset)
	if (length(wkContAttrs) > 0)
	{
		if (length(wkContAttrs) > 1) # TODO:(28/07/23 carlos) bug corrected?
		{
			means <- apply(dataset$frame[, wkContAttrs], 2, mean, na.rm=TRUE) # TODO:(28/07/23 carlos) bug corrected?
			sds <- apply(dataset$frame[, wkContAttrs], 2, sd, na.rm=TRUE) # TODO:(28/07/23 carlos) bug corrected?
		}
		else # TODO:(28/07/23 carlos) bug corrected?
		{
			means <- mean(dataset$frame[, wkContAttrs], na.rm=TRUE)
			sds <- sd(dataset$frame[, wkContAttrs], na.rm=TRUE)
		}
		s <- sapply(1:length(wkContAttrs), function(i, ca, df, m, s)
				{
					sum((df[, ca[i]] - m[i])^3)/(s[i]^3)
				}, wkContAttrs, dataset$frame, means, sds)
	}
	else
		s <- NA
	
	s/nrow(dataset$frame)
}

# ***
avg.skewness <- function(dataset, data.char)
# ***
{
mean(GetMeasure("skewness", data.char))
}

# ***
avg.abs.skewness <- function(dataset, data.char)
# ***
{
mean(abs(GetMeasure("skewness", data.char)))

}

# ***
class.skewness <- function(dataset, data.char)
# ***
{
wkContAttrs <- ContAttrs(dataset)
if (length(wkContAttrs) > 0)
  sapply(split(dataset$frame[, wkContAttrs], dataset$frame[,dataset$attributes$target.attr]), function(df)
         {
           apply(df, 2, function(col)
                 {
                   sd.col <- sd(col, na.rm=TRUE)
                   if(! is.na(sd.col) & (sd.col != 0 ) )
                     mean((col - mean(col, na.rm=TRUE)) ^ 3 / (sd.col ^ 3), na.rm=TRUE)
                   else
                     NA
                 })
         })
else
  NA
}


# ***
avg.class.skewness <- function(dataset, data.char)
# ***
{
mean(GetMeasure("class.skewness", data.char), na.rm=T)
}

# ***
avg.abs.class.skewness <- function(dataset, data.char)
# ***
{
mean(abs(GetMeasure("class.skewness", data.char)), na.rm=T)
}

# ***
kurtosis <- function(dataset, data.char)
# ***
{
	wkContAttrs <- ContAttrs(dataset)
	if (length(wkContAttrs) > 0)
	{
		if (length(wkContAttrs) > 1) # TODO:(28/07/23 carlos) bug corrected?
		{
			means <- apply(dataset$frame[, wkContAttrs], 2, mean, na.rm=TRUE) # TODO:(28/07/23 carlos) bug corrected?
			sds <- apply(dataset$frame[, wkContAttrs], 2, sd, na.rm=TRUE) # TODO:(28/07/23 carlos) bug corrected?
		}
		else # TODO:(28/07/23 carlos) bug corrected?
		{
			means <- mean(dataset$frame[, wkContAttrs], na.rm=TRUE)
			sds <- sd(dataset$frame[, wkContAttrs], na.rm=TRUE)
		}
		
		k <- sapply(1:length(wkContAttrs), function(i, ca, df, m, s)
				{
					sum((df[, ca[i]] - m[i])^4)/(s[i]^4)
				}, wkContAttrs, dataset$frame, means, sds)
	}
	else
		k <- NA
	
	k/nrow(dataset$frame)
}

# ***
avg.kurtosis <- function(dataset, data.char)
# ***
{
mean(GetMeasure("kurtosis", data.char))
}

# ***
class.kurtosis <- function(dataset, data.char)
# ***
{
wkContAttrs <- ContAttrs(dataset)
if (length(wkContAttrs) > 0)
  sapply(split(dataset$frame[, wkContAttrs], dataset$frame[,dataset$attributes$target.attr]), function(df)
         {
           apply(df, 2, function(col)
                 {
                   sd.col <- sd(col, na.rm=TRUE)
                   if(! is.na(sd.col) & (sd.col != 0 ) )
                     mean((col - mean(col, na.rm=TRUE)) ^ 4 / (sd.col ^ 4), na.rm=TRUE)
                   else
                     NA
                 })
         })
else
  NA
}

# ***
avg.class.kurtosis <- function(dataset, data.char)
# ***
{
mean(GetMeasure("class.kurtosis", data.char), na.rm=T)
}


### CLASSIFICATION ###

# ***
n.classes <- function(dataset, data.char)

  # -assumes classification problem
# ***
{
length(dataset$attributes$attr.type[[dataset$attributes$target.attr]])
}

# ***
class.abs.freq <- function(dataset, data.char)
# ***
{
c(table(dataset$frame[,dataset$attributes$target.attr])[dataset$attributes$attr.type[[dataset$attributes$target.attr]]])
}

# ***
class.rel.freq <- function(dataset, data.char)

# -changed implementation
# ***
{
n.examples <- GetMeasure("n.examples", data.char)

class.count <- table(dataset$frame[,dataset$attributes$target.attr])[dataset$attributes$attr.type[[dataset$attributes$target.attr]]]

class.freq <- alist()
class.freq <- as.list(class.count / n.examples)

class.freq
}


# ***
class.entropy <- function(dataset, data.char)
# ***
{
class.rel.freq <- GetMeasure("class.rel.freq", data.char)

-sum(as.numeric(subset(class.rel.freq, class.rel.freq!=0)) * Log2(as.numeric(subset(class.rel.freq, class.rel.freq!=0))))
}

# ***
class.symb.value.freq <- function(dataset, data.char)

# -changed implementation
# ***
{
n.examples <- GetMeasure("n.examples", data.char)
#symb.value.freq <- GetMeasure("symb.value.freq", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
	# Identify class values
	wkClassValues <- dataset$attributes$attr.type[[dataset$attributes$target.attr]]

	class.symb.value.freq <- list()
	for(wkAttr in wkSymbAttrs) 
	{
		# Identify attribute values
		wkSymbAttrValues <- dataset$attributes$attr.type[[wkAttr]]

		# Calculate frequencies of each attribute value per class
		class.symb.value.freq[[wkAttr]] <- as.table(table(dataset$frame[,c(wkAttr, dataset$attributes$target.attr)])[wkSymbAttrValues, wkClassValues]) / n.examples

		# Set correct names
		dimnames(class.symb.value.freq[[wkAttr]]) <- list(attr = wkSymbAttrValues, class = wkClassValues)
	}
}
else
	class.symb.value.freq <- NA

class.symb.value.freq
}

# **********************************************************************
total.sum.squares.cross.products <- function(dataset, data.char)
  
  #- covariance matrix * (n-1)!
{
#cont.attrs <- ContAttrs(dataset)

#if (length(cont.attrs) > 0)
#	{
#	classes <- dataset$attributes$attr.type[[dataset$attributes$target.attr]]

#	diff.mean <- dataset$frame[, cont.attrs] - matrix(mean(dataset$frame[, cont.attrs], na.rm=TRUE), nrow=nrow(dataset$frame), ncol=length(cont.attrs), byrow=TRUE)
#	T <- unlist(diff.mean[1, ]) %*% t(unlist(diff.mean[1, ]))
#	for (r in 2:nrow(dataset$frame))
#		T <- unlist(diff.mean[r, ]) %*% t(unlist(diff.mean[r, ]))
#	}
#else
#	T <- NA
#dimnames(T) <- list(cont.attrs, cont.attrs)
#T
cont.attrs <- ContAttrs(dataset)
if (length(cont.attrs) > 0)
  (nrow(dataset$frame) - 1) * cov(dataset$frame[,cont.attrs], use="pairwise.complete.obs")
else
  NA
}

# **********************************************************************
within.groups.sum.squares.cross.products <- function(dataset, data.char)
  
  #-pooled covariance matrix * (n_i - 1)!
#
{
#cont.attrs <- ContAttrs(dataset)

#if (length(cont.attrs) > 0)
#	{
#	classes <- dataset$attributes$attr.type[[dataset$attributes$target.attr]]

#	mean.by.class <- matrix(unlist(sapply(classes, function(class, df, tgt) {mean(df[tgt == class, cont.attrs], na.rm=TRUE)}, simplify=FALSE, dataset$frame, dataset$frame[[dataset$attributes$target.attr]])), nrow=length(classes), byrow=TRUE)
#	diff.mean <- dataset$frame[ , cont.attrs] - mean.by.class[dataset$frame[[dataset$attributes$target.attr]], ]
#	W <- unlist(diff.mean[1, ]) %*% t(unlist(diff.mean[1, ]))
#	for (r in 2:nrow(dataset$frame))
#		W <- unlist(diff.mean[r, ]) %*% t(unlist(diff.mean[r, ]))
#	}
#else
#	W <- NA
#dimnames(W) <- list(cont.attrs, cont.attrs)
#W
class.covar <- GetMeasure("class.covar", data.char)
class.abs.freq <- GetMeasure("class.abs.freq", data.char)

# Calculate pooled covariance matrix
if (length(ContAttrs(dataset)) > 0)
  apply(class.covar * (class.abs.freq - 1), 2:3, sum, na.rm=TRUE)
else
  NA
}

# **********************************************************************
between.groups.sum.squares.cross.products <- function(dataset, data.char)
#
{
GetMeasure("total.sum.squares.cross.products", data.char) - GetMeasure("within.groups.sum.squares.cross.products", data.char)
}

# **********************************************************************
eigen.values.linear.discriminant.functions <- function(dataset, data.char)

  #-converted to real because in some cases it may yield imaginary numbers
#
{
B <- GetMeasure("between.groups.sum.squares.cross.products", data.char)
W <- GetMeasure("within.groups.sum.squares.cross.products", data.char)

if (length(ContAttrs(dataset)) > 0)
  if (det(W, tol=1e-1000) != 0)
    {
      IW <- solve(W, tol=1e-1000)
      as.real(eigen(IW %*% B, only.values=TRUE)$values)
    }
  else
    NA
else
  NA
}

# **********************************************************************
canonical.correlation.best.linear.combination <- function(dataset, data.char)
  
  #-could be implemented with cancor(mva)
#
{
ev <- GetMeasure("eigen.values.linear.discriminant.functions", data.char) 
if (! is.na(ev))
  sqrt(ev[1] / (1 + ev[1]))
else
  NA
}

# **********************************************************************
relative.prop.best.linear.combination <- function(dataset, data.char)
#
{
ev <- GetMeasure("eigen.values.linear.discriminant.functions", data.char)
fract1 <- NA
if (! is.na(ev))
	fract1 <- ev[1]/sum(ev[ev>0])
fract1
}

# ***
mutual.information <- function(dataset, data.char)

# -changed implementation
# -checks frequencies (useful if used with samples)
# ***
{
class.rel.freq <- GetMeasure("class.rel.freq", data.char)
symb.value.freq <- GetMeasure("symb.value.freq", data.char)
class.symb.value.freq <- GetMeasure("class.symb.value.freq", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
	# Identify class values
	wkClassValues <- dataset$attributes$attr.type[[dataset$attributes$target.attr]]

	mutual.information <- list()	
	for(wkAttr in wkSymbAttrs) 
	{
		tmp <- 0

		for(wkClass in wkClassValues)
		{
			if (class.rel.freq[[wkClass]] > 0)
				tmp <- tmp + sum(class.symb.value.freq[[wkAttr]][symb.value.freq[[wkAttr]] > 0, wkClass] * Log2(class.symb.value.freq[[wkAttr]][symb.value.freq[[wkAttr]] > 0, wkClass]/(class.rel.freq[[wkClass]] * symb.value.freq[[wkAttr]][symb.value.freq[[wkAttr]] > 0])))
		}
		mutual.information[[wkAttr]] <- tmp
	}
}
else
	mutual.information <- NA

mutual.information
}

# ***
avg.mutual.information <- function(dataset, data.char)
# ***
{
mutual.information <- GetMeasure("mutual.information", data.char)
mean(as.numeric(mutual.information))
}

# ***
noise.signal.ratio <- function(dataset, data.char)
# ***
{
avg.attr.entropy <- GetMeasure("avg.attr.entropy", data.char)
avg.mutual.information <- GetMeasure("avg.mutual.information", data.char)

if (! is.na(avg.mutual.information) & (avg.mutual.information != 0) )
  (avg.attr.entropy - avg.mutual.information) / avg.mutual.information
else
  NA
}

# ***
equivalent.n.attrs <- function(dataset, data.char)
# ***
{
class.entropy <- GetMeasure("class.entropy", data.char)
avg.mutual.information <- GetMeasure("avg.mutual.information", data.char)

if (! is.na(avg.mutual.information) & (avg.mutual.information != 0) )
  class.entropy / avg.mutual.information
else
  NA
}

# ***
class.covar <- function(dataset, data.char)

  # Elements may be NA!!
# ***
{
# Determine continuous attributes
wkContAttrs <- ContAttrs(dataset)

if (length(wkContAttrs)>0)
{
  # Calculate covariance matrices of numeric attributes for each class
#do.call("abind", c(lapply(split(dataset$frame[, wkContAttrs], dataset$frame[,dataset$attributes$target.attr], drop=T), cov, use="pairwise.complete.obs")), along=0, use.anon.names=F)
do.call("abind", c(lapply(split(dataset$frame[, wkContAttrs], dataset$frame[,dataset$attributes$target.attr]), function(m) 
						{
							if (nrow(m)==0) 
							{
								dummy <- matrix(0, ncol=ncol(m), nrow=ncol(m))
								colnames (dummy) <- colnames(m)
								rownames (dummy) <- colnames(m)
								dummy
							} 
							else 
							{
							cov(m, use="pairwise.complete.obs")
							}
						}), along=0, make.names=F))
}
else
{
  NA
}
}


# ***
class.covar.tmp <- function(dataset, data.char)

# ***
{
# Determine continuous attributes
wkContAttrs <- ContAttrs(dataset)

if (length(wkContAttrs)>0)
  # Calculate covariance matrices of numeric attributes for each class
  do.call("abind", c(lapply(split(dataset$frame[, wkContAttrs], dataset$frame[,dataset$attributes$target.attr]), function(m)
                            {
                              means <- apply(m, 2, mean, na.rm=TRUE)
                              (t(as.matrix(m)) %*% as.matrix(m)) / (nrow(m) - 1) - means %*% t(means)
                            }), along=0, use.anon.names=F))
else
  NA
}

# ***
pooled.covar <- function(dataset, data.char)
  
  #-assumes classes are first dimension of class.covar
# ***
{
n.examples <- GetMeasure("n.examples", data.char)
n.classes <- GetMeasure("n.classes", data.char)
class.covar <- GetMeasure("class.covar", data.char)
class.abs.freq <- GetMeasure("class.abs.freq", data.char)

# Calculate pooled covariance matrix
if (! is.na(class.covar))
  apply(class.covar * (class.abs.freq - 1), 2:3, sum, na.rm=TRUE) / (n.examples - n.classes)
else
  NA
}

# ***
m.stat <- function(dataset, data.char)
# ***
{
n <- GetMeasure("n.examples", data.char)
q <- GetMeasure("n.classes", data.char)
p <- GetMeasure("n.continuous.attrs", data.char)
class.abs.freq <- GetMeasure("class.abs.freq", data.char)
class.covar <- GetMeasure("class.covar", data.char)
pooled.covar <- GetMeasure("pooled.covar", data.char)

if (! is.na(pooled.covar))
  {
    lambda <- 1 - (2 * p ^ 2 + 3 * p - 1) / (6 * (p + 1) * (q - 1)) * (sum(1 / (class.abs.freq - 1)) - 1 / (n - q))
    
#    lambda * sum((class.abs.freq - 1) * log(apply(class.covar, 1, function(m1, m2)
#                                                  {
#                                                    det(m2 / m1)
#                                                  }, pooled.covar)))

#    lambda * sum((class.abs.freq - 1) * log(apply(class.covar, 1, det) / det(pooled.covar)))
    
    lambda * sum((class.abs.freq - 1) * log(apply(class.covar, 1, function(m1, m2)
                                                  {
# a la DCT
#                                                    det(solve(m1) * m2)
                                                    if ((! any(apply(m1, 1, function(v) {all(is.na(v))}))) & (det(m1, tol=1e-1000) != 0))
                                                      det(solve(m1, tol=1e-1000) %*% m2, tol=1e-1000)
                                                    else
                                                      NA
                                                  }, pooled.covar)))
  }
else
  NA
}

# ***
sd.ratio <- function(dataset, data.char)
# ***
{
class.abs.freq <- GetMeasure("class.abs.freq", data.char)
p <- GetMeasure("n.continuous.attrs", data.char)
m.stat <- GetMeasure("m.stat", data.char)

if (! is.na(m.stat))
  exp(m.stat / (p * sum(class.abs.freq - 1)))
else
  NA
}


### REGRESSION ###

# ***
target.coefficient.variation <- function(dataset, data.char)
# ***
{
tgt.mean <- mean(dataset$frame[,dataset$attributes$target.attr])
tgt.sd <- sd(dataset$frame[,dataset$attributes$target.attr])

tgt.sd / tgt.mean
}

# ***
abs.target.coefficient.variation <- function(dataset, data.char)
# ***
{
tgt.cv <- GetMeasure("target.coefficient.variation", data.char)
abs(tgt.cv)
}

# ***
target.cv.sparsity <- function(dataset, data.char)
# ***
{
tgt.cv <- GetMeasure("target.coefficient.variation", data.char)
sparsity <- 2

if (tgt.cv < 0.2)
	sparsity <- 0
else if (tgt.cv < 0.5)
	sparsity <- 1
sparsity
}

# ***
target.abscv.sparsity <- function(dataset, data.char)
# ***
{
tgt.cv <- GetMeasure("abs.target.coefficient.variation", data.char)
sparsity <- 2

if (tgt.cv < 0.2)
	sparsity <- 0
else if (tgt.cv < 0.5)
	sparsity <- 1
sparsity
}

# ***
target.h.outlier <- function(dataset, data.char)
# ***
{
n.examples <- nrow(dataset$frame) # no need to call GetMeasure!

alpha <- 0.05 # should be parameter?

num.exs.alpha <- floor(n.examples * alpha)
upper.bound <- n.examples - num.exs.alpha
lower.bound <- num.exs.alpha + 1	                
v <- sort(dataset$frame[, dataset$attributes$target.attr])[lower.bound:upper.bound]
x.alpha <- mean(v)
sd <- sd(dataset$frame[ , dataset$attributes$target.attr])

if(sd == 0)
	h.outl <- 1
else 
	h.outl <- sqrt(sum((v-x.alpha)^2) / (n.examples - 1)) / sd

h.outl
}

# ***
target.has.outliers <- function(dataset, data.char)
# ***
{
h.outlier <- GetMeasure("target.h.outlier", data.char)

has.outliers <- 0
if (h.outlier < 0.7)
	has.outliers <- 1

has.outliers
}

# ***
target.stationarity <- function(dataset, data.char)
# ***
{
tgt.mean <- mean(dataset$frame[,dataset$attributes$target.attr])
tgt.sd <- sd(dataset$frame[,dataset$attributes$target.attr])

if (tgt.sd > tgt.mean)
	1
else
	0
}


# ***
r.squared <- function(dataset, data.char)

# -does changing the names implies copying the data frame?
# ***
{
rs <- NA
wkContAttrs <- ContAttrs(dataset)

if (length(wkContAttrs) > 0)
	{
	#summary(lm(dataset$frame[, dataset$attributes$target.attr] ~ matrix(apply(dataset$frame[, wkContAttrs], 2, as.double), ncol=length(wkContAttrs))))$r.squared # I think this copies the data frame

	names(dataset$frame) <- paste("v", names(dataset$frame), sep="")
	rs <- summary(lm(as.formula(paste("v", dataset$attributes$target.attr, " ~ ", paste(paste("v", wkContAttrs, sep=""), sep="", collapse="+"), sep="")), data=dataset$frame))$r.squared
	}
rs
}

# ***
r.squared.bin.symb <- function(dataset, data.char)
# ***
{
wkAttrs <- dataset$attributes$attr.name[dataset$attributes$attr.name != dataset$attributes$target.attr] # I don't know if -target.attr would work

summary(lm(dataset$frame[, dataset$attributes$target.attr] ~ matrix(apply(dataset$frame[, wkAttrs], 2, as.double), ncol=length(wkAttrs))))$r.squared
}

# ***
target.correlation <- function(dataset, data.char)
# ***
{
target.correlation <- NA

wkContAttrs <- ContAttrs(dataset)
if (length(wkContAttrs) > 0)
{
	target.correlation <- as.vector(cor(dataset$frame[, wkContAttrs], dataset$frame[,dataset$attributes$target.attr]))
	names(target.correlation) <- wkContAttrs
}

target.correlation
}

# ***
avg.abs.target.correlation <- function(dataset, data.char)
# ***
{
tc <- GetMeasure("target.correlation", data.char)
mean(abs(tc[! is.na(tc)]))
}

# ***
min.abs.target.correlation <- function(dataset, data.char)
# ***
{
tc <- GetMeasure("target.correlation", data.char)
min(abs(tc[! is.na(tc)]))
}

# ***
max.abs.target.correlation <- function(dataset, data.char)
# ***
{
tc <- GetMeasure("target.correlation", data.char)
max(abs(tc[! is.na(tc)]))
}

# ***
dispersion.gain <- function(dataset, data.char)

# -does changing the names implies copying the data frame?
# ***
{
wkAttrs <- dataset$attributes$attr.name[dataset$attributes$attr.name != dataset$attributes$target.attr]

names(dataset$frame) <- paste("v", names(dataset$frame), sep="")
rt <- rpart(as.formula(paste("v", dataset$attributes$target.attr, " ~ ", paste(paste("v", wkAttrs, sep=""), sep="", collapse="+"), sep="")), data=dataset$frame, control=rpart.control(maxdepth=1, maxcompete=length(wkAttrs)))
dispersion.gains <- rt$splits[rt$splits[,"adj"] == 0, "improve"][paste("v", wkAttrs, sep="")]
names(dispersion.gains) <- wkAttrs
dispersion.gains
}

# ***
avg.dispersion.gain <- function(dataset, data.char)
# ***
{
dg <- GetMeasure("dispersion.gain", data.char)
mean(dg[! is.na(dg)])
}

# ***
min.dispersion.gain <- function(dataset, data.char)
# ***
{
dg <- GetMeasure("dispersion.gain", data.char)
min(dg[! is.na(dg)])
}

# ***
max.dispersion.gain <- function(dataset, data.char)
# ***
{
dg <- GetMeasure("dispersion.gain", data.char)
max(dg[! is.na(dg)])
}

# ***
res.symb.value.target.range <- function(dataset, data.char)
# ***
{
wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
	value.target.range <- list()
	for(wkAttr in wkSymbAttrs)
	{
          value.target.range[[wkAttr]] <- tapply(dataset$frame[, dataset$attributes$target.attr], dataset$frame[, wkAttr], function(x)
                 {
                   max(x, na.rm=TRUE) - min(x, na.rm=TRUE)
                 })/(max(dataset$frame[, dataset$attributes$target.attr], na.rm=TRUE) - min(dataset$frame[, dataset$attributes$target.attr], na.rm=TRUE))
          value.target.range[[wkAttr]] [is.na(value.target.range[[wkAttr]])] <- 0
	}
}
else
	value.target.range <- NA

value.target.range
}

# ***
avg.res.symb.value.target.range <- function(dataset, data.char)
# ***
{
target.range <- GetMeasure("res.symb.value.target.range", data.char)

wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
  avg.target.range <-sapply(target.range, mean)
}
else
  avg.target.range <- NA

avg.target.range
}

# ***
min.avg.res.symb.value.target.range <- function(dataset, data.char)
# ***
{
avg.target.range <- GetMeasure("avg.res.symb.value.target.range", data.char)

min(avg.target.range)
}


# ***
target.hist.sparsity <- function(dataset, data.char)
# ***
{
sd(hist(breaks=100, dataset$frame[, dataset$attributes$target.attr],plot=FALSE)$counts/nrow(dataset$frame))
}

# ***
symb.target.fstat.pval <- function(dataset, data.char)
# ***
{
wkSymbAttrs <- SymbAttrs(dataset)
if (length(wkSymbAttrs) > 0)
{
  stfp <- sapply(wkSymbAttrs, function(s, df, t)
                 {
					 if (length(unique(df[!is.na(df[,s]),s])) > 1 )
                     {
                       summary(aov(t ~ df[, s], na.action="na.omit"))[[1]]$"Pr(>F)"[1]
                     }
                   else
                     {
                       1
                     }
                 }, dataset$frame,  dataset$frame[, dataset$attributes$target.attr])
}
else
  stfp <- NA

stfp
}

# ***
avg.symb.target.fstat.pval <- function(dataset, data.char)
# ***
{
mean(GetMeasure("symb.target.fstat.pval", data.char))
}

# ***
min.symb.target.fstat.pval <- function(dataset, data.char)
# ***
{
min(GetMeasure("symb.target.fstat.pval", data.char))
}

# ***
max.symb.target.fstat.pval <- function(dataset, data.char)
# ***
{
max(GetMeasure("symb.target.fstat.pval", data.char))
}

# ***
mean.res.dist.adjacent.target <- function(dataset, data.char)
# ***
{
st <- sort(dataset$frame[, dataset$attributes$target.attr])

dist <- c(st[2]-st[1], sapply(2:(nrow(dataset$frame)-1), function(i,st) { mean(st[i]-st[i-1], st[i+1]-st[i]) }, st), st[nrow(dataset$frame)]-st[nrow(dataset$frame)-1])
names(dist) <- as.character(1:nrow(dataset$frame))
dist / (st[length(st)] - st[1])
}

# ***
avg.mean.res.dist.adjacent.target <- function(dataset, data.char)
# ***
{
mean(GetMeasure("mean.res.dist.adjacent.target", data.char))
}

# ***
max.mean.res.dist.adjacent.target <- function(dataset, data.char)
# ***
{
max(GetMeasure("mean.res.dist.adjacent.target", data.char))
}

### OTHER ###

# ***
min.dist.case <- function(dataset, data.char)
# ***
{
library(mva)
d <- as.matrix(dist(dataset$frame[1:min(5000, nrow(dataset$frame)),ContAttrs(dataset)]))
diag(d) <- NA
apply(d, 2, min, na.rm=TRUE)
#wkContAttrs <- ContAttrs(dataset)
#wkSymbAttrs <- SymbAttrs(dataset)
#nas <- is.na(dataset$frame[, c(wkContAttrs, wkSymbAttrs)])
#ranges <- apply(dataset$frame[, c(wkContAttrs, wkSymbAttrs)], 2, function(c) {max(c, na.rm=TRUE) - min(c, na.rm=TRUE) })
#sapply(1:nrow(dataset$frame[, c(wkContAttrs, wkSymbAttrs)]), function(i, m, mna, mr) { min(sapply((1:nrow(m))[-i], function(j, i, m, mna, mr) { sqrt(sum((m[i,] - m[j,])^2, mr[xor(mna[i,], mna[j,])]^2, na.rm=TRUE)) }, i, m, mna, mr))}, dataset$frame[, c(wkContAttrs, wkSymbAttrs)], nas, ranges)
}

# ***
jaakkola <- function(dataset, data.char)
# ***
{
mean(GetMeasure("min.dist.case", data.char))
}



### KERNEL MATRIX ###

#***
gaussian.kernel.measures.old <- function(dataset, data.char)

  # -Assumes regression task
  # -Assumes all attributes are continuous
  # -Assumes there are no missing values
  # -Hardwired widths
  
#***
{
library(mva)
widths <- c(0.25, 1, 4, 16, 64, 256, 1000, 4000, 16000, 64000, 256000)

wkAttrs <- dataset$attributes$attr.name[dataset$attributes$attr.name != dataset$attributes$target.attr]

# Negative euclidean distances
ned <- lapply(1:(nrow(dataset$frame)-1), function(r, df)
            {
              l <- sapply((r+1):nrow(df), function(r2, r1, df)
                     {
                       - sum((df[r1, ] - df[r2, ]) ^ 2)
                     }, r, df)
              names(l) <- (r+1):nrow(df)
              l
            }, dataset$frame[, wkAttrs])
names(ned) <- 1:(nrow(dataset$frame)-1)

# Kernel "goodness" measures
mkv <- lapply(2*widths^2, function(d, ned)
             {
               mean(exp(unlist(ned) / d))
             }, ned)
names(mkv) <- widths
vkv <- lapply(2*widths^2, function(d, ned)
             {
               var(exp(unlist(ned) / d))
             }, ned)
names(vkv) <- widths

# Adjust target
target <- dataset$frame[ , dataset$attributes$target.attr] - mean(dataset$frame[ , dataset$attributes$target.attr])

# Ideal matrix
im <- lapply(1:(length(target) - 1), function(i, t)
            {
              l <- sapply((i+1):length(t), function(i2, i1, df)
                     {
                       t[i1] * t[i2]
                     }, i, t)
              names(l) <- (i+1):length(t)
              l
            }, target)
names(im) <- 1:(length(target) - 1)
#ted <- target %*% t(target)
#dimnames(ted) <- list(rows=1:nrow(ted), cols=1:ncol(ted))

# Kernel alignment to target
kat <- as.list(sapply(2*widths^2, function(d, ned, im, t)
              {
                2 * sum(exp(unlist(ned) / d) * unlist(im)) + sum(t * t)
              }, ned, im, target)
      /
        sapply(2*widths^2, function(d, ned, im, t)
               {
                 sqrt((2 * sum(exp(unlist(ned) / d) * exp(unlist(ned) / d)) + length(t)) * (2 * sum( unlist(im) * unlist(im)) + sum(t * t)))
               }, ned, im, target))
names(kat) <- widths

list(mean.kernel.value=mkv, var.kernel.value=vkv, kat=kat)
}



#***
gaussian.kernel.measures <- function(dataset, data.char)

  # -Assumes regression task
  # -Assumes all attributes are continuous
  # -Assumes there are no missing values
  # -Hardwired widths
  # -First 5000 cases only.
  
#***
{
widths <- c(0.25, 1, 4, 16, 64, 256, 1000, 4000, 16000, 64000, 256000)

wkAttrs <- dataset$attributes$attr.name[dataset$attributes$attr.name != dataset$attributes$target.attr]

# Negative, squared euclidean distances
ned <- -dist(dataset$frame[1:min(5000, nrow(dataset$frame)), wkAttrs])^2

# Kernel "goodness" measures
mkv <- lapply(2*widths^2, function(d, ned)
             {
               mean(exp(ned / d))
             }, ned)
names(mkv) <- widths
vkv <- lapply(2*widths^2, function(d, ned)
             {
               var(exp(ned / d))
             }, ned)
names(vkv) <- widths

# Adjust target
target <- dataset$frame[1:min(5000, nrow(dataset$frame)), dataset$attributes$target.attr] - mean(dataset$frame[1:min(5000, nrow(dataset$frame)), dataset$attributes$target.attr])

# Ideal matrix
im.tmp <- target %*% t(target)
im <- im.tmp[lower.tri(im.tmp)]
rm(im.tmp)

# Kernel alignment to target
kat <- as.list(sapply(2*widths^2, function(d, ned, im, t)
              {
                2 * sum(exp(ned / d) * im) + sum(t * t)
              }, ned, im, target)
      /
        sapply(2*widths^2, function(d, ned, im, t)
               {
                 sqrt((2 * sum(exp(ned / d) * exp(ned / d)) + length(t)) * (2 * sum( im * im) + sum(t * t)))
               }, ned, im, target))
names(kat) <- widths

list(mean.kernel.value=mkv, var.kernel.value=vkv, kat=kat)
}

# ********************* OTHER AUXILIARY FUNCTIONS *********************

Log2 <- function(values) sapply(values,log0) 

log0 <- function(x) 
{
	if (x == 0.0 || is.na(x)) 

		0.0
	else
		log2(x)
}

gaussian.dist <- function(x,y,w)
{
  exp(-sum((x-y)^2)/(2*w^2))
}



