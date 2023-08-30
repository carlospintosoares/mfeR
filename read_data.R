#
# read_data.R
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
# This file contains R functions to read in a METAL dataset (data and names 
# files)
#
# ----------------------------------------------------------------------

# BY: Carlos Soares

# COMMENTS:
# -assumes two files to describe the data set
# -assumes data is in correct format (I think this currently means Luis
#  Torgo's format)

# CHANGES:
# (29/08/23 carlos) corrected warning in ReadDataFrame that, from v4.3.0, 
# is an error
# (27/07/23 carlos) corrected error in regexpr
# (07/05/04 carlos) v2.1: removed variables with underscores which are
# not allowed anymore (Rita's fault!! :-)
# (14/03/02 carlos) v2.0: changed structure for names file representation; 
# provided option to anonymize.
# (24/02/02 carlos) changed structure to simplify extensibility: reading 
# dataset and characterizing it are now separate (actually in different 
# files) and the latter uses the former.

# TODO:
# -names file can also be given separately (e.g. ignore)
# -is attributes$attr.name necessary? (== names(attr.type))
#
# To have a fast trial on it just do:
# shell> R
# > source("read_data.R")
# > df <- ReadDataSet("housing")
#
# This small example assumes the existance of files "housing.data" and
# "housing.names", that should conform METAL sintax (actually it accepts
# some liberal interpretations of this sintax, but as long as you stick
# to METAL conventions you should be safe).
# The function "ReadDataSet" will return a list with several elements
# containing the information about the data set and the data itself
# ######################################################################

# ======================================================================
ReadDataSet <- function(probName,
 			attrs.exts=c(".domain",".names"),
 			data.exts=c(".data"),
			target.declared=TRUE,
			data.file=NULL,
			anonymize = FALSE)

# This is the top level function that you should use to read in a dataset

# IN:
# -probName
# Filestem of files with the data.
#
# -attrs.exts
# This function will assume that the names file will have as extension 
# either ".names" or ".domain". If other extension is used we may specify
# it through the "attrs.exts" parameter.
# e.g. : > d <- ReadDataSet("housing",attrs.exts=".declarations")
#
# -data.exts
# It defaults to the extension ".data" but you may likewise specify other
# alternative using the same strategy as for the "attrs.exts" parameter.
#
# -target.declared
# By default it is assumed that the last declaration of the "names" file
# corresponds to the type declaration of the target variable. However, there
# are some formats where this target variable is not declared in the "names"
# file. For these cases you should explicitly use the following:
# > d <- ReadDataSet("housing",target.declared=FALSE)
#
# -data.file
# By default it is assumed that the "data" filename is the same of the "names" 
# file. However, other filename can be specified if you want to run a 
# different dataset for the same problem.
# e.g. : > d <- ReadDataSet("housing",data.file="housing2.data")
#
# -anonymize
# Replace class values, attribute names and attribute values
# with corresponding indexes in the names file

# OUT:
# a list containing components that describe the names (see ReadtAttrsInfo) 
# and the data (see ReadData) files

# ***

{
# Setup data structures
wkDataset <- alist()
class(wkDataset) <- "dataset" # why is this necessary?
  
# Read names file
wkProbDescription <- ReadAttrsInfo(probName,attrs.exts,target.declared)

# Anonymize, creating a translation table for data anonymization purposes
wkAnonymizationTable <- NULL
if (anonymize)
  {
    wkAnonymizationTable <- list()
    
    # Anonymize target attribute
    wkProbDescription$attributes$target.attr <- as.character(match(wkProbDescription$attributes$target.attr, wkProbDescription$attributes$attr.name))

    # Anonymize attribute names
    wkProbDescription$attributes$attr.name <- as.character(1:length(wkProbDescription$attributes$attr.name))
    wkProbDescription$attributes$original.name <- wkProbDescription$attributes$attr.name
    names(wkProbDescription$attributes$attr.type) <- wkProbDescription$attributes$attr.name

    # Anonymize symbolic attribute values
    for (wkAttr in wkProbDescription$attributes$attr.name)
      {
		  if (wkProbDescription$attributes$attr.type[[wkAttr]][1] != "continuous") # TODO:(28/07/23 carlos) bug corrected? 
		  #if (wkProbDescription$attributes$attr.type[[wkAttr]] != "continuous")
          {
            wkAnonymizationTable[[wkAttr]] <- as.character(1:length(wkProbDescription$attributes$attr.type[[wkAttr]]))
            names(wkAnonymizationTable[[wkAttr]]) <- as.character(wkProbDescription$attributes$attr.type[[wkAttr]])
            wkProbDescription$attributes$attr.type[[wkAttr]] <- as.character(1:length(wkProbDescription$attributes$attr.type[[wkAttr]]))
          }
      }
  }

# Read data file
# -maybe this will force  several copies of the data frame to be made :-(
c(wkProbDescription, ReadData(probName, wkProbDescription, data.file, data.exts, anonymization.table = wkAnonymizationTable))
}

# ======================================================================
ReadAttrsInfo <- function(fileName,
			attrs.exts=c(".domain",".names"),
			target.declared=TRUE) 

# Auxiliary function used to read in the "names" file information (i.e.
# the information on the attributes of the domain under consideration).

# IN:
# -fileName: filestem of the files containing the data
# -attrs.ext
# -target.declared: see ReadDataSet

# OUT:
# a list containing the name of the names file and another list containing:
# -target.attr
#	the name of the target attribute.
# -problem.type
#	classification or regression. 
# -attr.name
#	the attribute names assigned by R
# -attr.original.name
#	the original attribute names
# - "attr.type"
#	attribute type ("continuous" or list of symbolic values)

# ***

{
# Read .names file
attrsFile <- FindFile(fileName,attrs.exts,"variables")
   
lines <- readLines(attrsFile)

# Identify names and types of attributes	        
original.attr.name <- sapply(subset(lines,regexpr(":",lines)!=-1),
                             FUN="ExtractAttrName",USE.NAMES= FALSE)

attr.name <- make.names(original.attr.name)

attr.type <- sapply(subset(lines,regexpr(":",lines)!=-1),
                    FUN="ExtractAttrType",USE.NAMES=FALSE)

# analysing the first line of the domain file in order to determine 
# the target attribute and the problem type .
   
if(regexpr(":",lines[1]) == -1 && regexpr(",",lines[1]) != -1)
  {
    if(target.declared) 
      target.attr <- attr.name[length(attr.name)]
    else
      {
        target.attr <- "class"                  

        original.attr.name <- c(original.attr.name,"class")
        
        attr.name <- c(attr.name,"class")
        
#			attr.type <- c(attr.type, list(ExtractAttrType(lines[1]))) #why?
        attr.type <- c(attr.type, ExtractAttrType(lines[1]))
      }
             
    problem.type <- "classification"               

  }
else
  {

    if(regexpr(":",lines[1]) != -1)
      target.attr <- attr.name[length(attr.name)] 
    else 
      target.attr <- make.names(substr(lines[1],1,
                                       attr(regexpr("[^\\.]+",lines[1]),
                                            "match.length")))
    
    if(attr.type[
                 match(target.attr,attr.name)] == "continuous")
      problem.type <- "regression"
    else
      problem.type <- "classification"
  }

# Set names of type slots to attribute name
names(attr.type) <- attr.name
	
list(namesfile=attrsFile, attributes=list(target.attr = target.attr, problem.type = problem.type, attr.name = attr.name, original.name = original.attr.name, attr.type = attr.type))
}

# ======================================================================
ReadData <- function(probName,
			dataset,
			data.file=NULL,
			data.exts=c(".data"),
			anonymization.table)

# Auxiliary function used to read in the "data" file information, taking 
# the basic features of the problem obtained from the "names" file.

# IN:
# -dataset: 
# List containing the basic features of the problem returned by 
# the ReadAttrsInfo function.
#
# -anonymization.table:
# Vector with translation from real to anonymous values. If null, no
# translation is done
#
# -probName
# -data.exts
# -data.file: see ReadDataSet

# OUT:
# -list with the name of the data file and the data frame containing the data
# ***
{
if(is.null(data.file))
  {
    data.file <- FindFile(probName,data.exts,"data")
  }
else
  {
    fileName <- unlist(strsplit(data.file,"\\."))
    data.file <- FindFile(fileName[1],paste(".",fileName[2],
                                            sep=""),"data")
  }


frame <- ReadDataFrame(data.file, dataset$attributes, anonymization.table)

list(data.file = data.file, frame = frame)
}

# ************ OTHER AUXILIARY FUNCTIONS ************

FindFile <- function(fileName,fileExts,fileType) 
{
ext.ind <- match(TRUE,file.exists(paste(fileName,fileExts,sep="")))

if(is.na(ext.ind)) stop(paste(fileType,"file missing"))
             
paste(fileName,fileExts[ext.ind],sep="")
}

ExtractAttrName <- function(line)
  substr(line,1,regexpr("([ \t]*:)+",line)-1)


ExtractAttrType <- function(line) 
{
# get a more specific description of the attribute
attrDescr <- unlist(strsplit(line,":"))

attrType <- attrDescr[length(attrDescr)]

r <- regexpr("[^\t:.]+",attrType)

list(gsub("([ ])","",unlist(strsplit(
                                     substr(attrType,r,attr(r,"match.length")+r-1),","))))
}

GetSeparator <- function(string) 
{
if(regexpr("\"",string)==-1 && regexpr("\'",string)!=-1)  
  GetSeparatorR(string,"\'")
else 
  GetSeparatorR(string)
}

GetSeparatorR <- function(string,charStr="\"") 
{
if(string=="") return("")

ind <- regexpr(",",string)

# If there are no commas, then the separator would be space or tab.
if(ind == -1) return("")
             
# check if the comma is not within a string
if((length(unlist(strsplit(substr(string,1,ind),charStr)))%%2))
  return(",")

ind1 <- regexpr(charStr,substr(string,ind+1,nchar(string)))

GetSeparatorR(substr(string,ind+ind1+1,nchar(string)),charStr)
}

ReadDataFrame <- function(data.file, attr.info, anonymization.table) 
{
	
dataSep <- GetSeparator(scan(file=data.file,what=character(0),
                             nmax=1,,sep="\n",quiet=TRUE))
	
if (! is.null(anonymization.table))
  data <- read.table(data.file, check.names = FALSE, col.names=attr.info$attr.name, sep=dataSep, strip.white=TRUE, na.strings="?")
else
  data <- read.table(data.file, col.names=attr.info$attr.name, sep=dataSep, strip.white=TRUE, na.strings="?")

for(wkAttr in attr.info$attr.name) 
  {
#    if((attr.info$attr.type[[wkAttr]] != "continuous" && attr.info$attr.type[[wkAttr]] != "nominal")) # What is nominal? # NOTE: fails in R v4.3.1
	if (length(attr.info$attr.type[[wkAttr]]) > 1) # assumes "nominal" doesn't appear 
		{
        # Anonymize
        if (! is.null(anonymization.table))
          data[, wkAttr] <- anonymization.table[[wkAttr]][as.vector(data[, wkAttr])] # I don't understand the need for as.vector but it works :-)

	# Factors should include all values in the .names file (e.g. samples)
        f <- factor(data[, wkAttr], levels = attr.info$attr.type[[wkAttr]])
        data[ , wkAttr] <- f	
      }
  }
data
}

