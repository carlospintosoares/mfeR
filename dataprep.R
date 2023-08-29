# 
# Author: csoares
###############################################################################


x <- df$horsepower
n <- 4

# simple discretization
nbins.discretize <- function(x, n)
{
	domain <- range(x, na.rm=T)
	bin.size <- (domain[2] - domain[1]) / n
	less.than.upper.bin.limit <- sapply(1:n, function(i) # TODO: ugly code
			{
				x <= domain[1] + i * bin.size
			})
	bins <- apply(less.than.upper.bin.limit, 1, function(bins)
			{
				min((1:n)[bins])
			})
	new.x <- paste0("v",bins)
	new.x[is.na(x)] <- NA
	as.factor(new.x)
}



###
### AUTP-MPG
###

# NOTE: create .names files accordingly

setwd("auto+mpg")

# read data
df <- read.csv(file="auto-mpg.data-tmp", header=F, na.strings=c("NA", "?"))

# add column names
names(df) <- c("mpg", "cylinders", "displacement", "horsepower", "weight", "acceleration", "model.year", "origin")

# write file
write.table(df, file="auto-mpg-ncont-0nom.data", quote=F, sep=",", na="?", row.names=F, col.names=F)

# represent origin as nominal
new.origin.values <- c("a","b","c")
new.origin <- new.origin.values[df$origin]
df$origin <- as.factor(new.origin)

# write file
write.table(df, file="auto-mpg-ncont-1nom.data", quote=F, sep=",", na="?", row.names=F, col.names=F)

# represent horsepower and discplacement as nominal
new.horsepower <- nbins.discretize(df$horsepower, 4)
new.displacement <- nbins.discretize(df$displacement, 4)

df$horsepower <- new.horsepower
df$displacement <- new.displacement

# write file
write.table(df, file="auto-mpg-ncont-nnom.data", quote=F, sep=",", na="?", row.names=F, col.names=F)

# represent all but one variables as nominal
new.cylinders <- nbins.discretize(df$cylinders, 5)
new.weight <- nbins.discretize(df$weight, 2)
new.model.year <- nbins.discretize(df$model.year, 2)

df$cylinders <- new.cylinders
df$weight <- new.weight
df$model.year <- new.model.year

# write file
write.table(df, file="auto-mpg-1cont-nnom.data", quote=F, sep=",", na="?", row.names=F, col.names=F)

# represent remaining variable as nominal
new.acceleration <- nbins.discretize(df$acceleration, 3)

df$acceleration <- new.acceleration

# write file
write.table(df, file="auto-mpg-0cont-nnom.data", quote=F, sep=",", na="?", row.names=F, col.names=F)

