pkgname <- "random.polychor.pa"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('random.polychor.pa')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("random.polychor.pa")
### * random.polychor.pa

flush(stderr()); flush(stdout())

### Name: random.polychor.pa
### Title: A Parallel Analysis With Randomly Generated Polychoric
###   Correlation Matrices
### Aliases: random.polychor.pa
### Keywords: PARALLEL ANALISYS POLYCHORIC CORRELATION EXPLORATORY FACTOR
###   ANALYSIS

### ** Examples

### EXAMPLE 1:
### example data
raw.data<-data.frame(ss=1:10, v1=c(1,2,2,1,2,2,2,1,2,1), 
                              v2=c(2,3,3,2,3,3,1,2,3,1), 
                              v3=c(2,4,2,3,3,2,1,3,2,4),
                              v4=c(3,1,3,2,3,2,3,2,2,3),
                              v5=c(3,1,4,5,3,4,3,4,2,5))
raw.item.data <- (raw.data[,2:6])
summary (raw.item.data)
cor(raw.item.data)
eigen(cor(raw.item.data))

random.polychor.pa(nrep=5, data.matrix=raw.item.data, q.eigen=.99)


### EXAMPLE 2a:
### this example is particularly instructive on how the solution may
### change by changing the type of correlation, method of extraction and
### method of selection.
### Before launching the example consider that the
### ESTIMATED TIME TO COMPLETE THE SIMULATION IS ABOUT: 10 MIN.
require(psych)
data(bfi)
raw.data<-as.matrix(bfi)
raw.data <- (raw.data[1:100,2:6])
test.1<-random.polychor.pa(nrep=3, data.matrix=raw.data, q.eigen=.99)
test.1

### EXAMPLE 2b:
### in this example one of the categories of item1 is recoded: 2=1
### so this item has 5 categories: 1 (2) 3 4 5 6
### category 1 is within brackets as it has frequency=0
### so this is a case where empirical data (0 2 3 4 5 6) diverge from
### theorethical data (0 1 2 3 4 5 6)
require(psych)
data(bfi)
raw.data.1<-as.matrix(bfi)
raw.data.1 <- (raw.data.1[1:100,1:25])
for(i in 1:nrow(raw.data.1)) { if(raw.data.1[i,1]==2) raw.data.1[i,1]<-1} 
test.2<-random.polychor.pa(nrep=3, data.matrix=raw.data.1, q.eigen=.99)
test.2

### EXAMPLE 2c:
### in this example one of the categories of item1 is recoded: 1=0
### so this item has one of its categories coded as 0 
### this will cause polychoric() function to stop with error
### and the random.polychoric.pa will prompt a warning message
#require(psych)
#data(bfi)
#raw.data.2<-as.matrix(bfi)
#raw.data.2 <- (raw.data.2[1:100,1:25])
#for(i in 1:nrow(raw.data.2)) { if(raw.data.2[i,1]==1) raw.data.2[i,1]<-0} # recode 1=0
# random.polychor.pa(nrep=3, data.matrix=raw.data.2, q.eigen=.99)

### EXAMPLE 3:
######## for SPSS users ####
### the following instructions can used to load a SPSS data file (.sav).
### 1) load the library to read external datafile (e.g., SPSS datafile)
### 2) choose the SPSS datafile by pointing directly in the folder on your hard-disk
### 3) select only the variables (i.e., the items) needed to for Parallel Analysis
# library(foreign) ### load the needed library
# raw.data <- read.spss(choose.files(), use.value.labels=TRUE,
#                       max.value.labels=Inf, to.data.frame=TRUE)
# raw.spss.item <- na.exclude(raw.data[,2:4])
# summary (raw.spss.item)
# random.polychor.pa(nrep=5, data.matrix=raw.spss.item, q.eigen=.99)




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
