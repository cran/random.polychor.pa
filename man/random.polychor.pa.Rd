\name{random.polychor.pa}
\alias{random.polychor.pa}
\title{A Parallel Analysis With Randomly Generated Polychoric
  Correlation Matrices
}
\description{The function perform a parallel analysis using simulated polychoric
  correlation matrices. The function will extract the eigenvalues from
  each random generated polychoric correlation matrix and from the
  polychoric correlation matrix of real data. A plot comparing
  eigenvalues extracted from the specified real data and simulated data
  will help determine which of real eigenvalue outperform random data.
}
\usage{random.polychor.pa(nvar, n.ss, nrep, nstep,
                          data.matrix, q.eigen, r.seed = "NULL")
}

\arguments{
  \item{nvar}{Number of variables (items) in the raw data matrix
}
  \item{n.ss}{
    Number of participants of the raw data matrix
}
  \item{nrep}{
Number of random samples that should simulated
}
  \item{nstep}{
Number of ordered categories of the generic item (e.g., Likert-like 3
  ordered category item). The function does not allow to vary the number
  of steps within the pool of items. So it is assumed that all the items
  have the same number of ordered category.
}
  \item{data.matrix}{
the name of raw data matrix
}
  \item{q.eigen}{
a number comprised within the interval of 0 and 1 and
indicating the quantile that is used to choose the number of non-random factors
}
  \item{r.seed}{
eventually, a preferred number that will be used to initialize the
random generator. Default value: 1335031435
}
}
\details{
The function perform a parallel analysis (Horn, 1976) using randomly
  simulated polychoric correlations. Generates nrep random samples of
  simulated data having specified the number of participants and of
  variables. The number of ordered categories in the generic item should
  be declared (i.e., two-categories, three-categories, etc.). This
  information is used to generate the nrep random raw datasets with
  the same characteristics of the original real data set. With the
  exclusion of nseed, all the remaining informations should be provided
  and specific checks are implemented to verify the effective provision
  and their coherence. A check for missing values within the real
  dataset is present and if present will be treated LISTWISE. In this
  case a warning message will prompt the user singalling how NA were
  treated (LISTWISE is by now the only treatment considered) and the new
  sample size. No further checks are made on the raw data. So
  out-of-range values are not detected and it is on the behalf of the
  user to make a preliminary check on the reliability of data. The
  function will extract the eigenvalues from each randomly generated
  polychoric matrices and the requested percentile is
  returned. Eigenvalues from polychoric correlation matrix obtained from
  real data is also compute and compared, in a (scree) plot, with the
  eigenvalues extracted from the simulation (Polychoric
  matrices). Recently, Cho, Li & Bandalos (2009) showed that, in using
  PA method, it is important to match the type of the correlation matrix
  used to recover the eigenvalues from real data with the type of
  correlation matrix used to estimate random eigenvalues. Crossing the
  type of correlations (using Polychoric correlation matrix to estimate
  real eigenvalues and random simulated Pearson correlation matrice) may
  result in a wrong decision (i.e., retaining more non-random factors
  than the needed). A comparison with eigenvalues extracted from both
  randomly simulated pearson correlation matrix and real data is also
  included. Finally, for both type of correlation matrix (Polychoric vs
  Pearson), the two versions (the classic squared coefficient and the
  4th power coefficient) of Velicer's MAP criterion are calculated
  (Velicer, 1976; Velicer, Eaton, & Fava, 2000) by implementing under R
  the code released by O'Connor (2000) for SPSS, SAS and MATLAB.
}
\value{
  The function returns the number of factors for
  Polychoric and Pearson Correlation PA methods ans also the number of
  fatctors chosen by the two Velicer's MAP criterion (original and 4th
  power) for both Polychoric and Pearson correlation
  matrix. Furthermore, the function will return the (scree) plot of the
  eigenvalues for real (Polychoric vs Pearson correlation matrices) and
  simulated data (Polychoric vs Pearson correlation matrices). Finally
  the following LIST of matrices will be printed:
  \item{$MAP.selection}{Returns a matrix with five columns (variables) and
  with as many rows as the number of selected factors (by the Velicer's
  MAP method) plus 1: Factor (i.e., the number of factors);
  POLY.MAP.squared (classic, squared MAP coefficient calcutated on the
  polychoric correlation matrix); POLY.MAP.4th (the modern, 4th power,
  of the MAP coefficient calculated on the polychoric correlation
  matrix); CORR.MAP.squared (classic, squared MAP coefficient calcutated
  on the Pearson correlation matrix); CORR.MAP.4th (the modern, 4th
  power, of the MAP coefficient calculated on the Pearson correlation
  matrix)}
  \item{$POLYCHORIC}{Returns a matrix with five columns (variables) and as
  many rows as the number of selected factors (by the Polychoric PA
  method) plus 1: Factors (number of factors); Emp.Polyc.Eigen
  (eigenvalues extracted from the empirical polychoric correlation
  matrix through the corFA function of nFactors package, i.e. by
  substituting the item communalities along the main diagonal of the
  correlation matrix); P.SimMeanEigen (the average n-th eigenvalue, extracted
  from Polycoric correlation matrix, of the nrep simulated random
  samples); P.SimSDEigen (the standard deviation for the n-th
  eigenvalue, extracted from Polycoric correlation matrix, of the nrep
  simulated random samples); P.SimQuant (the q.eigen*100 Percentile of
  the distribution of eigenvalues, extracted from the Polychoric
  correlation matrix, of the nrep simulated random samples)}
  \item{$PEARSON}{Returns a matrix with five columns (variables) and as
  many rows as the number of selected factors (by the Pearson correlation PA
  method) plus 1: Factors (number of factors); Emp.Pears.Eigen
  (eigenvalues extracted from the empirical Pearson correlation
  matrix through the corFA function of nFactors package, i.e. by
  substituting the item communalities along the main diagonal of the
  correlation matrix); C.SimMeanEigen (the average n-th eigenvalue, extracted
  from Pearson correlation matrix, of the nrep simulated random
  samples); C.SimSDEigen (the standard deviation for the n-th
  eigenvalue, extracted from Pearson correlation matrix, of the nrep
  simulated random samples); C.SimQuant (the q.eigen*100 Quantile
  eigenvalue, extracted from the Pearson correlation matrix, of the
  nrep simulated random samples)}
}
\references{
Cho, S.J., Li, F., & Bandalos, D., (2009). Accuracy of the Parallel
  Analysis Procedure With Polychoric Correlations. Educational and Psychological
  Measurement, 69, 748-759.

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. Psychometrika, 32, 179-185.

O\kbd{'}Connor, B. P. (2000). SPSS and SAS programs for determining the number
  of components using parallel analysis and Velicer\kbd{'}s MAP test. Behavior
  Research Methods, Instrumentation, and Computers, 32, 396-402.

Velicer, W. F. (1976). Determining the number of factors from the matrix
of partial correlations. Psychometrika, 41, 321-327.

Velicer, W. F., Eaton, C. A., & Fava, J. L. (2000). Construct
  explication through factor or component analysis: A review and
  evaluation of alternative procedures for determining the number of
  factors or components. In R. D. Goffin & E. Helmes (Eds.), Problems
  and solutions in human assessment: Honoring Douglas N. Jackson at
  seventy (pp. 41-72). Norwell, MA: Kluwer Academic.
}
\author{
  Fabio Presaghi \email{fabio.presaghi@uniroma1.it}  and

Marta Desimoni \email{marta.desimoni@uniroma1.it}
}
\note{
In running the random.polychor.pa function it should be reminded that it
  may take a lot of time to complete the simulation. This is due in part
  to the fact that the estimation of the polychoric correlation matrix
  is cumbersome and in part to the fact that the code is not optimized,
  but simply it does the work.

  Occasionally, in calculating the polychoric correlation matrix it may
  occur an error when the matrix is non-positive definite. In this case
  you have to re-run the simulation.
}

\seealso{
nFactors, psych, paran.
}
\examples{
### example data
raw.data<-data.frame(ss=1:10, v1=c(1,3,2,3,2,2,2,3,2,3),
                              v2=c(2,3,3,3,3,3,1,2,3,3),
                              v3=c(2,3,2,3,3,2,1,3,2,3))
raw.item.data <- (raw.data[,2:4])
summary (raw.item.data)
cor(raw.item.data)
eigen(cor(raw.item.data))

random.polychor.pa(nvar=3, n.ss=10, nrep=5, nstep=3,
                   data.matrix=raw.item.data, q.eigen=.99)



######## for SPSS users ####
#### the following instructions can used to load a SPSS data file (.sav).
#### Prepare a spss datafile with only the items needed to perform the
#### Parallel Analysis and then run R
#### 1) within R, load the needed library to import spss datafile
#### 2) choose the SPSS datafile by pointing directly in the folder on your
#### hard-disk
#### 3) check if the imported data is correct 
#
# library(foreign) ### load the needed library
# raw.data <- read.spss(choose.files(), use.value.labels=TRUE,
#                       max.value.labels=Inf, to.data.frame=TRUE)
# summary (raw.data)

}
\keyword{ PARALLEL ANALISYS }
\keyword{ POLYCHORIC CORRELATION }
