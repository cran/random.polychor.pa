\name{Parallel-Analysis-of-Polychoric-Correlations}
\alias{random.polychor.pa-package}
\docType{package}
\title{
  A Parallel Analysis with Random Polychoric Correlation Matrices
}
\description{The function performs a parallel analysis using simulated
  polychoric correlation matrices. The function will extract the
  eigenvalues from each random generated polychoric correlation matrix
  and from the polychoric correlation matrix of real data. A plot
  comparing eigenvalues extracted from the specified real data with
  simulated data will help determine which of real eigenvalue outperform
  random data. A series of matrices comparing MAP vs PA-Polychoric vs
  PA-Pearson correlations methods, FA vs PCA solutions are finally presented.
}
\details{
\tabular{ll}{
Package: \tab random.polychor.pa\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2010-04-27\cr
License: \tab GPL Version 2 or later\cr
LazyLoad: \tab yes\cr
}
The function perform a parallel analysis (Horn, 1965) using randomly
simulated polychoric correlations. Generates nrep random samples of the
same dimension of the empirical provided data.matrix(i.e, with the same
number of participants and of variables). Items are allowed to have
varying number of categories. The function will extract the eigenvalues
from each randomly generated polychoric matrices and the declared
quantile will be extracted. Eigenvalues from polychoric correlation
matrix obtained from real data are also computed and compared, in a
(scree) plot, with the eigenvalues extracted from the simulation (Polychoric
matrices). Recently Cho, Li & Bandalos (2009), showed that in using PA
method, it is important to match the type of the correlation matrix used
to recover the eigenvalues from real data with the type of correlation
matrix used to estimate random eigenvalues. Crossing the type of
correlations (using Polychoric correlation matrix to estimate real
eigenvalues and random simulated Pearson correlation matrice) may result
in a wrong decision (i.e., retaining less non-random factors than the
needed). A comparison with eigenvalues extracted from both randomly
simulated Pearson correlation matrices and real data is also also
included. Finally, for both type of correlation matrix (Polychoric vs
Pearson), the two versions (the classic squared coefficient and the 4th
power coefficient) of Velicer's MAP criterion are calculated.  
} 
\author{
 Fabio Presaghi \email{fabio.presaghi@uniroma1.it} and Marta Desimoni
 \email{marta.desimoni@uniroma1.it} 
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
\keyword{PARALLEL ANALYSIS}
\keyword{POLYCHORIC CORRELATION MATRIX}
\keyword{EXPLORATORY FACTOR ANALYSIS}
