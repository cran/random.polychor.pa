\name{Parallel-Analysis-of-Polychoric-Correlations}
\alias{random.polychor.pa-package}
\docType{package}
\title{
  A Parallel Analysis with Random Polychoric Correlation Matrices
}
\description{
  The function performs a parallel analysis using simulated polychoric correlation matrices. The function will extract the eigenvalues from each random generated polychoric correlation matrix and from the polychoric correlation matrix of real data. A plot comparing eigenvalues   extracted from the specified real data with simulated data will help determine which of real eigenvalue outperform random data. A series of matrices comparing MAP vs PA-Polychoric vs PA-Pearson correlations methods, FA vs PCA solutions are finally presented.
}
\details{
  \tabular{ll}{
    Package: \tab random.polychor.pa\cr
    Type: \tab Package\cr
    Version: \tab 1.1.4-1\cr
    Date: \tab 2014-08-12\cr
    License: \tab GPL Version 2 or later\cr
    LazyLoad: \tab yes\cr
}
The function perform a parallel analysis (Horn, 1965) using randomly simulated polychoric correlations and generates \code{nrep} random samples of the same dimension of the empirical provided data.matrix(i.e, with the same number of participants and of variables). Items are allowed to have varying number of categories. The function will extract the eigenvalues from each randomly generated polychoric matrices and the declared quantile will be extracted. Eigenvalues from polychoric correlation matrix obtained from real data are also computed and compared, in a (scree) plot, with the eigenvalues extracted from the simulation (Polychoric matrices). Recently Cho, Li & Bandalos (2009), showed that in using PA method, it is important to match the type of the correlation matrix used to recover the eigenvalues from real data with the type of correlation matrix used to estimate random eigenvalues. Crossing the type of correlations (using Polychoric correlation matrix to estimate real eigenvalues and random simulated Pearson correlation matrices) may result in a wrong decision (i.e., retaining less non-random factors than the needed). A comparison with eigenvalues extracted from both randomly simulated Pearson correlation matrices and real data is also also included. Finally, for both type of correlation matrix (Polychoric vs Pearson), the two versions (the classic squared coefficient and the 4th power coefficient) of Velicer's MAP criterion are calculated. In version 1.1.1, a minor bug in the regarding the estimated time needed to complete the simulation is fixed. Also in this version, the function is now able to manage supplied data.matrix in which variables representing factors (i.e., variables with ordered categories) are present and may cause an error when the Pearson correlation matrix is calculated. As the \code{poly.mat()} function used to calculate the polychoric correlation matrix is going to be deprecated in favor of poly.mat(polychoric()) function, the random.polychor.pa was consequently updated (version 1.1.2) to account for changes in \code{psych()} package.  Version 1.1.3 tackles two problems signaled by users: 1) the possibility to make available the results of simulation for plotting them in other software. Now the random.polychor.pa will show, upon request, all the data used in the scree-plot. 2) The function \code{polichoric()} of the \code{psych()} package does not handle data matrices that include 0 as possible category and will cause the function to stop with error. So a check for the detection of the 0 code within the provided data.matrix is now added and will cause the random.polychor.pa function to stop with a warning message. In version 1.1.3.5 a parameter was added, \code{diff.fact} in order to simulate random data set with the same probability of observing each category for each variable as that observed in the provided (empirical) data set. This parameter allows to simulate data sets that have the same item difficulties distribution as well as the same difficulty factors of the real data (Ledesma and Valero-Mora, 2007). Finally the search for zeroes within the provided data file was removed, so data with zeroes are now accepted. In version 1.1.3.6 a check for the range of quantile (between 0 and 1) was added. 
In version 1.1.4 was added the possibility to choose between uniform and multinomial distribution, between random and boostrap distribution and between single or multiple sample parallel analysis. Also two switches were added: one allows to print or a default output (including only the number of factors to retain following Parallel Analysis and the scree-plot) or a full output (including the eigenvalues distributions of simulated and empirical data sets); the other switch allows to print a further Fit Statistics for the factor solutions indicated by the Parallel Analysis. 

The default value of the \code{comparison} is set to \code{random} and this means that random uniform samples (\code{runif()}) will be computed. Two types of distribution are available: the uniform distribution and the multinomial distribution (Liu and Rijmen, 2008). To switch between them just set accordingly the parameter \code{distr}. As method of resampling is now also available the bootstrap method that is based on random permutations of cases obtained via the \code{boot()} function. 

Multigroup version of random samples may be obtained by setting the parameter \code{comparison} to \code{random-mg}. To perform multigroup Parallel Analysis simulating random samples with multinomial distributions the parameter  \code{distr} should be set to \code{multinomial}. Multigroup Parallel Analysis is also available for the bootstrap method by setting the \code{comparison} parameter to \code{bootstrap-mg}. When Multigroup Parallel Analysis is used, the first variable (column) of the dataset should be reserved to the factor used to identify the different samples (i.e., gender, agre groups, nationalities, etc.).  

Fit statistics (Chi-squared, TLI, RMSEA, RMR, BIC) for all factor solutions indicated by Parallel Analysis are available when the parameter \code{fit.pa} is set to \code{TRUE}. Consider that in estimating the these Fit Statistic indexes not allways converge and this may break the computations for the parallel analysis. In this case switch the parameter off. Moreover the values of these indexes are not allways within the expected range. 

The \code{print.all} parameter when set to  \code{TRUE} will show the table of eigenvalue distributions for simulated (random or boostrap) and empirical data. For boostrap simulation, also the bias and standard error will be printed. If the \code{diff.fact} is set to \code{TRUE} when random simulations are requested, then also the weighting matrix used to reproduce the item frequencies in random data set will be printed.
}
\author{
  Fabio Presaghi \email{fabio.presaghi@uniroma1.it} and Marta Desimoni \email{marta.desimoni@uniroma1.it}
}
\references{
  Cho, S.J., Li, F., & Bandalos, D., (2009). Accuracy of the Parallel Analysis Procedure With Polychoric Correlations. Educational and Psychological Measurement, 69, 748-759.

  Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 32, 179-185.
  
  Ledesma RD, Valero-Mora P (2007) Determining the number of factors to retain in EFA: an  easy to use computer program for carrying out parallel analysis. Practical Assessment, Research & Evaluation 12:1–11

  O\kbd{'}Connor, B. P. (2000). SPSS and SAS programs for determining the number of components using parallel analysis and Velicer\kbd{'}s MAP test. Behavior Research Methods, Instrumentation, and Computers, 32, 396-402.
  
  Reckase, M.D. (2009). Multidimensional Item Response Theory. Springer. 

  Velicer, W. F. (1976). Determining the number of factors from the matrix of partial correlations. Psychometrika, 41, 321-327.

  Velicer, W. F., Eaton, C. A., & Fava, J. L. (2000). Construct explication through factor or component analysis: A review and evaluation of alternative procedures for determining the number of factors or components. In R. D. Goffin & E. Helmes (Eds.), Problems and solutions in human assessment: Honoring Douglas N. Jackson at seventy (pp. 41-72). Norwell, MA: Kluwer Academic.
}

\keyword{PARALLEL ANALYSIS}
\keyword{POLYCHORIC CORRELATION MATRIX}
\keyword{UNIFORM DISTRIBUTION}
\keyword{MULTINOMIAL DISTRIBUTION}
\keyword{MULTIGROUP}
\keyword{BOOTSTRAP}
\keyword{FIT INDEXES}
\keyword{EXPLORATORY FACTOR ANALYSIS}