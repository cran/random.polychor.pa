\name{NEWS}
\title{News for Package \pkg{random.polychor.pa}}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\description{
  The function performs a parallel analysis using simulated polychoric correlation matrices. The function will extract the eigenvalues from each random generated polychoric correlation matrix and from the polychoric correlation matrix of real data. A plot comparing eigenvalues   extracted from the specified real data with simulated data will help determine which of real eigenvalue outperform random data. A series of matrices comparing MAP vs PA-Polychoric vs PA-Pearson correlations methods, FA vs PCA solutions are finally presented. Random data sets are simulated assuming or a uniform or a multinomial distribution or via the bootstrap method of resampling (i.e., random permutations of cases). Also Multigroup Parallel analysis is made available for random (uniform and multinomial distribution and with or without difficulty factor) and bootstrap methods. An option to choose between default or full output is also available as well as a parameter to print Fit Statistics (Chi-squared, TLI, RMSEA, RMR and BIC) for the factor solutions indicated by the Parallel Analysis.
}

\section{Changes in version 1.1.4}{
\itemize{
Version 1.1.4 added a number of cahnges: 
\item a parameter \code{distr} allows to shift the simulation from uniform distribtuion to multinomial distribution
\item Bootstrap method (with random permutations of cases) was added
\item Multi-Group for random and bootstrap version of the Parallel Analysis is made available
\item Fit statistics (Chi-squared, TLI, RMSEA, RMR, BIC) for all factor solution indicated by Parallel Analysis
\item option to print a default output (number of factors indicated by Parallel Analysis) or a full output (adds: matrices of simulated and empirical eigenvalues for random, bootstrap, and multigroup)
}
}

\section{Changes in version 1.1.3.6}{
\itemize{
\item In version 1.1.3.6 a check for the range of quantile (between 0 and 1) was added.
}
}

\section{Changes in version 1.1.3.5}{
\itemize{
\item The search for zeroes within the provided datafile was removed, so data with zeroes are now accepted.
}
}

\section{Changes in version 1.1.3.5}{
\itemize{
\item In version 1.1.3.5 a paramether was added, \code{diff.fact}, in order to simulate random dataset with the same probability of observing each category for each variable as that         observed in the provided (empirical) dataset.
}
}

\section{Changes in version 1.1.3}{
\itemize{
Version 1.1.3 tackles two problems signalled by users: 
\item the possibility to make available the results of simulation for plotting them in other software. Now the \code{random.polychor.pa()} will show, upon request, all the data used in the scree-plot.
\item The function \code{polichoric()} of the \code{psych()} package does not handle data matrices that include 0 as possible category and will cause the function to stop with error. So a check for the detection of the 0 code within the provided data.matrix is now added and will cause the random.polychor.pa function to stop with a warning message.
}
}

\section{Changes in version 1.1.2}{
\itemize{
\item Version 1.1.2 simply has updated the function that calculates the polychoric correlation matrix due to changes in the \code{psych()} package.
}
}

\section{Changes in version 1.1.1}{
\itemize{
\item Version 1.1.1, fixed a minor bug in the regarding the estimated time needed to complete the simulation. 
\item Also in this version, the function is now able to manage supplied data.matrix in which variables representing factors (i.e., variables with ordered categories) are present and may cause an error when the Pearson correlation matrix is calculated. 

}
}



