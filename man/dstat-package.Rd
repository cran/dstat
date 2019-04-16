\name{dstat-package}
\alias{dstat-package}
\docType{package}
\title{
\packageTitle{dstat}
}
\description{
\packageDescription{dstat}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{dstat}
\packageIndices{dstat}
The package provides a sensitivity analysis for a conditional test of the null hypothesis of no treatment effect in a matched observational study in which the unmeasured bias in treatment assignment is quantified by a sensitivity parameter gamma>=1.  The test uses only those categories of pairs that demonstrate insensitivity to a bias of magnitude gamma, correcting for data-dependent selection of categories by conditional inference.  The main function in the package is dstat().
}
\author{
\packageAuthor{dstat}

Maintainer: \packageMaintainer{dstat}
}
\references{
Rosenbaum, P. R. (1999). Using quantile averages in matched observational studies. Journal of the Royal Statistical Society: Series C (Applied Statistics), 48(1), 63-78. <doi.org/10.1111/1467-9876.00140>
}
\keyword{ package }

\examples{
data("dental")
attach(dental)
head(dental)
dstat(y,gamma=4.1,f=dose:age,fscore=c(1,1,2,2))
amplify(4,c(5,6,7))
detach(dental)
}
