\name{lalive}
\alias{lalive}
\concept{Matched pairs}
\concept{Sensitivity analysis}
\concept{Observational study}
\docType{data}
\title{
Unemployment Duration Following an Increase in Unemployment Benefits
}
\description{
Data from a study by Lalive, van Ours and Zweimüller (2006) concerning
the duration of unemployment before and after an increase in
unemployment benefits, both the benefit amount and the duration
of benefits.  The original study takes account of many relevant
considerations not included in the current subset of the data
The data were used as methodological example in Rosenbaum (2019).
}
\usage{data("lalive")}
\format{
  A data frame with 2782 observations on the following 17 variables.
  \describe{
    \item{\code{id}}{ID number}
    \item{\code{mset}}{Matched pair, 1,2,...,1391.}
    \item{\code{after}}{Treatment indicator, 1=after
          benefits increase, 0=before benefits increase}
    \item{\code{type}}{a factor with levels \code{PBD and RR}}
    \item{\code{dur}}{Duration of unemployment in weeks.}
    \item{\code{bdur}}{Duration of unemployment benefits in weeks}
    \item{\code{e3_5}}{1 if worked for at least 3 of the past 5 years,
             0 otherwise.}
    \item{\code{lehre}}{1 if apprenticeship, 0 otherwise}
    \item{\code{married}}{1 if married, 0 otherwise}
    \item{\code{divorced}}{1 if divorced, 0 otherwise}
    \item{\code{bc}}{1 if lost a blue colar job, 0 otherwise}
    \item{\code{seasonal}}{1 if lost a seasonal job, 0 otherwise}
    \item{\code{manuf}}{1 if lost a manufacturing job, 0 otherwise}
    \item{\code{age}}{Age in years}
    \item{\code{nwage_pj}}{Wage in the prior job in Austrian schillings}
    \item{\code{educ}}{0 if primary education, 1 if secondary education, 2 if tertiary education}
    \item{\code{propensity}}{An estimated propensity score}
  }
}
\details{
The data are from Lalive, van Ours and Zweimüller (2006), by way of the web-page for
the textbook Cahuc, P., Carcillo, S. and Zylberberg, A. (2014).

In August 1989, Austria increased its unemployment benefits for certain categories
of workers.  The category considered here, type=PBD and RR, had an increase
in the duration of unemployment benefits and an increase in unemployment
compensation.  There are two groups, those unemployed in the two years
before the benefit increase, after=0, and those unemployed in the two
years after the increase, after=1.

The data are 1391 matched pairs, matached for e3_5, lehre, married, divorced, bc,
seasonal, manuf, age, nwage_pj, and educ, with fine balance for quintiles of
the propensity score.  All are men, and none were temporarily laid off.  The
matching used a simplified version of the method in Rosenbaum (2017).

The original study by Lalive et al. (2006) sensibly
takes account of many relevant
considerations not included in the current subset of the data.
The limited data available here were used to illustrate certain
methodological issues in Rosenbaum (2019).
}
\source{
Lalive, R., Van Ours, J., & Zweimüller, J. (2006).
}
\references{
Cahuc, P., Carcillo, S. and Zylberberg, A. (2014). Labor Economics,
Second Edition.  Cambridge, MA: MIT Press.
https://mitpress.mit.edu/books/labor-economics-second-edition

Lalive, R., Van Ours, J., & Zweimüller, J. (2006). How changes
in financial incentives affect the duration of unemployment.
The Review of Economic Studies, 73, 1009-1038.
<doi:10.1111/j.1467-937X.2006.00406.x>

Rosenbaum, P. R. (2017). Imposing minimax and quantile
constraints on optimal matching in observational studies.
Journal of Computational and Graphical Statistics, 26, 66-78.

Rosenbaum, P. R. (2019). A highly adaptive test for matched
observational studies.  Manuscript.
}
\examples{
data(lalive)
attach(lalive)
# covariate balance
boxplot(propensity~after,names=c("Before","After"),ylab="Propensity Score")
boxplot(age~after,names=c("Before","After"),ylab="Age")
boxplot(nwage_pj~after,names=c("Before","After"),ylab="Prior Wage")
table(after,seasonal)
# outcome
y<-log2((1+dur[after==1])/52)-log2((1+dur[after==0])/52)
boxplot(y,ylab="Pair Difference in base 2 logs",
    main="Unemployment Duration")
abline(h=c(-1,0,1),lty=2)
rm(y)
detach(lalive)
}
\keyword{datasets}
