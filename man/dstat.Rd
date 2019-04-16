\name{dstat}
\alias{dstat}
\concept{Sensitivity analysis}
\concept{Causal inference}
\concept{Adaptive inference}
\concept{Observational studies}
\concept{Rbounds}
\concept{d-statistic}
\concept{Effect modification}

\title{
Sensitivity Analysis Focusing on Subgroups with Demonstrated Insensitivity to Unmeasured Bias
}
\description{
Sensitivity analysis using a d-statistic employing conditional inference to focus on those subgroups with demonstrated insensitivity to unmeasured biases.
}
\usage{
dstat(y, qs = c(1/3, 2/3), gamma = 1, f = NULL, fscore = NULL, fr = 1, alpha = 0.05)
}
\arguments{
  \item{y}{
A numeric vector of treated-minus-control matched pair differences in outcomes.
}
  \item{qs}{
Quantiles of |y| that partly define the d-statistic.    Each coordinate of qs must be a number
strictly between 0 and 1; otherwise, an error will result.  See Details.
}
  \item{gamma}{
The sensitivity parameter, a number gamma>=1.
}
  \item{f}{
If f is not NULL, then it must be a factor that further subdivides y beyond the subdivisions
implied by qs.  If f is not NULL, then the length of f must equal the length of y; otherwise,
an error will result.
}
  \item{fscore}{
If fscore is not NULL, then fscore contains integer scores to be attached to the levels of f.  If f is not NULL but fscore is NULL, then the levels of f are viewed as nominal with equal emphasis.
An error will result if fscore is not NULL but: (i) the scores are not integers, (ii) f is NULL, or (iii) the number of scores does not equal the number levels of f.
}
  \item{fr}{
A nonnegative number.  If fr=0, then the test is simply a group-rank test, using every category, without conditional inference.  The recommended default of fr=1 uses a category only if the proportion of positive y's in this category is at least equal to gamma/(1+gamma), and the conditional inference corrects for selection of categories based on y.  In general, a category is used if the proportion of positive y is at least fr*gamma/(1+gamma), reducing to all categories if fr=0.
}
  \item{alpha}{
Of limited importance, a text message interprets numerical results in terms of rejection or not
of the null hypothesis of no treatment effect in a one-sided, level-alpha test in the presence of a bias in treatment assignment of at most gamma>=1.
}
}
\details{
The method is from Rosenbaum (2019).  The example reproduces aspects of this
manuscript.

The default values of qs, 1/3 and 2/3, are from Brown (1981)'s test.
See Markowski and Hettmansperger (1982) for discussion of other choices.
See Rosenbaum (2015) for comparisons of performance of
different fixed choices of qs; here, a fixed choice is obtained by setting fr=0.

If a pair difference in y is zero, it falls in the lowest quantile of pairs and
therefore receives weight zero along with other pair differences with small |y|.
}
\value{
\item{T }{The test statistic}
\item{comp2 }{The sharp upper bound on the one-sided, exact P-value testing the null hypothesis of no treatment effect in the presence of a bias in treatment assignment of at most gamma.}
\item{scores}{A vector reminding you of the scores, fscore, that you may have attached to the levels of f.}
\item{table}{A table showing how individual categories contribute to the overall test.  The notation in this table is from Rosenbaum (2019).}
\item{summary}{A text summary of the conclusion.}
}
\references{

Brown, B. M. (1981). Symmetric quantile averages and related estimators. Biometrika, 68(1), 235-242.

Lalive, R., Van Ours, J., & Zweimüller, J. (2006). How changes
in financial incentives affect the duration of unemployment.
The Review of Economic Studies, 73, 1009-1038.

Markowski, E. P., Hettmansperger, T. P. (1982). Inference based on simple rank step score statistics for the location model. Journal of the American Statistical Association, 77(380), 901-907.

Noether, G. E. (1973). Some simple distribution-free confidence intervals for the center of a symmetric distribution. Journal of the American Statistical Association, 68(343), 716-719.

Rosenbaum, P. R. (1999). Using quantile averages in matched observational studies. Journal of the Royal Statistical Society: Series C (Applied Statistics), 48(1), 63-78. <doi.org/10.1111/1467-9876.00140>

Rosenbaum, P. R. and Silber, J. H. (2009) Amplification of sensitivity analysis in observational studies.  Journal of the American Statistical Association, 104, 1398-1405.
<doi:10.1198/jasa.2009.tm08470>

Rosenbaum, P. R. (2010) The power of a sensitivity analysis and its limit.  Chapter 14 of
Design of Observational Studies.  NY: Springer. <doi:10.1007/978-1-4419-1213-8_14>

Rosenbaum, P. R. (2015). Bahadur efficiency of sensitivity analyses in observational studies. Journal of the American Statistical Association, 110(509), 205-217. <doi.org/10.1080/01621459.2014.960968>

Rosenbaum, P. R. (2017). Observation and Experiment: An Introduction to Causal Inference.
Cambridge, MA: Harvard University Press.

Rosenbaum, P. R. (2019). A highly adaptive test for matched
observational studies.  Manuscript.
}
\author{
Paul R. Rosenbaum
}


\examples{
# First example is from Rosenbaum (2019)
data("lalive")
attach(lalive)
y<-log2((1+dur[after==1])/52)-log2((1+dur[after==0])/52)
dstat(y,qs=c(1/3,2/3),fr=0,gamma=1.15) # Brown's (1981) test
dstat(y,qs=c(2/3),fr=0,gamma=1.25) # Noether's (1973) Test
#Amplification: see Rosenbaum and Silber (2009), Rosenbaum (2017, Table 9.1)
amplify(1.25,2)

bothseasonal<-(2==(seasonal[after==1]+seasonal[after==0]))*1
bothseasonal<-factor(bothseasonal,levels=1:0,
   labels=c("S","O"),ordered=TRUE)
straddle<-1*((pmin(dur[after==1],dur[after==0])<=(bdur[after==0]))
   &(pmax(dur[after==1],dur[after==0])>(bdur[after==0])))
straddle<-factor(straddle,levels=c(0,1),
   labels=c("N","Y"),ordered=TRUE)
dose<-1*((bdur[after==1]-bdur[after==0])>9.5)
dose<-factor(dose,levels=c(0,1),
   labels=c("L","H"),ordered=TRUE)
f<-bothseasonal:dose:straddle
dstat(y,qs=c(1/3,2/3),f=f,gamma=1.25)
# Reproduces Table 2 in Rosenbaum (2019)
dstat(y,qs=c(1/3,2/3),f=f,gamma=1.45)
amplify(1.45,c(2,2.5,3,4))

# Doubling the weight for high-dose matched pairs
levels(f)
fs<-c(1,1,2,2,1,1,2,2)
dstat(y,qs=c(1/3,2/3),f=f,fscore=fs,gamma=1.45)

rm(y,f,dose,straddle,bothseasonal)
detach(lalive)
rm(lalive)

# Second example uses dental data
data(dental)
attach(dental)
f<-age:dose
levels(f)
# Doubles the weight at high dose using fscore
# For qs = (.4,.8), see Markowski and Hettmansperger (1982)
dstat(y,qs=c(.4,.8),gamma=4.25,f=f,fscore=c(1,2,1,2))
rm(f)
detach(dental)
rm(dental)
}

\keyword{htest }
