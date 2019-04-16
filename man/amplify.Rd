\name{amplify}
\alias{amplify}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Amplification}
\concept{Amplify}
\title{
Amplification of sensitivity analysis in observational studies.
}
\description{
Uses the method in Rosenbaum and Silber (2009) to interpret a value of the sensitivity parameter gamma.  Each value of gamma amplifies to a curve (lambda,delta) in a two-dimensional sensitivity analysis, the inference being the same for all points on the curve.
That is, a one-dimensional sensitivity analysis in terms of gamma has a two-dimensional interpretation in terms of (lambda,delta).
}
\usage{
amplify(gamma, lambda)
}
\arguments{
  \item{gamma}{
	gamma > 1 is the value of the sensitivity parameter, for instance the parameter in senmv.  length(gamma)>1 will generate an error.
}
  \item{lambda}{
	lambda is a vector of values > gamma.  An error will result unless lambda[i] > gamma > 1 for every i.
}
}
\details{
A single value of gamma, say gamma = 2.2 in the example, corresponds to a curve of values of (lambda, delta), including (3, 7),
(4, 4.33), (5, 3.57), and (7, 3) in the example.  An unobserved covariate that is associated with a lambda = 3 fold increase in the odds of treatment and a delta = 7 fold increase in the odds of a positive pair difference is equivalent to gamma = 2.2.

The curve is gamma = (lambda*delta + 1)/(lambda+delta).  Amplify is given one gamma and a vector of lambdas and solves for the vector of deltas.  The calculation is elementary.

This interpretation of gamma is developed in detail in Rosenbaum and Silber (2009), and it makes use of Wolfe's (1974) family of semiparametric deformations of an arbitrary symmetric distribuiton.

Strictly speaking, the amplification describes matched pairs, not matched sets.  The senm function views a k-to-1 matched set with k controls matched to one treated individual as a collection of k correlated treated-minus-control matched pair differences; see Rosenbaum (2007).  For matched sets, it is natural to think of the amplification as describing any one of the k matched pair differences in a k-to-1 matched set.

The curve has asymptotes that the function amplify does not compute: gamma corresponds with (lambda,delta) = (gamma, Inf) and (Inf, gamma).

A related though distict idea is developed in Gastwirth et al (1998).  The two approaches agree when the outcome is binary, that is, for McNemar's test.
}
\value{
Returns a vector of values of delta of length(lambda) with names lambda.
}
\references{
Gastwirth, J. L., Krieger, A. M., Rosenbaum, P. R. (1998) Dual and
simultaneous sensitivity analysis for matched pairs. Biometrika, 85, 907-920.

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (R package sensitivitymv) <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2016) Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. <doi:10.1214/16-AOAS942>.

Rosenbaum, P. R. and Silber, J. H. (2009) Amplification of sensitivity analysis in observational studies.  Journal of the American Statistical Association, 104, 1398-1405.
<doi:10.1198/jasa.2009.tm08470>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)

Wolfe, D. A. (1974) A charaterization of population weighted symmetry and related results.  Journal of the American Statistical Association, 69, 819-822.
}
\author{
Paul R. Rosenbaum
}
\note{
The amplify function is also in the sensitivitymv package where a different example is used.
}

\examples{
amplify(4,c(5,6,7))
}

