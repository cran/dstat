\name{dental}
\alias{dental}
\docType{data}
\concept{Observational study}
\concept{Causal inference}
\concept{Sensitivity analysis}
\concept{Matched pairs}
\title{
Dental Problems Caused by Smoking
}
\description{
Data from NHANES 2011-2012 containing 441 matched pairs of a daily cigarette smoker and a never smoker, recording the extent of periodontal disease. Pairs were matched for sex, age, black race, education in five categories, and ratio of family income to the poverty level.
}
\usage{data("dental")}
\format{
  A data frame with 441 observations on the following 5 variables.
  \describe{
    \item{\code{smoker}}{Periodontal disease in the daily-smoker}
    \item{\code{control}}{Periodontal disease in the never smoker}
    \item{\code{y}}{Smoker-minus-control pair difference}
    \item{\code{age}}{Age \code{<= 50} \code{> 50}}
    \item{\code{dose}}{Cigarettes smoked per day by the smoker \code{<10} \code{>=10}}
  }
}
\details{
Excluding wisdom teeth, 6 measurements are taken for each tooth that is present, up to 28 teeth. Following Tomar and Asma (2000), a measurement indicates periodontal disease if either there is a loss of attachment of at least 4mm or a pocket depth of at least 4mm. The first individual has 11 measurements indicative of periodontal disease, out of 106 measurements, so pcteither is 100*11/106 = 10.38 percent.


}
\source{
Data are from the National Health and Nutrition Examination Survey 2011-2012 and were used as an example in Rosenbaum (2017). In the second edition of Design of Observational Studies, these data are discussed in the chapter entitled Evidence Factors. Although the same 2x441 individuals were used here and in Rosenbaum (2017), the pairing was slightly changed to be exact for age>50.
}
\references{
Rosenbaum, P. R. (2017) The general structure of evidence factors in observational studies.
Statist Sci 32, 514-530. <doi:10.1214/17-STS621>

Tomar, S. L. and Asma, S. (2000) Smoking attributable periodontitis in the US: Findings from NHANES III.
J Periodont 71, 743-751.

US National Health and Nutrition Examination Survey 2011-2012. www.cdc.gov/nchs/nhanes/index.htm
}
\examples{
data(dental)
attach(dental)
boxplot(y~dose:age)
abline(h=0)
detach(dental)
rm(dental)
}
\keyword{datasets}
