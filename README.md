# Details

Diagnostic measures try to find whether a set of observations is influencing disproportionately the linear regression estimator. Diagnostic measures for independent observations (GLM) are vastly studied and statistical packages exist in Python and R to compute them. In the case of correlated responses (see [Correlated\_responses\\_repo](https://github.com/lokopobit/Correlated_responses) simulation studies asserting that GEE return a less-biased estimator than the classical GLM) no open source code exists to perform these measures.  

In this repository the next diagnostic measures for Generalized Estimating Equations estimator are provided:  

- Leverages for each ID
- Hat matrix
- Exact and approximate cooks distance
- DFFITS
- DFBETAS (not implemented)

A complete explanation of the diagnostic measures covered along with an application to a provided longitudinal dataset can be found in the diagnostic\_measures\_study.pdf file. 


The requried packages are the next ones (3):  

- dplyr
- gepack
- readr

**Note**: the language of the study is spanish.
