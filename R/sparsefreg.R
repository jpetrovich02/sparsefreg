#'sparsefreg: A package for performing scalar-on-function regression with sparsely-observed functional covariates
#'
#'The sparsefreg package implements the MISFIT (Multiple Imputation for Sparsely-sampled
#'Functions at Irregular Times) approach, used for scalar-on-function regression when
#'the functional covariate is observed sparsely. MISFIT takes a missing data approach to
#'imputing (over a denser grid) a functional covariate which is only observed sparsely on
#'a grid. More specifically, the curves are represented by their Karhunen-Loeve expansion
#'where the scores are imputed and can be used to reconstruct the curves, and/or directly
#'used to fit an FPC regression model.
#'
#'
#'@docType package
#'@name sparsefreg

NULL
