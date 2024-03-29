#' Drug consumption data set
#'
#' A dataset collected by Fehrman et al. (2017), freely available on the UCI Machine Learning Repository (Lichman, 2013) containing records of 1885 respondents regarding their use of 18 types of drugs, and their measurements on 12 predictors.
#' #' All predictors were originally categorical and were quantified by Fehrman et al. (2017). The meaning of the values can be found on
#' \url{https://archive.ics.uci.edu/dataset/373/drug+consumption+quantified}.
#' The original response categories for each drug were: never used the drug, used it over a decade ago, or in the last decade, year, month, week, or day.
#' We transformed these into binary response categories, where 0 (non-user) consists of the categories never used the drug and used it over a decade ago and 1 (user) consists of all other categories.
#' @format A data frame with 1185 rows and 32 variables:
#' \describe{
#'   \item{ID}{Respondent ID}
#'   \item{Age}{Age of respondent}
#'   \item{Gender}{Gender of respondent, where 0.48 denotes female and -0.48 denotes male}
#'   \item{Edu}{Level of education of participant}
#'   \item{Country}{Country of current residence of participant}
#'   \item{Ethn}{Ethnicity of participant}
#'   \item{Neuro}{NEO-FFI-R Neuroticism score}
#'   \item{Extr}{NEO-FFI-R Extraversion score}
#'   \item{Open}{NEO-FFI-R Openness to experience score}
#'   \item{Agree}{NEO-FFI-R Agreeableness score}
#'   \item{Consc}{NEO-FFI-R Conscientiousness score}
#'   \item{Impul}{Impulsiveness score measured by BIS-11}
#'   \item{Sensat}{Sensation seeking score measured by ImpSS}
#'   \item{Alc}{Alcohol user (1) or non-user (0)}
#'   \item{Amphet}{Amphetamine user (1) or non-user (0)}
#'   \item{Amyl}{Amyl nitrite user (1) or non-user (0)}
#'   \item{Benzos}{Benzodiazepine user (1) or non-user (0)}
#'   \item{Caff}{Caffeine user (1) or non-user (0)}
#'   \item{Can}{Cannabis user (1) or non-user (0)}
#'   \item{Choco}{Chocolate user (1) or non-user (0)}
#'   \item{Coke}{Coke user (1) or non-user (0)}
#'   \item{Crack}{Crack user (1) or non-user (0)}
#'   \item{Ecst}{Ecstacy user (1) or non-user (0)}
#'   \item{Her}{Heroin user (1) or non-user (0)}
#'   \item{Ket}{Ketamine user (1) or non-user (0)}
#'   \item{Leghighs}{Legal Highs user (1) or non-user (0)}
#'   \item{LSD}{LSD user (1) or non-user (0)}
#'   \item{Meth}{Methadone user (1) or non-user (0)}
#'   \item{Mush}{Magical Mushroom user (1) or non-user (0)}
#'   \item{Nico}{Nicotine user (1) or non-user (0)}
#'   \item{Semeron}{Semeron user (1) or non-user (0), fictitious drug to identify over-claimers}
#'   \item{VSA}{volatile substance abuse user(1) or non-user (0)}
#'
#' }
#' @source \url{https://archive.ics.uci.edu/dataset/373/drug+consumption+quantified}
#' @references \cite{Fehrman, E., Muhammad, A. K., Mirkes, E. M., Egan, V., & Gorban, A. N. (2017). The Five Factor Model of personality and evaluation of drug consumption risk. In Data Science (pp. 231-242). Springer, Cham.}
#' \cite{Lichman, M. (2013). UCI machine learning repository.}
"drugs"


