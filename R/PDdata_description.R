#' Description of PD data
#' @docType data
#' 
#' @usage data(PDdata)
#' 
#' ID: subject ID
#'
#' Tensor responses (teeth x surfaces x biomarkers) ##
#' Tooth -- 
#' 1-2, 13-16, 27-28 : Molars (total 48 surfaces)
#' 3-4, 11-12, 17-18, 25-26 : Pre-molars (total 48 surfaces)
#' 5, 10, 19, 24 : Canines (total 24 surfaces)
#' 6-9, 20-23 : Incisors (48 surfaces)

#' Side & Surface --
#' 1 : Buccal - Distal (Disto-Buccal)
#' 2 : Buccal - Buccal (mid-buccal)
#' 3 : Buccal(facial) - Mesial (Mesio-Buccal)
#' 4 : Lingual - Distal (Disto-Lingual)
#' 5 : Lingual - Lingual (mid-lingual)
#' 6 : Lingual(tongue) - Mesial (Mesio-Lingual)

#' Two Biomarkers 
#' PPD: Periodontal Pocket Depth
#' CAL: Clinical Attachment Level

#' Covariates
#' Age: in years (continuous variable)
#' Gender: 1 = Female/0 = Male (binary variable)
#' BMI: Body Mass Index, measure of obesity with obese when BMI >= 30 (binary variable)
#' Smoker: 1 = Current or previous smoker, 0 = Never smoker (binary variable)
#' HbA1c: measure of Type-2 diabetes, 1 = Uncontrolled, 0 = Controlled (binary variable)
#' 
#' 
#' @example 
#' data(PDdata)
