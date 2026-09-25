###############################################################################
# Figure 4. Discrimination of prevalent cardiovascular disease (CVD), fatty
# liver disease (FLD), and chronic kidney disease (CKD) by 27
# cardiometabolic indices, overall and by BMI category
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Analysis: unadjusted area under the receiver-operating-characteristic
# curve (AUC) of each index for each outcome, with 95% CIs estimated by the
# DeLong method, overall and by BMI category. The direction of comparison
# is determined from the median values in participants with and without the
# outcome (pROC default, direction = "auto").
#
# Indices (n = 27), grouped as:
#   - Biochemical (n = 8)
#   - Anthropometric and body-composition (n = 11)
#   - Composite, integrating both (n = 8)
# Formulas and references are provided in Table 2 of the manuscript.
#
# Output (printed to the graphics device and console; no files are written):
#   - Figure 4: AUC (95% CI) of each index by BMI category (points) and
#     outcome (panels CVD, FLD, CKD); indices ordered within each group by
#     overall AUC for CVD
#   - Console summaries: AUC (95% CI) of each index by outcome and BMI
#     category, and the six indices with the highest overall AUC per outcome
#
# Analytical population: BMI >= 18.5 kg/m2, sex recorded as female or male,
# and metabolic health status resolvable under both Definition A and
# Definition B (same as Table 1). Each AUC uses participants with available
# data for the index and the outcome.
#
# R version 4.6.0
# Packages: dplyr, readr, stringr, tidyr, pROC, ggplot2
###############################################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(pROC)
library(ggplot2)

# Path to the oriGen questionnaire data file (tab-separated).
# Column names correspond to the original oriGen data release.
data_path <- "path/to/oriGen_questionnaire_data.tsv"
raw_data  <- read_tsv(data_path, show_col_types = FALSE)

# -- Helper functions ---------------------------------------------------------

# Harmonises free-text responses: upper case, no accents, single spaces.
normalise_text <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- trimws(x)
  x <- str_to_upper(x)
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  gsub("[[:space:]]+", " ", x)
}

# Converts yes/no responses (Spanish or English) to 1/0.
# Responses that cannot be interpreted are set to NA.
safe_bool <- function(x) {
  y <- normalise_text(x)
  case_when(
    y %in% c("TRUE", "VERDADERO", "SI", "1", "YES") ~ 1L,
    y %in% c("FALSE", "FALSO", "NO", "0")           ~ 0L,
    TRUE                                            ~ NA_integer_
  )
}

# Values outside the prespecified physiological range [lo, hi] are set to NA.
# Applied to the raw variables before any index is calculated.
clip_range <- function(x, lo, hi) if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)

# -- Data preparation ---------------------------------------------------------

dat <- raw_data %>%
  transmute(
    SEX_RAW = str_to_upper(trimws(as.character(PARTICIPANTE.SEXO))),
    SEX = case_when(
      SEX_RAW %in% c("H", "HOMBRE", "MASCULINO", "MALE", "MAN")         ~ "Male",
      SEX_RAW %in% c("M", "MUJER", "F", "FEMENINO", "FEMALE", "WOMAN") ~ "Female",
      TRUE ~ NA_character_
    ),
    BMI        = suppressWarnings(as.numeric(ANTROPOMETRIA.IMC)),
    WEIGHT     = suppressWarnings(as.numeric(INBODY.Weight)),
    WAIST      = suppressWarnings(as.numeric(ANTROPOMETRIA.CIRCUNFERENCIA_CINTURA)),
    HEIGHT_RAW = suppressWarnings(as.numeric(ANTROPOMETRIA.ESTATURA)),
    # Height recorded in cm or m: values > 3 are taken as cm
    HEIGHT_M   = if_else(HEIGHT_RAW > 3, HEIGHT_RAW / 100, HEIGHT_RAW),
    HEIGHT_CM  = if_else(HEIGHT_RAW > 3, HEIGHT_RAW, HEIGHT_RAW * 100),
    TG         = suppressWarnings(as.numeric(METABOLICOS.TRIGLICERIDOS)),
    HDL        = suppressWarnings(as.numeric(METABOLICOS.HDL)),
    LDL        = suppressWarnings(as.numeric(METABOLICOS.LDL)),
    TC         = suppressWarnings(as.numeric(METABOLICOS.COLESTEROL)),
    SBP        = suppressWarnings(as.numeric(METABOLICOS.SISTOLICA)),
    DBP        = suppressWarnings(as.numeric(METABOLICOS.DIASTOLICA)),
    WHR        = suppressWarnings(as.numeric(INBODY.WHR_Waist_Hip_Ratio)),
    FASTING_HOURS = suppressWarnings(as.numeric(ANTROPOMETRIA.HORAS_AYUNO)),
    # Glucose is considered a fasting value only after >= 8 h of fasting.
    GLUCOSE    = if_else(
      FASTING_HOURS >= 8,
      suppressWarnings(as.numeric(METABOLICOS.GLUCOSA)),
      NA_real_
    ),
    BFM        = suppressWarnings(as.numeric(INBODY.BFM_Body_Fat_Mass)),
    SMM        = suppressWarnings(as.numeric(INBODY.SMM_Skeletal_Muscle_Mass)),

    # Self-reported diabetes and glucose-lowering treatment
    diabetes                 = safe_bool(HISTORIA_MEDICA.DIABETES),
    diabetes_type1           = safe_bool(HISTORIA_MEDICA.DIABETES_TIPO1),
    diabetes_type2           = safe_bool(HISTORIA_MEDICA.DIABETES_TIPO2),
    glucose_lowering_oral    = safe_bool(HISTORIA_MEDICA.PASTILL_CONTROL_AZUCAR),
    glucose_lowering_insulin = safe_bool(HISTORIA_MEDICA.INSULINA_CONTROL_AZUCAR),
    glucose_lowering_any     = safe_bool(HISTORIA_MEDICA.TRATAMIENTO_AZUCAR),

    # Antihypertensive and lipid-lowering treatment
    antihypertensive = safe_bool(HISTORIA_MEDICA.PASTILLAS_PRESION_ALTA),
    lipid_lowering   = safe_bool(HISTORIA_MEDICA.TRATA_MEDICAMENTO),

    # Outcome items (self-reported physician diagnosis)
    cvd_mi         = safe_bool(HISTORIA_MEDICA.INFARTO_AGUDO_MIOCARDIO),
    cvd_stroke     = safe_bool(HISTORIA_MEDICA.INFARTO_CEREBRAL),
    cvd_ischaemic  = safe_bool(HISTORIA_MEDICA.ENF_CORONARIA_ISQ),
    cvd_pad        = safe_bool(HISTORIA_MEDICA.ENF_ARTERIAL_PERIF),
    cvd_athero     = safe_bool(HISTORIA_MEDICA.ATEROESCLEROSIS),
    cvd_hf         = safe_bool(HISTORIA_MEDICA.INSUFICIENCIA_CARDIACA),
    fatty_liver    = safe_bool(HISTORIA_MEDICA.HIGADO_GRASO),
    ckd_diagnosis  = safe_bool(HISTORIA_MEDICA.ENFERMEDAD_RENAL_CRONICA),
    ckd_renal_insf = safe_bool(HISTORIA_MEDICA.INSUFICIENCIA_RENAL),
    ckd_dialysis   = safe_bool(HISTORIA_MEDICA.DIALISIS_POR_RENAL),
    ckd_haemodial  = safe_bool(HISTORIA_MEDICA.HEMODIALISIS)
  ) %>%

  # Prespecified physiological ranges (same as Table 1), applied before
  # calculating the indices
  mutate(
    WEIGHT    = clip_range(WEIGHT,     30,  300),
    HEIGHT_CM = clip_range(HEIGHT_CM, 130,  220),
    HEIGHT_M  = clip_range(HEIGHT_M, 1.30, 2.20),
    BMI       = clip_range(BMI,        12,   70),
    WAIST     = clip_range(WAIST,      50,  200),
    GLUCOSE   = clip_range(GLUCOSE,    40,  600),
    TG        = clip_range(TG,         20, 2000),
    TC        = clip_range(TC,         50,  500),
    HDL       = clip_range(HDL,        10,  150),
    LDL       = clip_range(LDL,        10,  400),
    SBP       = clip_range(SBP,        70,  260),
    DBP       = clip_range(DBP,        40,  160),
    WHR       = clip_range(WHR,       0.6,  1.5),
    SMM       = clip_range(SMM,         5,   70),
    BFM       = clip_range(BFM,         2,  120)
  ) %>%
  mutate(
    # Composite outcomes: 1 if any item is positive; 0 if none is positive
    # and at least one item is evaluable; NA only if all items are missing.
    CVD = case_when(
      cvd_mi == 1L | cvd_stroke == 1L | cvd_ischaemic == 1L |
        cvd_pad == 1L | cvd_athero == 1L | cvd_hf == 1L ~ 1L,
      !is.na(cvd_mi) | !is.na(cvd_stroke) | !is.na(cvd_ischaemic) |
        !is.na(cvd_pad) | !is.na(cvd_athero) | !is.na(cvd_hf) ~ 0L,
      TRUE ~ NA_integer_
    ),

    FLD = case_when(
      fatty_liver == 1L ~ 1L,
      fatty_liver == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),

    CKD = case_when(
      ckd_diagnosis == 1L | ckd_renal_insf == 1L |
        ckd_dialysis == 1L | ckd_haemodial == 1L ~ 1L,
      !is.na(ckd_diagnosis) | !is.na(ckd_renal_insf) |
        !is.na(ckd_dialysis) | !is.na(ckd_haemodial) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # BMI category (NA for BMI < 18.5 kg/m2 or missing BMI)
    BMI_CATEGORY = case_when(
      BMI >= 18.5 & BMI < 25 ~ "Normal weight",
      BMI >= 25   & BMI < 30 ~ "Overweight",
      BMI >= 30              ~ "Obesity",
      TRUE                   ~ NA_character_
    ),

    # -------------------------------------------------------------------------
    # METABOLIC HEALTH DEFINITIONS (same as Table 1)
    # Used only to restrict the analysis to the same analytical population
    # as Table 1 and Figures 1-3.
    #
    # Component coding: 1 = present; 0 = absent.
    #   - 1 if any of its evaluable criteria is met (measured value meets the
    #     cut-off, or the corresponding diagnosis or treatment is reported);
    #   - 0 if none is met and at least one criterion is evaluable;
    #   - NA only if all criteria are missing, uninterpretable, or out of range.
    # A definition is resolvable only when all of its components are non-missing.
    # -------------------------------------------------------------------------

    # Definition B (MUH if 1 or more of 3 components)
    defB_diabetes = case_when(
      diabetes == 1L ~ 1L,
      diabetes == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),

    defB_bp = case_when(
      (!is.na(SBP) & SBP >= 130) |
        antihypertensive == 1L ~ 1L,
      !is.na(SBP) |
        !is.na(antihypertensive) ~ 0L,
      TRUE ~ NA_integer_
    ),

    defB_whr = case_when(
      SEX == "Female" & !is.na(WHR) & WHR >= 0.95 ~ 1L,
      SEX == "Female" & !is.na(WHR) & WHR <  0.95 ~ 0L,
      SEX == "Male"   & !is.na(WHR) & WHR >= 1.03 ~ 1L,
      SEX == "Male"   & !is.na(WHR) & WHR <  1.03 ~ 0L,
      TRUE ~ NA_integer_
    ),

    defB_valid_components = rowSums(!is.na(cbind(
      defB_diabetes, defB_bp, defB_whr
    ))),

    defB_n = if_else(
      defB_valid_components == 3L,
      defB_diabetes + defB_bp + defB_whr,
      NA_integer_
    ),

    defB_muh = case_when(
      !is.na(defB_n) & defB_n == 0L ~ 0L,
      !is.na(defB_n) & defB_n >= 1L ~ 1L,
      TRUE ~ NA_integer_
    ),

    # Definition A (MUH if 2 or more of 5 components)
    defA_waist = case_when(
      SEX == "Male"   & !is.na(WAIST) & WAIST >= 90 ~ 1L,
      SEX == "Male"   & !is.na(WAIST) & WAIST <  90 ~ 0L,
      SEX == "Female" & !is.na(WAIST) & WAIST >= 80 ~ 1L,
      SEX == "Female" & !is.na(WAIST) & WAIST <  80 ~ 0L,
      TRUE ~ NA_integer_
    ),

    defA_tg = case_when(
      (!is.na(TG) & TG >= 150) |
        lipid_lowering == 1L ~ 1L,
      !is.na(TG) |
        !is.na(lipid_lowering) ~ 0L,
      TRUE ~ NA_integer_
    ),

    defA_hdl = case_when(
      (SEX == "Male"   & !is.na(HDL) & HDL < 40) |
        (SEX == "Female" & !is.na(HDL) & HDL < 50) |
        lipid_lowering == 1L ~ 1L,
      !is.na(HDL) |
        !is.na(lipid_lowering) ~ 0L,
      TRUE ~ NA_integer_
    ),

    defA_bp = case_when(
      (!is.na(SBP) & SBP >= 130) |
        (!is.na(DBP) & DBP >= 85) |
        antihypertensive == 1L ~ 1L,
      !is.na(SBP) |
        !is.na(DBP) |
        !is.na(antihypertensive) ~ 0L,
      TRUE ~ NA_integer_
    ),

    defA_glucose = case_when(
      (!is.na(GLUCOSE) & GLUCOSE >= 100) |
        diabetes == 1L |
        diabetes_type1 == 1L |
        diabetes_type2 == 1L |
        glucose_lowering_oral == 1L |
        glucose_lowering_insulin == 1L |
        glucose_lowering_any == 1L ~ 1L,
      !is.na(GLUCOSE) |
        !is.na(diabetes) |
        !is.na(diabetes_type1) |
        !is.na(diabetes_type2) |
        !is.na(glucose_lowering_oral) |
        !is.na(glucose_lowering_insulin) |
        !is.na(glucose_lowering_any) ~ 0L,
      TRUE ~ NA_integer_
    ),

    defA_valid_components = rowSums(!is.na(cbind(
      defA_waist, defA_tg, defA_hdl, defA_bp, defA_glucose
    ))),

    defA_n = if_else(
      defA_valid_components == 5L,
      defA_waist + defA_tg + defA_hdl + defA_bp + defA_glucose,
      NA_integer_
    ),

    defA_muh = case_when(
      !is.na(defA_n) & defA_n <= 1L ~ 0L,
      !is.na(defA_n) & defA_n >= 2L ~ 1L,
      TRUE ~ NA_integer_
    ),

    # -------------------------------------------------------------------------
    # CARDIOMETABOLIC INDICES
    # -------------------------------------------------------------------------

    # Biochemical
    TG_HDL    = TG / HDL,
    TC_HDL    = TC / HDL,
    LDL_HDL   = LDL / HDL,
    HDL_LDL   = HDL / LDL,
    NON_HDL   = TC - HDL,
    AIP       = log10(TG / HDL),
    TyG       = if_else(!is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0,
                        log((TG * GLUCOSE) / 2), NA_real_),
    LCI       = if_else(!is.na(HDL) & HDL > 0 & !is.na(TC) & !is.na(TG) & !is.na(LDL),
                        (TC * TG * LDL) / HDL, NA_real_),

    # Anthropometric and body composition (BMI is used as recorded)
    WHtR      = if_else(!is.na(WAIST) & !is.na(HEIGHT_CM) & HEIGHT_CM > 0,
                        WAIST / HEIGHT_CM, NA_real_),
    BRI       = if_else(!is.na(WAIST) & !is.na(HEIGHT_CM) & HEIGHT_CM > 0,
                        364.2 - 365.5 * sqrt(pmax(0,
                          1 - ((WAIST / (2 * pi))^2 / ((HEIGHT_CM / 2)^2)))),
                        NA_real_),
    RFM       = case_when(
      SEX == "Male"   & !is.na(HEIGHT_M) & !is.na(WAIST) & WAIST > 0 ~
        64 - 20 * (HEIGHT_M / (WAIST / 100)),
      SEX == "Female" & !is.na(HEIGHT_M) & !is.na(WAIST) & WAIST > 0 ~
        76 - 20 * (HEIGHT_M / (WAIST / 100)),
      TRUE ~ NA_real_
    ),
    FMI       = if_else(!is.na(BFM) & !is.na(HEIGHT_M) & HEIGHT_M > 0,
                        BFM / HEIGHT_M^2, NA_real_),
    SMI       = if_else(!is.na(SMM) & !is.na(HEIGHT_M) & HEIGHT_M > 0,
                        SMM / HEIGHT_M^2, NA_real_),
    FAT_MUSCLE = if_else(!is.na(BFM) & !is.na(SMM) & SMM > 0,
                         BFM / SMM, NA_real_),
    WWI       = if_else(!is.na(WAIST) & !is.na(WEIGHT) & WEIGHT > 0,
                        WAIST / sqrt(WEIGHT), NA_real_),
    ABSI      = if_else(!is.na(BMI) & BMI > 0 & !is.na(HEIGHT_M) & HEIGHT_M > 0 & !is.na(WAIST),
                        WAIST / (BMI^(2/3) * HEIGHT_M^(1/2)), NA_real_),
    C_INDEX   = if_else(!is.na(WEIGHT) & WEIGHT > 0 & !is.na(HEIGHT_M) & HEIGHT_M > 0 & !is.na(WAIST),
                        (WAIST / 100) / (0.109 * sqrt(WEIGHT / HEIGHT_M)), NA_real_),
    # FFMI = (weight - body fat mass) / height^2 (same as Table 1)
    FFMI      = if_else(
      !is.na(WEIGHT) & !is.na(BFM) & !is.na(HEIGHT_M) & HEIGHT_M > 0 & WEIGHT > BFM,
      (WEIGHT - BFM) / HEIGHT_M^2,
      NA_real_
    ),

    # Composite
    TyG_BMI   = if_else(!is.na(TyG) & !is.na(BMI), TyG * BMI, NA_real_),
    TyG_WC    = if_else(!is.na(TyG) & !is.na(WAIST), TyG * WAIST, NA_real_),
    TyG_WHtR  = if_else(!is.na(TyG) & !is.na(WHtR), TyG * WHtR, NA_real_),
    CMI       = if_else(!is.na(WHtR) & !is.na(TG) & !is.na(HDL) & HDL > 0,
                        (TG / HDL) * WHtR, NA_real_),
    LAP       = case_when(
      SEX == "Male"   & !is.na(WAIST) & WAIST > 65 & !is.na(TG) ~ (WAIST - 65) * TG,
      SEX == "Female" & !is.na(WAIST) & WAIST > 58 & !is.na(TG) ~ (WAIST - 58) * TG,
      TRUE ~ NA_real_
    ),
    VAI       = case_when(
      SEX == "Male"   & !is.na(HDL) & HDL > 0 & !is.na(BMI) & !is.na(TG) & !is.na(WAIST) ~
        (WAIST / (39.68 + 1.88 * BMI)) * (TG / 1.03) * (1.31 / HDL),
      SEX == "Female" & !is.na(HDL) & HDL > 0 & !is.na(BMI) & !is.na(TG) & !is.na(WAIST) ~
        (WAIST / (36.58 + 1.89 * BMI)) * (TG / 0.81) * (1.52 / HDL),
      TRUE ~ NA_real_
    ),
    MCMI      = if_else(
      !is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0 &
        !is.na(HDL) & HDL > 0 & !is.na(WHtR),
      log((TG * GLUCOSE) / HDL) * WHtR, NA_real_
    ),
    METS_IR   = if_else(
      !is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0 &
        !is.na(HDL) & HDL > 0 & !is.na(BMI) & BMI > 0,
      (log((2 * GLUCOSE) + TG) * BMI) / log(HDL), NA_real_
    )
  ) %>%

  # Analytical population (same as Table 1); a non-missing BMI category
  # implies valid BMI >= 18.5 kg/m2
  filter(
    !is.na(BMI_CATEGORY),
    !is.na(SEX),
    !is.na(defA_muh),
    !is.na(defB_muh)
  )

# -- Index definitions and groups ---------------------------------------------

index_def <- tibble(
  var     = c("TG_HDL", "TC_HDL", "LDL_HDL", "HDL_LDL", "NON_HDL", "AIP", "TyG", "LCI",
              "BMI", "WHtR", "BRI", "RFM", "FMI", "SMI", "FAT_MUSCLE", "WWI", "ABSI", "C_INDEX", "FFMI",
              "TyG_BMI", "TyG_WC", "TyG_WHtR", "CMI", "LAP", "VAI", "MCMI", "METS_IR"),
  label   = c("TG/HDL ratio", "TC/HDL ratio", "LDL/HDL ratio", "HDL/LDL ratio", "Non-HDL cholesterol",
              "AIP", "TyG index", "LCI",
              "BMI", "WHtR", "BRI", "RFM", "FMI", "SMI", "Fat-to-muscle ratio", "WWI", "ABSI",
              "Conicity index", "FFMI",
              "TyG\u2013BMI", "TyG\u2013waist circumference", "TyG\u2013WHtR", "CMI", "LAP", "VAI", "MCMI", "METS-IR"),
  section = c(rep("Biochemical", 8),
              rep("Anthropometric / Body composition", 11),
              rep("Composite", 8))
)

# -- Strata and outcomes ------------------------------------------------------

strata <- list(
  list(label = "Total",       df = dat),
  list(label = "Normal weight", df = dat %>% filter(BMI_CATEGORY == "Normal weight")),
  list(label = "Overweight",    df = dat %>% filter(BMI_CATEGORY == "Overweight")),
  list(label = "Obesity",       df = dat %>% filter(BMI_CATEGORY == "Obesity"))
)

outcome_list <- list(
  list(var = "CVD", label = "Cardiovascular disease (CVD)"),
  list(var = "FLD", label = "Fatty liver disease (FLD)"),
  list(var = "CKD", label = "Chronic kidney disease (CKD)")
)

# -- AUC with 95% CI (DeLong) -------------------------------------------------

calc_auc_ci <- function(x_vec, y_vec) {
  df <- data.frame(x = x_vec, y = y_vec) %>%
    filter(!is.na(x), !is.na(y), is.finite(x))
  if (n_distinct(df$y) < 2 | nrow(df) < 20)
    return(tibble(AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_))
  roc_obj <- tryCatch(roc(df$y, df$x, quiet = TRUE, ci = TRUE), error = function(e) NULL)
  if (is.null(roc_obj))
    return(tibble(AUC = NA_real_, CI_low = NA_real_, CI_high = NA_real_))
  ci_vals <- tryCatch(ci.auc(roc_obj, method = "delong"), error = function(e) NULL)
  tibble(
    AUC     = round(as.numeric(auc(roc_obj)), 3),
    CI_low  = if (!is.null(ci_vals)) round(ci_vals[1], 3) else NA_real_,
    CI_high = if (!is.null(ci_vals)) round(ci_vals[3], 3) else NA_real_
  )
}

# -- Full AUC table -----------------------------------------------------------

auc_table <- bind_rows(lapply(strata, function(s) {
  bind_rows(lapply(outcome_list, function(o) {
    bind_rows(lapply(seq_len(nrow(index_def)), function(i) {
      res <- calc_auc_ci(s$df[[index_def$var[i]]], s$df[[o$var]])
      tibble(
        Stratum = s$label,
        Outcome = o$label,
        var     = index_def$var[i],
        label   = index_def$label[i],
        AUC     = res$AUC,
        CI_low  = res$CI_low,
        CI_high = res$CI_high
      )
    }))
  }))
})) %>%
  mutate(
    Stratum = factor(Stratum, levels = c("Total", "Normal weight", "Overweight", "Obesity"))
  )

# -- Console summaries --------------------------------------------------------

fmt_auc <- function(a, lo, hi) ifelse(is.na(a), "-", sprintf("%.3f (%.3f-%.3f)", a, lo, hi))

auc_summary <- auc_table %>%
  left_join(index_def %>% select(var, section), by = "var") %>%
  mutate(AUC_txt = fmt_auc(AUC, CI_low, CI_high))

section_order  <- c("Biochemical", "Anthropometric / Body composition", "Composite")
outcome_labels <- sapply(outcome_list, function(o) o$label)

# (1) By outcome: index by BMI category (ordered by overall AUC within group)
for (oc in outcome_labels) {
  d_ov <- auc_summary %>% filter(Outcome == oc, Stratum == "Total") %>%
    select(var, label, section, Total = AUC_txt, AUC_ov = AUC)
  d_nw <- auc_summary %>% filter(Outcome == oc, Stratum == "Normal weight") %>% select(var, NW = AUC_txt)
  d_ow <- auc_summary %>% filter(Outcome == oc, Stratum == "Overweight")    %>% select(var, OW = AUC_txt)
  d_ob <- auc_summary %>% filter(Outcome == oc, Stratum == "Obesity")       %>% select(var, OB = AUC_txt)

  wide <- d_ov %>%
    left_join(d_nw, by = "var") %>% left_join(d_ow, by = "var") %>% left_join(d_ob, by = "var") %>%
    arrange(factor(section, levels = section_order), desc(AUC_ov)) %>%
    select(Section = section, Index = label, Total,
           `Normal weight` = NW, Overweight = OW, Obesity = OB)

  cat("\n=== AUC (95% CI) of cardiometabolic indices for", oc, "by BMI category ===\n")
  print(as.data.frame(wide), row.names = FALSE)
}

# (2) Six indices with the highest overall AUC per outcome
for (oc in outcome_labels) {
  cat("\n=== Six indices with the highest overall AUC -", oc, "===\n")
  top <- auc_summary %>%
    filter(Outcome == oc, Stratum == "Total", !is.na(AUC)) %>%
    arrange(desc(AUC)) %>% slice_head(n = 6) %>%
    select(Index = label, Section = section, AUC, CI_low, CI_high)
  print(as.data.frame(top), row.names = FALSE)
}

# -- Figure 4 -----------------------------------------------------------------

stratum_colours <- c(
  "Total"       = "black",
  "Normal weight" = "#2171B5",
  "Overweight"    = "#238B45",
  "Obesity"       = "#CB181D"
)
stratum_shapes <- c(
  "Total"       = 18,   # diamond
  "Normal weight" = 16,   # circle
  "Overweight"    = 17,   # triangle
  "Obesity"       = 15    # square
)
stratum_offset <- c(
  "Total"       =  0.22,
  "Normal weight" =  0.07,
  "Overweight"    = -0.07,
  "Obesity"       = -0.22
)

# Common x-axis range (AUC) for the three panels. Applied with
# coord_cartesian: points outside the range are only clipped visually and
# CI bars are truncated at the limits; no estimate is removed.
auc_lims   <- c(0.40, 0.80)
auc_breaks <- seq(auc_lims[1], auc_lims[2], by = 0.05)

outcome_abbrev <- c(
  "Cardiovascular disease (CVD)" = "CVD",
  "Fatty liver disease (FLD)"    = "FLD",
  "Chronic kidney disease (CKD)" = "CKD"
)

# y positions: groups from top to bottom (Biochemical, Anthropometric /
# Body composition, Composite); within each group, indices ordered by
# overall AUC (highest at the top), with a gap between groups.
build_y_map <- function(df_overall, gap = 1.5) {
  section_order_y <- c("Composite", "Anthropometric / Body composition", "Biochemical")
  df_ordered <- df_overall %>%
    left_join(index_def %>% select(var, section), by = "var") %>%
    mutate(section = factor(section, levels = section_order_y)) %>%
    arrange(section, AUC)

  y_pos <- numeric(nrow(df_ordered))
  current_y <- 1
  current_section <- df_ordered$section[1]

  for (i in seq_len(nrow(df_ordered))) {
    if (df_ordered$section[i] != current_section) {
      current_y <- current_y + gap
      current_section <- df_ordered$section[i]
    }
    y_pos[i] <- current_y
    current_y <- current_y + 1
  }

  df_ordered %>%
    mutate(y_base = y_pos) %>%
    select(var, section, y_base)
}

make_figure_4 <- function(tab) {

  # Index order based on the overall AUC for CVD (same order in all panels)
  df_overall_cvd <- tab %>%
    filter(Stratum == "Total", Outcome == outcome_list[[1]]$label)
  y_map <- build_y_map(df_overall_cvd)

  # Group labels and separator lines
  section_labels <- y_map %>%
    group_by(section) %>%
    summarise(y_sec = max(y_base) + 1.15, .groups = "drop")

  separator_lines <- y_map %>%
    group_by(section) %>%
    summarise(y_max = max(y_base), .groups = "drop") %>%
    arrange(section) %>%
    filter(row_number() < n()) %>%
    mutate(y_sep = y_max + 1.45)

  df_plot <- tab %>%
    left_join(y_map, by = "var") %>%
    mutate(
      y_pos     = y_base + stratum_offset[as.character(Stratum)],
      Stratum   = factor(Stratum, levels = c("Total", "Normal weight", "Overweight", "Obesity")),
      Outcome   = factor(Outcome, levels = sapply(outcome_list, `[[`, "label")),
      CI_low_p  = pmax(CI_low,  auc_lims[1]),
      CI_high_p = pmin(CI_high, auc_lims[2])
    )

  y_labels <- df_plot %>%
    group_by(var, label) %>%
    summarise(y_mid = mean(y_pos), .groups = "drop")

  # Group labels only in the first panel
  first_outcome <- levels(df_plot$Outcome)[1]
  section_text_df <- section_labels %>%
    mutate(
      Outcome = factor(first_outcome, levels = levels(df_plot$Outcome)),
      label   = as.character(section)
    )

  ggplot(df_plot, aes(x = AUC, y = y_pos,
                      color = Stratum, shape = Stratum)) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               color = "gray60", linewidth = 0.4) +
    geom_hline(data = separator_lines, aes(yintercept = y_sep),
               color = "gray75", linewidth = 0.35, linetype = "solid",
               inherit.aes = FALSE) +
    geom_errorbarh(aes(xmin = CI_low_p, xmax = CI_high_p),
                   height = 0.08, linewidth = 0.4, na.rm = TRUE) +
    geom_point(aes(size = Stratum), na.rm = TRUE) +
    geom_text(data = section_text_df,
              aes(y = y_sec, label = label), x = -Inf,
              hjust = 1, size = 3.0, fontface = "bold",
              color = "gray30", inherit.aes = FALSE) +
    facet_wrap(~ Outcome, nrow = 1, scales = "fixed",
               strip.position = "bottom",
               labeller = as_labeller(outcome_abbrev)) +
    scale_size_manual(values = c("Total" = 2.5, "Normal weight" = 1.8,
                                 "Overweight" = 1.8, "Obesity" = 1.8),
                      guide = "none") +
    scale_color_manual(values = stratum_colours, name = "BMI category") +
    scale_shape_manual(values = stratum_shapes,  name = "BMI category") +
    scale_x_continuous(
      breaks = auc_breaks,
      labels = function(x) sprintf("%.2f", x)
    ) +
    coord_cartesian(xlim = auc_lims, clip = "off") +
    scale_y_continuous(
      breaks = y_labels$y_mid,
      labels = y_labels$label
    ) +
    labs(
      title = NULL,
      x     = "AUC (95% CI)",
      y     = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title         = element_blank(),
      strip.text         = element_text(face = "bold", size = 9, color = "gray20"),
      strip.background   = element_rect(fill = "gray95", color = NA),
      strip.placement    = "outside",
      axis.text.x        = element_text(size = 7,  color = "gray20", angle = 30, hjust = 1),
      axis.text.y        = element_text(size = 7.5, color = "gray20"),
      axis.title.x       = element_text(size = 8,  margin = margin(t = 6)),
      legend.title       = element_text(face = "bold", size = 10.4),
      legend.text        = element_text(size = 9.8),
      legend.position    = "bottom",
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.margin        = margin(15, 15, 15, 35)
    )
}

figure_4 <- make_figure_4(auc_table)
print(figure_4)
