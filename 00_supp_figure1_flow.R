###############################################################################
# Supplementary Figure 1. Selection of participants (flow diagram counts)
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Exclusion criteria are applied sequentially, in the order shown in
# Supplementary Figure 1:
#   1. Missing or implausible BMI (outside 12-70 kg/m2)
#   2. Underweight (BMI < 18.5 kg/m2)
#   3. Sex not recorded as female or male
#   4. Metabolic health status not classifiable under Definition A
#   5. Metabolic health status not classifiable under Definition B
#
# Output (printed to the console; no files are written):
#   - Number of participants at each step of the flow diagram
#   - Missing components among participants not classifiable under each
#     definition (categories not mutually exclusive)
#   - Analytical sample by BMI category
#   - Participants with available data for each outcome
#
# Data cleaning, physiological ranges, and definitions are identical to
# those in 01_table1_characteristics.R.
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs; MUH = 2 or more of 5 components
# Definition B: empirically derived definition (Zembic and colleagues);
#               MUH = 1 or more of 3 components
#
# R version 4.6.0
# Packages: dplyr, readr, stringr, tibble
###############################################################################

library(dplyr)
library(readr)
library(stringr)
library(tibble)

# Path to the oriGen questionnaire data file (tab-separated).
# Column names correspond to the original oriGen data release.
data_path <- "path/to/oriGen_questionnaire_data.tsv"
raw_data  <- read_tsv(data_path, show_col_types = FALSE)

# == Helper functions ==========================================================

# Harmonises free-text responses: upper case, no accents, single spaces.
normalise_text <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- trimws(x)
  x <- str_to_upper(x)
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  gsub("[[:space:]]+", " ", x)
}

# Converts yes/no responses (Spanish or English) to 1/0. normalise_text()
# removes accents, so "SI" also covers the accented variant.
safe_bool <- function(x) {
  y <- normalise_text(x)
  case_when(
    y %in% c("TRUE", "VERDADERO", "SI", "1", "YES") ~ 1L,
    y %in% c("FALSE", "FALSO", "NO", "0")           ~ 0L,
    TRUE                                            ~ NA_integer_
  )
}

# Values outside the prespecified physiological range [lo, hi] are set to NA.
clip_range <- function(x, lo, hi) if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

# == Data preparation (same as 01_table1_characteristics.R) ====================

dat_all <- raw_data %>%
  transmute(

    SEX_RAW = str_to_upper(trimws(as.character(PARTICIPANTE.SEXO))),
    Sex = case_when(
      SEX_RAW %in% c("H", "HOMBRE", "MASCULINO", "MALE", "MAN")         ~ "Male",
      SEX_RAW %in% c("M", "MUJER", "F", "FEMENINO", "FEMALE", "WOMAN") ~ "Female",
      TRUE ~ NA_character_
    ),

    BMI_RAW   = suppressWarnings(as.numeric(ANTROPOMETRIA.IMC)),
    WAIST_RAW = suppressWarnings(as.numeric(ANTROPOMETRIA.CIRCUNFERENCIA_CINTURA)),
    WHR_RAW   = suppressWarnings(as.numeric(INBODY.WHR_Waist_Hip_Ratio)),

    FASTING_HOURS = suppressWarnings(as.numeric(ANTROPOMETRIA.HORAS_AYUNO)),
    GLUCOSE_RAW   = suppressWarnings(as.numeric(METABOLICOS.GLUCOSA)),
    TG_RAW        = suppressWarnings(as.numeric(METABOLICOS.TRIGLICERIDOS)),
    HDL_RAW       = suppressWarnings(as.numeric(METABOLICOS.HDL)),
    SBP_RAW       = suppressWarnings(as.numeric(METABOLICOS.SISTOLICA)),
    DBP_RAW       = suppressWarnings(as.numeric(METABOLICOS.DIASTOLICA)),

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
    FLD            = safe_bool(HISTORIA_MEDICA.HIGADO_GRASO),
    ckd_diagnosis  = safe_bool(HISTORIA_MEDICA.ENFERMEDAD_RENAL_CRONICA),
    ckd_renal_insf = safe_bool(HISTORIA_MEDICA.INSUFICIENCIA_RENAL),
    ckd_dialysis   = safe_bool(HISTORIA_MEDICA.DIALISIS_POR_RENAL),
    ckd_haemodial  = safe_bool(HISTORIA_MEDICA.HEMODIALISIS)

  ) %>%
  # Prespecified physiological ranges (same as the other scripts)
  mutate(
    BMI         = clip_range(BMI_RAW,       12,   70),
    WAIST       = clip_range(WAIST_RAW,     50,  200),
    WHR         = clip_range(WHR_RAW,      0.6,  1.5),
    GLUCOSE_RAW = clip_range(GLUCOSE_RAW,   40,  600),
    TG          = clip_range(TG_RAW,        20, 2000),
    HDL         = clip_range(HDL_RAW,       10,  150),
    SBP         = clip_range(SBP_RAW,       70,  260),
    DBP         = clip_range(DBP_RAW,       40,  160)
  ) %>%
  mutate(

    # Glucose is considered a fasting value only after >= 8 h of fasting.
    GLUCOSE = if_else(FASTING_HOURS >= 8, GLUCOSE_RAW, NA_real_),

    # Composite outcomes: 1 if any item is positive; 0 if none is positive
    # and at least one item is evaluable; NA only if all items are missing.
    CVD = case_when(
      cvd_mi == 1L | cvd_stroke == 1L | cvd_ischaemic == 1L |
        cvd_pad == 1L | cvd_athero == 1L | cvd_hf == 1L ~ 1L,
      !is.na(cvd_mi) | !is.na(cvd_stroke) | !is.na(cvd_ischaemic) |
        !is.na(cvd_pad) | !is.na(cvd_athero) | !is.na(cvd_hf) ~ 0L,
      TRUE ~ NA_integer_
    ),
    CKD = case_when(
      ckd_diagnosis == 1L | ckd_renal_insf == 1L |
        ckd_dialysis == 1L | ckd_haemodial == 1L ~ 1L,
      !is.na(ckd_diagnosis) | !is.na(ckd_renal_insf) |
        !is.na(ckd_dialysis) | !is.na(ckd_haemodial) ~ 0L,
      TRUE ~ NA_integer_
    ),

    BMI_CATEGORY = case_when(
      !is.na(BMI) & BMI >= 18.5 & BMI < 25 ~ "Normal weight",
      !is.na(BMI) & BMI >= 25   & BMI < 30 ~ "Overweight",
      !is.na(BMI) & BMI >= 30              ~ "Obesity",
      TRUE ~ NA_character_
    ),

    # -------------------------------------------------------------------------
    # METABOLIC HEALTH DEFINITIONS (same coding as 01_table1_characteristics.R)
    #   - 1 if any of its evaluable criteria is met;
    #   - 0 if none is met and at least one criterion is evaluable;
    #   - NA only if all criteria are missing, uninterpretable, or out of range.
    # A definition is resolvable only when all of its components are non-missing.
    # -------------------------------------------------------------------------

    # ============================ DEFINITION B ==============================
    defB_diabetes = case_when(
      diabetes == 1L ~ 1L,
      diabetes == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),
    defB_bp = case_when(
      (!is.na(SBP) & SBP >= 130) | antihypertensive == 1L ~ 1L,
      !is.na(SBP) | !is.na(antihypertensive)              ~ 0L,
      TRUE ~ NA_integer_
    ),
    defB_whr = case_when(
      Sex == "Female" & !is.na(WHR) & WHR >= 0.95 ~ 1L,
      Sex == "Female" & !is.na(WHR) & WHR <  0.95 ~ 0L,
      Sex == "Male"   & !is.na(WHR) & WHR >= 1.03 ~ 1L,
      Sex == "Male"   & !is.na(WHR) & WHR <  1.03 ~ 0L,
      TRUE ~ NA_integer_
    ),
    defB_resolvable = !is.na(defB_diabetes) & !is.na(defB_bp) & !is.na(defB_whr),

    # ============================ DEFINITION A ==============================
    defA_waist = case_when(
      Sex == "Male"   & !is.na(WAIST) & WAIST >= 90 ~ 1L,
      Sex == "Male"   & !is.na(WAIST) & WAIST <  90 ~ 0L,
      Sex == "Female" & !is.na(WAIST) & WAIST >= 80 ~ 1L,
      Sex == "Female" & !is.na(WAIST) & WAIST <  80 ~ 0L,
      TRUE ~ NA_integer_
    ),
    defA_tg = case_when(
      (!is.na(TG) & TG >= 150) | lipid_lowering == 1L ~ 1L,
      !is.na(TG) | !is.na(lipid_lowering)             ~ 0L,
      TRUE ~ NA_integer_
    ),
    defA_hdl = case_when(
      (Sex == "Male"   & !is.na(HDL) & HDL < 40) |
        (Sex == "Female" & !is.na(HDL) & HDL < 50) |
        lipid_lowering == 1L ~ 1L,
      !is.na(HDL) | !is.na(lipid_lowering) ~ 0L,
      TRUE ~ NA_integer_
    ),
    defA_bp = case_when(
      (!is.na(SBP) & SBP >= 130) | (!is.na(DBP) & DBP >= 85) |
        antihypertensive == 1L ~ 1L,
      !is.na(SBP) | !is.na(DBP) | !is.na(antihypertensive) ~ 0L,
      TRUE ~ NA_integer_
    ),
    defA_glucose = case_when(
      (!is.na(GLUCOSE) & GLUCOSE >= 100) |
        diabetes == 1L | diabetes_type1 == 1L | diabetes_type2 == 1L |
        glucose_lowering_oral == 1L | glucose_lowering_insulin == 1L |
        glucose_lowering_any == 1L ~ 1L,
      !is.na(GLUCOSE) | !is.na(diabetes) | !is.na(diabetes_type1) |
        !is.na(diabetes_type2) | !is.na(glucose_lowering_oral) |
        !is.na(glucose_lowering_insulin) | !is.na(glucose_lowering_any) ~ 0L,
      TRUE ~ NA_integer_
    ),
    defA_resolvable = !is.na(defA_waist) & !is.na(defA_tg) & !is.na(defA_hdl) &
      !is.na(defA_bp) & !is.na(defA_glucose)
  )

# == Sequential exclusions =====================================================

n_total <- nrow(dat_all)

step1 <- dat_all %>% filter(!is.na(BMI))       # valid BMI (12-70 kg/m2)
step2 <- step1   %>% filter(BMI >= 18.5)       # not underweight
step3 <- step2   %>% filter(!is.na(Sex))       # sex recorded as female or male
step4 <- step3   %>% filter(defA_resolvable)   # classifiable under Definition A
dat   <- step4   %>% filter(defB_resolvable)   # classifiable under Definition B

flow <- tibble(
  Step = c(
    "Participants in the oriGen data release",
    "Excluded: missing or implausible BMI (outside 12-70 kg/m2)",
    "Excluded: underweight (BMI < 18.5 kg/m2)",
    "Excluded: sex not recorded as female or male",
    "Excluded: metabolic health status not classifiable under Definition A",
    "Excluded: metabolic health status not classifiable under Definition B",
    "Total excluded",
    "Analytical sample"
  ),
  N = c(
    n_total,
    n_total     - nrow(step1),
    nrow(step1) - nrow(step2),
    nrow(step2) - nrow(step3),
    nrow(step3) - nrow(step4),
    nrow(step4) - nrow(dat),
    n_total     - nrow(dat),
    nrow(dat)
  )
) %>%
  mutate(N = fmt_n(N))

cat("\n=== Supplementary Figure 1: selection of participants ===\n")
cat("(exclusion criteria applied sequentially)\n\n")
print(as.data.frame(flow), row.names = FALSE, right = FALSE)

# Internal check: the steps must add up to the analytical sample.
stopifnot(
  n_total - (n_total - nrow(step1)) - (nrow(step1) - nrow(step2)) -
    (nrow(step2) - nrow(step3)) - (nrow(step3) - nrow(step4)) -
    (nrow(step4) - nrow(dat)) == nrow(dat)
)

# == Missing components among non-classifiable participants ====================
# Categories are not mutually exclusive.

not_A <- step3 %>% filter(!defA_resolvable)

missing_A <- tibble(
  Component = c(
    "Triglycerides",
    "HDL cholesterol",
    "Waist circumference",
    "Blood pressure",
    "Glucose",
    "Two or more components missing"
  ),
  N = c(
    sum(is.na(not_A$defA_tg)),
    sum(is.na(not_A$defA_hdl)),
    sum(is.na(not_A$defA_waist)),
    sum(is.na(not_A$defA_bp)),
    sum(is.na(not_A$defA_glucose)),
    sum(rowSums(is.na(cbind(not_A$defA_waist, not_A$defA_tg, not_A$defA_hdl,
                            not_A$defA_bp, not_A$defA_glucose))) >= 2)
  )
) %>%
  mutate(N = fmt_n(N))

cat("\n=== Missing components: not classifiable under Definition A (n = ",
    fmt_n(nrow(not_A)), ") ===\n", sep = "")
print(as.data.frame(missing_A), row.names = FALSE, right = FALSE)

not_B <- step4 %>% filter(!defB_resolvable)

missing_B <- tibble(
  Component = c(
    "Blood pressure",
    "Waist-to-hip ratio, outside the plausible range (0.6-1.5)",
    "Waist-to-hip ratio, not recorded",
    "Diabetes",
    "Two or more components missing"
  ),
  N = c(
    sum(is.na(not_B$defB_bp)),
    sum(is.na(not_B$defB_whr) & !is.na(not_B$WHR_RAW)),
    sum(is.na(not_B$defB_whr) &  is.na(not_B$WHR_RAW)),
    sum(is.na(not_B$defB_diabetes)),
    sum(rowSums(is.na(cbind(not_B$defB_diabetes, not_B$defB_bp,
                            not_B$defB_whr))) >= 2)
  )
) %>%
  mutate(N = fmt_n(N))

cat("\n=== Missing components: not classifiable under Definition B (n = ",
    fmt_n(nrow(not_B)), ") ===\n", sep = "")
print(as.data.frame(missing_B), row.names = FALSE, right = FALSE)

# == Analytical sample by BMI category =========================================

bmi_counts <- dat %>%
  count(BMI_CATEGORY) %>%
  mutate(
    BMI_CATEGORY = factor(BMI_CATEGORY,
                          levels = c("Normal weight", "Overweight", "Obesity")),
    Percent = sprintf("%.1f%%", 100 * n / sum(n))
  ) %>%
  arrange(BMI_CATEGORY) %>%
  transmute(`BMI category` = as.character(BMI_CATEGORY),
            N = fmt_n(n), Percent)

cat("\n=== Analytical sample by BMI category ===\n")
print(as.data.frame(bmi_counts), row.names = FALSE, right = FALSE)

stopifnot(sum(!is.na(dat$BMI_CATEGORY)) == nrow(dat))

# == Outcome data availability in the analytical sample ========================

any_missing_outcome <- is.na(dat$CVD) | is.na(dat$FLD) | is.na(dat$CKD)

outcome_counts <- tibble(
  Item = c(
    "CVD, participants with available data",
    "FLD, participants with available data",
    "CKD, participants with available data",
    "Missing data for at least one outcome (retained in the sample)"
  ),
  N = c(
    sum(!is.na(dat$CVD)),
    sum(!is.na(dat$FLD)),
    sum(!is.na(dat$CKD)),
    sum(any_missing_outcome)
  )
) %>%
  mutate(
    Percent = c(rep("", 3),
                sprintf("%.1f%%", 100 * sum(any_missing_outcome) / nrow(dat))),
    N = fmt_n(N)
  )

cat("\n=== Outcome data in the analytical sample (n = ",
    fmt_n(nrow(dat)), ") ===\n", sep = "")
print(as.data.frame(outcome_counts), row.names = FALSE, right = FALSE)
