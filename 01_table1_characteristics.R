###############################################################################
# Table 1. Characteristics of the study population, overall and by
# metabolic health status under Definition A and Definition B
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Columns: Variable | Overall | MH A | MUH A | p* | MH B | MUH B | p** |
#          p*** | p****
#   p*    : MH vs MUH within Definition A
#   p**   : MH vs MUH within Definition B
#   p***  : MH under Definition A vs MH under Definition B
#   p**** : MUH under Definition A vs MUH under Definition B
# Because both definitions were applied to the same participants, p*** and
# p**** compare overlapping groups and are interpreted descriptively.
#
# Tests: Kruskal-Wallis test (continuous variables) and chi-squared test
# (categorical variables). Data are mean (SD) or n (%).
#
# Output (printed to the console; no files are written): Table 1 and its
# footnote.
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs; MUH = 2 or more of 5 components
# Definition B: empirically derived definition (Zembic et al.);
#               MUH = 1 or more of 3 components
#
# Analytical population: BMI >= 18.5 kg/m2, sex recorded as female or male,
# and metabolic health status resolvable under both definitions. Each
# variable is summarised in participants with available data.
#
# Non-ASCII characters in table labels are written as Unicode escapes
# (\uXXXX) so that this file remains pure ASCII.
#
# R version 4.6.0
# Packages: dplyr, readr, stringr, tidyr, tibble
###############################################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)

# Path to the oriGen questionnaire data file (tab-separated).
# Column names correspond to the original oriGen data release.
data_path <- "path/to/oriGen_questionnaire_data.tsv"
raw_data  <- read_tsv(data_path, show_col_types = FALSE)

# Reference year used to estimate age from year of birth
REFERENCE_YEAR <- 2025L

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

fmt_mean_sd <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  paste0(round(mean(x), 1), " (", round(sd(x), 1), ")")
}

fmt_n_pct <- function(x) {
  n   <- sum(x == 1L, na.rm = TRUE)
  tot <- sum(!is.na(x))
  if (tot == 0) return(NA_character_)
  paste0(n, " (", round(100 * n / tot, 1), "%)")
}

# MH vs MUH within one definition
p_continuous <- function(x, grp) {
  df <- data.frame(x = x, g = grp) %>% filter(!is.na(x), is.finite(x), !is.na(g))
  if (nrow(df) < 5 || n_distinct(df$g) < 2) return(NA_real_)
  tryCatch(kruskal.test(x ~ g, data = df)$p.value, error = function(e) NA_real_)
}

p_categorical <- function(x, grp) {
  df <- data.frame(x = x, g = grp) %>% filter(!is.na(x), !is.na(g))
  if (nrow(df) < 5) return(NA_real_)
  tt <- table(df$x, df$g)
  if (nrow(tt) < 2 || ncol(tt) < 2) return(NA_real_)
  tryCatch(chisq.test(tt)$p.value, error = function(e) NA_real_)
}

# Comparison between two subgroups across definitions (eg, MH under
# Definition A vs MH under Definition B). The subgroups overlap (same
# participants classified by two definitions); groups are treated as
# independent and the comparison is descriptive, not paired.
p_continuous_groups <- function(x_a, x_b) {
  df <- data.frame(
    x = c(x_a, x_b),
    g = c(rep("A", length(x_a)), rep("B", length(x_b)))
  ) %>% filter(!is.na(x), is.finite(x))
  if (nrow(df) < 5 || n_distinct(df$g) < 2) return(NA_real_)
  tryCatch(kruskal.test(x ~ g, data = df)$p.value, error = function(e) NA_real_)
}

p_categorical_groups <- function(x_a, x_b) {
  df <- data.frame(
    x = c(x_a, x_b),
    g = c(rep("A", length(x_a)), rep("B", length(x_b)))
  ) %>% filter(!is.na(x))
  if (nrow(df) < 5 || n_distinct(df$g) < 2) return(NA_real_)
  tt <- table(df$x, df$g)
  if (nrow(tt) < 2) return(NA_real_)
  tryCatch(chisq.test(tt)$p.value, error = function(e) NA_real_)
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  as.character(round(p, 3))
}

# == Data preparation ==========================================================

dat_unfiltered <- raw_data %>%
  transmute(

    # Age sources
    AGE_REGISTRATION_RAW = suppressWarnings(as.numeric(PARTICIPANTE.EDAD_REGISTRO)),
    BIRTH_YEAR           = suppressWarnings(as.numeric(PARTICIPANTE.YYYY_NACIMIENTO)),
    AGE_INBODY_RAW       = suppressWarnings(as.numeric(INBODY.Age)),

    SEX_RAW = str_to_upper(trimws(as.character(PARTICIPANTE.SEXO))),
    Sex = case_when(
      SEX_RAW %in% c("H", "HOMBRE", "MASCULINO", "MALE", "MAN")         ~ "Male",
      SEX_RAW %in% c("M", "MUJER", "F", "FEMENINO", "FEMALE", "WOMAN") ~ "Female",
      TRUE ~ NA_character_
    ),

    BMI        = suppressWarnings(as.numeric(ANTROPOMETRIA.IMC)),
    WEIGHT     = suppressWarnings(as.numeric(INBODY.Weight)),
    WAIST      = suppressWarnings(as.numeric(ANTROPOMETRIA.CIRCUNFERENCIA_CINTURA)),
    WHR        = suppressWarnings(as.numeric(INBODY.WHR_Waist_Hip_Ratio)),
    SMM        = suppressWarnings(as.numeric(INBODY.SMM_Skeletal_Muscle_Mass)),
    BFM        = suppressWarnings(as.numeric(INBODY.BFM_Body_Fat_Mass)),
    PBF        = suppressWarnings(as.numeric(INBODY.PBF_Percent_Body_Fat)),
    HEART_RATE = suppressWarnings(as.numeric(METABOLICOS.LATIDOS_POR_MIN)),

    HEIGHT_RAW = suppressWarnings(as.numeric(ANTROPOMETRIA.ESTATURA)),
    # Height recorded in cm or m: values > 3 are taken as cm
    HEIGHT_CM  = if_else(HEIGHT_RAW > 3, HEIGHT_RAW, HEIGHT_RAW * 100),
    HEIGHT_M   = if_else(HEIGHT_RAW > 3, HEIGHT_RAW / 100, HEIGHT_RAW),

    FASTING_HOURS = suppressWarnings(as.numeric(ANTROPOMETRIA.HORAS_AYUNO)),
    GLUCOSE_RAW   = suppressWarnings(as.numeric(METABOLICOS.GLUCOSA)),
    TG            = suppressWarnings(as.numeric(METABOLICOS.TRIGLICERIDOS)),
    TC            = suppressWarnings(as.numeric(METABOLICOS.COLESTEROL)),
    HDL           = suppressWarnings(as.numeric(METABOLICOS.HDL)),
    LDL           = suppressWarnings(as.numeric(METABOLICOS.LDL)),
    SBP           = suppressWarnings(as.numeric(METABOLICOS.SISTOLICA)),
    DBP           = suppressWarnings(as.numeric(METABOLICOS.DIASTOLICA)),

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
  # Prespecified physiological ranges, including the three age sources
  mutate(
    WEIGHT           = clip_range(WEIGHT,       30,  300),
    HEIGHT_CM        = clip_range(HEIGHT_CM,   130,  220),
    HEIGHT_M         = clip_range(HEIGHT_M,   1.30, 2.20),
    BMI              = clip_range(BMI,          12,   70),
    WAIST            = clip_range(WAIST,        50,  200),
    WHR              = clip_range(WHR,         0.6,  1.5),
    GLUCOSE_RAW      = clip_range(GLUCOSE_RAW,  40,  600),
    TG               = clip_range(TG,           20, 2000),
    TC               = clip_range(TC,           50,  500),
    HDL              = clip_range(HDL,          10,  150),
    LDL              = clip_range(LDL,          10,  400),
    SBP              = clip_range(SBP,          70,  260),
    DBP              = clip_range(DBP,          40,  160),
    PBF              = clip_range(PBF,           3,   70),
    AGE_REGISTRATION = clip_range(AGE_REGISTRATION_RAW, 18, 110),
    BIRTH_YEAR       = clip_range(BIRTH_YEAR, 1900, REFERENCE_YEAR),
    AGE_BIRTH_YEAR   = clip_range(REFERENCE_YEAR - BIRTH_YEAR, 18, 110),
    AGE_INBODY       = clip_range(AGE_INBODY_RAW, 18, 110),
    # Age: age at registration; if missing, age estimated from year of
    # birth; if missing, age recorded by the bioimpedance analyser
    Age              = coalesce(AGE_REGISTRATION, AGE_BIRTH_YEAR, AGE_INBODY),
    SMM              = clip_range(SMM, 5,  70),
    BFM              = clip_range(BFM, 2, 120)
  ) %>%
  mutate(

    # Glucose is considered a fasting value only after >= 8 h of fasting.
    GLUCOSE = if_else(FASTING_HOURS >= 8, GLUCOSE_RAW, NA_real_),
    WHtR    = if_else(!is.na(WAIST) & !is.na(HEIGHT_CM) & HEIGHT_CM >= 130 & HEIGHT_CM <= 220,
                      WAIST / HEIGHT_CM, NA_real_),
    SMI     = if_else(!is.na(SMM) & !is.na(HEIGHT_M) & HEIGHT_M >= 1.30 & HEIGHT_M <= 2.20,
                      SMM / HEIGHT_M^2, NA_real_),
    FMI     = if_else(!is.na(BFM) & !is.na(HEIGHT_M) & HEIGHT_M >= 1.30 & HEIGHT_M <= 2.20,
                      BFM / HEIGHT_M^2, NA_real_),

    TyG = if_else(
      !is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0,
      log((TG * GLUCOSE) / 2), NA_real_
    ),
    VAI = case_when(
      Sex == "Male" & !is.na(HDL) & HDL > 0 &
        !is.na(BMI) & !is.na(TG) & !is.na(WAIST) ~
        (WAIST / (39.68 + 1.88 * BMI)) * (TG / 1.03) * (1.31 / HDL),
      Sex == "Female" & !is.na(HDL) & HDL > 0 &
        !is.na(BMI) & !is.na(TG) & !is.na(WAIST) ~
        (WAIST / (36.58 + 1.89 * BMI)) * (TG / 0.81) * (1.52 / HDL),
      TRUE ~ NA_real_
    ),
    TyG_WC   = if_else(!is.na(TyG) & !is.na(WAIST), TyG * WAIST, NA_real_),
    TyG_WHtR = if_else(!is.na(TyG) & !is.na(WHtR),  TyG * WHtR,  NA_real_),
    MCMI     = if_else(
      !is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0 &
        !is.na(HDL) & HDL > 0 & !is.na(WHtR),
      log((TG * GLUCOSE) / HDL) * WHtR, NA_real_
    ),
    FAT_MUSCLE = if_else(!is.na(BFM) & !is.na(SMM) & SMM > 0, BFM / SMM, NA_real_),
    WWI        = if_else(!is.na(WAIST) & !is.na(WEIGHT) & WEIGHT > 0,
                         WAIST / sqrt(WEIGHT), NA_real_),
    BRI        = if_else(!is.na(WAIST) & !is.na(HEIGHT_CM) & HEIGHT_CM > 0,
                         364.2 - 365.5 * sqrt(pmax(0,
                           1 - ((WAIST / (2 * pi))^2 / ((HEIGHT_CM / 2)^2)))),
                         NA_real_),
    ABSI       = if_else(!is.na(BMI) & BMI > 0 & !is.na(HEIGHT_M) & HEIGHT_M > 0 & !is.na(WAIST),
                         WAIST / (BMI^(2/3) * HEIGHT_M^(1/2)), NA_real_),

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

    SEX_MALE = if_else(Sex == "Male", 1L, 0L),

    # BMI categories (BMI < 18.5 kg/m2 is excluded from the analytical population)
    NORMAL_WEIGHT = if_else(!is.na(BMI) & BMI < 25,             1L, 0L),
    OVERWEIGHT    = if_else(!is.na(BMI) & BMI >= 25 & BMI < 30, 1L, 0L),
    OBESITY       = if_else(!is.na(BMI) & BMI >= 30,            1L, 0L),

    # Age groups
    AGE_LT40  = if_else(!is.na(Age) & Age < 40,             1L, 0L),
    AGE_40_59 = if_else(!is.na(Age) & Age >= 40 & Age < 60, 1L, 0L),
    AGE_GE60  = if_else(!is.na(Age) & Age >= 60,            1L, 0L),

    # -------------------------------------------------------------------------
    # METABOLIC HEALTH DEFINITIONS
    #
    # Component coding: 1 = present; 0 = absent.
    #   - 1 if any of its evaluable criteria is met (measured value meets the
    #     cut-off, or the corresponding diagnosis or treatment is reported);
    #   - 0 if none is met and at least one criterion is evaluable;
    #   - NA only if all criteria are missing, uninterpretable, or out of range.
    # A definition is resolvable only when all of its components are non-missing.
    # -------------------------------------------------------------------------

    # ============================ DEFINITION B ==============================
    # MUH if 1 or more of 3 components are present.

    # B1) Self-reported diabetes
    defB_diabetes = case_when(
      diabetes == 1L ~ 1L,
      diabetes == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),
    # B2) Systolic blood pressure >= 130 mm Hg or antihypertensive treatment
    defB_bp = case_when(
      (!is.na(SBP) & SBP >= 130) | antihypertensive == 1L ~ 1L,
      !is.na(SBP) | !is.na(antihypertensive)              ~ 0L,
      TRUE ~ NA_integer_
    ),
    # B3) Waist-to-hip ratio >= 0.95 in women or >= 1.03 in men
    defB_whr = case_when(
      Sex == "Female" & !is.na(WHR) & WHR >= 0.95 ~ 1L,
      Sex == "Female" & !is.na(WHR) & WHR <  0.95 ~ 0L,
      Sex == "Male"   & !is.na(WHR) & WHR >= 1.03 ~ 1L,
      Sex == "Male"   & !is.na(WHR) & WHR <  1.03 ~ 0L,
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
    defB_status = case_when(
      !is.na(defB_n) & defB_n == 0L ~ "MH",
      !is.na(defB_n) & defB_n >= 1L ~ "MUH",
      TRUE ~ NA_character_
    ),

    # ============================ DEFINITION A ==============================
    # MH = 0-1 components; MUH = 2 or more of 5 components.

    # A1) Waist circumference >= 90 cm in men or >= 80 cm in women
    defA_waist = case_when(
      Sex == "Male"   & !is.na(WAIST) & WAIST >= 90 ~ 1L,
      Sex == "Male"   & !is.na(WAIST) & WAIST <  90 ~ 0L,
      Sex == "Female" & !is.na(WAIST) & WAIST >= 80 ~ 1L,
      Sex == "Female" & !is.na(WAIST) & WAIST <  80 ~ 0L,
      TRUE ~ NA_integer_
    ),
    # A2) Triglycerides >= 150 mg/dL or lipid-lowering treatment
    defA_tg = case_when(
      (!is.na(TG) & TG >= 150) | lipid_lowering == 1L ~ 1L,
      !is.na(TG) | !is.na(lipid_lowering)             ~ 0L,
      TRUE ~ NA_integer_
    ),
    # A3) HDL cholesterol < 40 mg/dL in men or < 50 mg/dL in women,
    #     or lipid-lowering treatment
    defA_hdl = case_when(
      (Sex == "Male"   & !is.na(HDL) & HDL < 40) |
        (Sex == "Female" & !is.na(HDL) & HDL < 50) |
        lipid_lowering == 1L ~ 1L,
      !is.na(HDL) | !is.na(lipid_lowering) ~ 0L,
      TRUE ~ NA_integer_
    ),
    # A4) Systolic >= 130 mm Hg, diastolic >= 85 mm Hg, or antihypertensive
    #     treatment
    defA_bp = case_when(
      (!is.na(SBP) & SBP >= 130) | (!is.na(DBP) & DBP >= 85) |
        antihypertensive == 1L ~ 1L,
      !is.na(SBP) | !is.na(DBP) | !is.na(antihypertensive) ~ 0L,
      TRUE ~ NA_integer_
    ),
    # A5) Fasting glucose >= 100 mg/dL, self-reported diabetes (any type),
    #     or glucose-lowering treatment (oral agents, insulin, or any)
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
    defA_valid_components = rowSums(!is.na(cbind(
      defA_waist, defA_tg, defA_hdl, defA_bp, defA_glucose
    ))),
    defA_n = if_else(
      defA_valid_components == 5L,
      defA_waist + defA_tg + defA_hdl + defA_bp + defA_glucose,
      NA_integer_
    ),
    defA_status = case_when(
      !is.na(defA_n) & defA_n <= 1L ~ "MH",
      !is.na(defA_n) & defA_n >= 2L ~ "MUH",
      TRUE ~ NA_character_
    )
  )

# == Analytical population =====================================================
# Valid BMI >= 18.5 kg/m2, sex recorded as female or male, and metabolic
# health status resolvable under both definitions.

dat <- dat_unfiltered %>%
  filter(
    !is.na(BMI),
    BMI >= 18.5,
    !is.na(Sex),
    !is.na(defA_status),
    !is.na(defB_status)
  )

# == Additional indices for Table 1 ============================================
# Lipids in mg/dL; waist circumference and height in cm.

dat <- dat %>%
  mutate(
    # Biochemical
    TG_HDL  = if_else(!is.na(TG) & !is.na(HDL) & HDL > 0, TG / HDL, NA_real_),
    TC_HDL  = if_else(!is.na(TC) & !is.na(HDL) & HDL > 0, TC / HDL, NA_real_),
    LDL_HDL = if_else(!is.na(LDL) & !is.na(HDL) & HDL > 0, LDL / HDL, NA_real_),
    HDL_LDL = if_else(!is.na(HDL) & !is.na(LDL) & LDL > 0, HDL / LDL, NA_real_),
    NON_HDL = if_else(!is.na(TC) & !is.na(HDL), TC - HDL, NA_real_),
    AIP     = if_else(!is.na(TG_HDL) & TG_HDL > 0, log10(TG_HDL), NA_real_),
    LCI     = if_else(!is.na(TC) & !is.na(TG) & !is.na(LDL) &
                        !is.na(HDL) & HDL > 0,
                      TC * TG * LDL / HDL, NA_real_),

    # Anthropometric and body composition
    RFM = case_when(
      Sex == "Male"   & !is.na(HEIGHT_CM) & !is.na(WAIST) & WAIST > 0 ~
        64 - 20 * (HEIGHT_CM / WAIST),
      Sex == "Female" & !is.na(HEIGHT_CM) & !is.na(WAIST) & WAIST > 0 ~
        76 - 20 * (HEIGHT_CM / WAIST),
      TRUE ~ NA_real_
    ),
    # FFMI = (weight - body fat mass) / height^2
    FFMI    = if_else(!is.na(WEIGHT) & !is.na(BFM) & !is.na(HEIGHT_M) & WEIGHT > BFM,
                      (WEIGHT - BFM) / HEIGHT_M^2, NA_real_),
    C_INDEX = if_else(!is.na(WAIST) & !is.na(WEIGHT) & WEIGHT > 0 &
                        !is.na(HEIGHT_M) & HEIGHT_M > 0,
                      (WAIST / 100) / (0.109 * sqrt(WEIGHT / HEIGHT_M)), NA_real_),

    # Composite
    TyG_BMI = if_else(!is.na(TyG) & !is.na(BMI), TyG * BMI, NA_real_),
    CMI     = if_else(!is.na(TG_HDL) & !is.na(WHtR), TG_HDL * WHtR, NA_real_),
    LAP     = case_when(
      Sex == "Male"   & !is.na(WAIST) & !is.na(TG) ~ (WAIST - 65) * TG,
      Sex == "Female" & !is.na(WAIST) & !is.na(TG) ~ (WAIST - 58) * TG,
      TRUE ~ NA_real_
    ),
    METS_IR = if_else(!is.na(GLUCOSE) & !is.na(TG) & !is.na(HDL) & HDL > 1 &
                        !is.na(BMI),
                      log(2 * GLUCOSE + TG) * BMI / log(HDL), NA_real_)
  )

# == Row specification =========================================================

table_spec <- tribble(
  ~Section, ~Label, ~Variable, ~Type,
  "SOCIODEMOGRAPHIC", "Age, years",                    "Age",       "cont",
  "SOCIODEMOGRAPHIC", "Age <40 years",                 "AGE_LT40",  "cat",
  "SOCIODEMOGRAPHIC", "Age 40\u201359 years",          "AGE_40_59", "cat",
  "SOCIODEMOGRAPHIC", "Age \u226560 years",            "AGE_GE60",  "cat",
  "SOCIODEMOGRAPHIC", "Male sex",                      "SEX_MALE",  "cat",

  "ANTHROPOMETRIC MEASURES", "Height, cm",               "HEIGHT_CM", "cont",
  "ANTHROPOMETRIC MEASURES", "Weight, kg",               "WEIGHT",    "cont",
  "ANTHROPOMETRIC MEASURES", "Waist circumference, cm",  "WAIST",     "cont",
  "ANTHROPOMETRIC MEASURES", "Skeletal muscle mass, kg", "SMM",       "cont",
  "ANTHROPOMETRIC MEASURES", "Body fat mass, kg",        "BFM",       "cont",
  "ANTHROPOMETRIC MEASURES", "Body fat percentage, %",   "PBF",       "cont",

  "BIOCHEMICAL / CLINICAL PARAMETERS", "Fasting glucose, mg/dL",         "GLUCOSE",    "cont",
  "BIOCHEMICAL / CLINICAL PARAMETERS", "Triglycerides, mg/dL",           "TG",         "cont",
  "BIOCHEMICAL / CLINICAL PARAMETERS", "Total cholesterol, mg/dL",       "TC",         "cont",
  "BIOCHEMICAL / CLINICAL PARAMETERS", "HDL cholesterol, mg/dL",         "HDL",        "cont",
  "BIOCHEMICAL / CLINICAL PARAMETERS", "LDL cholesterol, mg/dL",         "LDL",        "cont",
  "BIOCHEMICAL / CLINICAL PARAMETERS", "Systolic blood pressure, mm Hg", "SBP",        "cont",
  "BIOCHEMICAL / CLINICAL PARAMETERS", "Diastolic blood pressure, mm Hg","DBP",        "cont",
  "BIOCHEMICAL / CLINICAL PARAMETERS", "Heart rate, bpm",                "HEART_RATE", "cont",

  "BIOCHEMICAL INDICES", "TG/HDL ratio",                      "TG_HDL",  "cont",
  "BIOCHEMICAL INDICES", "TC/HDL ratio",                      "TC_HDL",  "cont",
  "BIOCHEMICAL INDICES", "LDL/HDL ratio",                     "LDL_HDL", "cont",
  "BIOCHEMICAL INDICES", "HDL/LDL ratio",                     "HDL_LDL", "cont",
  "BIOCHEMICAL INDICES", "Non-HDL cholesterol",               "NON_HDL", "cont",
  "BIOCHEMICAL INDICES", "Atherogenic index of plasma (AIP)", "AIP",     "cont",
  "BIOCHEMICAL INDICES", "Triglyceride-glucose index (TyG)",  "TyG",     "cont",
  "BIOCHEMICAL INDICES", "Lipid comprehensive index (LCI)",   "LCI",     "cont",

  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "BMI, kg/m\u00b2",                         "BMI",           "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Normal weight",                           "NORMAL_WEIGHT", "cat",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Overweight",                              "OVERWEIGHT",    "cat",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Obesity",                                 "OBESITY",       "cat",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Waist-to-height ratio (WHtR)",            "WHtR",          "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Body roundness index (BRI)",              "BRI",           "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Relative fat mass (RFM)",                 "RFM",           "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Fat mass index (FMI), kg/m\u00b2",        "FMI",           "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Skeletal muscle index (SMI), kg/m\u00b2", "SMI",           "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Fat-to-muscle ratio",                     "FAT_MUSCLE",    "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Weight-adjusted waist index (WWI)",       "WWI",           "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "A body shape index (ABSI)",               "ABSI",          "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Conicity index",                          "C_INDEX",       "cont",
  "ANTHROPOMETRIC / BODY COMPOSITION INDICES", "Fat-free mass index (FFMI), kg/m\u00b2",  "FFMI",          "cont",

  "COMPOSITE INDICES", "TyG\u2013BMI",                 "TyG_BMI",  "cont",
  "COMPOSITE INDICES", "TyG\u2013waist circumference", "TyG_WC",   "cont",
  "COMPOSITE INDICES", "TyG\u2013WHtR",                "TyG_WHtR", "cont",
  "COMPOSITE INDICES", "CMI index",                      "CMI",      "cont",
  "COMPOSITE INDICES", "LAP index",                      "LAP",      "cont",
  "COMPOSITE INDICES", "VAI index",                      "VAI",      "cont",
  "COMPOSITE INDICES", "MCMI index",                     "MCMI",     "cont",
  "COMPOSITE INDICES", "METS-IR index",                  "METS_IR",  "cont",

  "OUTCOMES", "Cardiovascular disease", "CVD", "cat",
  "OUTCOMES", "Fatty liver disease",    "FLD", "cat",
  "OUTCOMES", "Chronic kidney disease", "CKD", "cat"
)

# == Groups ====================================================================

idx_mh_a  <- dat$defA_status == "MH"
idx_muh_a <- dat$defA_status == "MUH"
idx_mh_b  <- dat$defB_status == "MH"
idx_muh_b <- dat$defB_status == "MUH"

summarise_row <- function(var, label, type) {
  x <- dat[[var]]

  if (type == "cont") {
    fmt <- fmt_mean_sd
    p1  <- p_continuous(x, dat$defA_status)
    p2  <- p_continuous(x, dat$defB_status)
    p3  <- p_continuous_groups(x[idx_mh_a],  x[idx_mh_b])
    p4  <- p_continuous_groups(x[idx_muh_a], x[idx_muh_b])
  } else {
    fmt <- fmt_n_pct
    p1  <- p_categorical(x, dat$defA_status)
    p2  <- p_categorical(x, dat$defB_status)
    p3  <- p_categorical_groups(x[idx_mh_a],  x[idx_mh_b])
    p4  <- p_categorical_groups(x[idx_muh_a], x[idx_muh_b])
  }

  tibble(
    Variable = label,
    Overall  = fmt(x),
    MH_A     = fmt(x[idx_mh_a]),
    MUH_A    = fmt(x[idx_muh_a]),
    p1       = fmt_p(p1),
    MH_B     = fmt(x[idx_mh_b]),
    MUH_B    = fmt(x[idx_muh_b]),
    p2       = fmt_p(p2),
    p3       = fmt_p(p3),
    p4       = fmt_p(p4)
  )
}

section_row <- function(section) {
  tibble(
    Variable = section, Overall = "", MH_A = "", MUH_A = "", p1 = "",
    MH_B = "", MUH_B = "", p2 = "", p3 = "", p4 = ""
  )
}

table_1 <- bind_rows(
  lapply(unique(table_spec$Section), function(sec) {
    sub <- table_spec %>% filter(Section == sec)
    bind_rows(
      section_row(sec),
      bind_rows(lapply(seq_len(nrow(sub)), function(i) {
        summarise_row(sub$Variable[i], sub$Label[i], sub$Type[i])
      }))
    )
  })
)

# == Column headers with sample sizes ==========================================

fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

names(table_1) <- c(
  "Variable",
  paste0("Overall\n(N=", fmt_n(nrow(dat)), ")"),
  paste0("MH Definition A\n(n=", fmt_n(sum(idx_mh_a)), ")"),
  paste0("MUH Definition A\n(n=", fmt_n(sum(idx_muh_a)), ")"),
  "p*",
  paste0("MH Definition B\n(n=", fmt_n(sum(idx_mh_b)), ")"),
  paste0("MUH Definition B\n(n=", fmt_n(sum(idx_muh_b)), ")"),
  "p**",
  "p***",
  "p****"
)

# == Output ====================================================================

print(table_1, n = Inf, width = Inf)

table_1_footnote <- paste0(
  "Data are mean (SD) or n (%). Metabolic health status was classified as ",
  "metabolically healthy (MH) or metabolically unhealthy (MUH) under two ",
  "operational definitions: Definition A, based on the harmonised metabolic ",
  "syndrome criteria (MUH: \u22652 of 5 components), and Definition B, the ",
  "empirically derived definition of Zembic and colleagues (MUH: \u22651 of 3 ",
  "components). p* compares MH versus MUH participants under Definition A. ",
  "p** compares MH versus MUH participants under Definition B. p*** compares ",
  "participants classified as MH under Definition A with those classified as ",
  "MH under Definition B. p**** compares participants classified as MUH under ",
  "Definition A with those classified as MUH under Definition B. Because both ",
  "definitions were applied to the same participants, cross-definition ",
  "comparisons involve overlapping groups and are descriptive. p values were ",
  "obtained with the Kruskal\u2013Wallis test for continuous variables and the ",
  "\u03c7\u00b2 test for categorical variables. Cardiometabolic indices are grouped ",
  "as biochemical, anthropometric and body-composition, and composite indices ",
  "(formulas in table 2). Fasting glucose was analysed only in ",
  "participants reporting \u22658 h of fasting. AIP=atherogenic index of plasma. ",
  "TyG=triglyceride\u2013glucose index. LCI=lipid comprehensive index. ",
  "BMI=body-mass index. WHtR=waist-to-height ratio. BRI=body roundness index. ",
  "RFM=relative fat mass. FMI=fat mass index. SMI=skeletal muscle index. ",
  "WWI=weight-adjusted waist index. ABSI=a body shape index. FFMI=fat-free mass ",
  "index. CMI=cardiometabolic index. LAP=lipid accumulation product. ",
  "VAI=visceral adiposity index. MCMI=modified cardiometabolic index. ",
  "METS-IR=metabolic score for insulin resistance."
)

cat("\n", table_1_footnote, "\n")
