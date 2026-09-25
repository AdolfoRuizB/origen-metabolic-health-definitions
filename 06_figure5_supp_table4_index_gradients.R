###############################################################################
# Figure 5 and Supplementary Table 4.
#   Figure 5: adiposity-related gradients in the predicted probability of
#   fatty liver disease (FLD) within metabolic phenotypes, under Definition A
#   and Definition B.
#   Supplementary Table 4: index-by-phenotype interaction and MH versus MUH
#   contrasts for all index-outcome pairs under both definitions.
# Both outputs come from the same models, so they are produced by this
# single script.
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Analysis (for each of the 26 non-BMI cardiometabolic indices and each
# outcome, separately under each definition):
#   - Logistic regression: outcome ~ z(index) * phenotype + age + sex,
#     where z(index) is the index standardised to mean 0 and SD 1 and
#     phenotype is one of six categories (MH/MUH x normal weight /
#     overweight / obesity). Phenotypes with fewer than 30 participants or
#     fewer than five events in a given model are excluded from that model.
#   - Phenotype-specific ORs per 1 SD increase in the index (emmeans).
#   - Heterogeneity across phenotypes: global likelihood-ratio test of the
#     index-by-phenotype interaction.
#   - Differences in per-SD associations between MH and MUH participants,
#     expressed as ratios of ORs, within each BMI category and overall; the
#     overall contrast is estimated from a model including the
#     index-by-metabolic-status interaction, age, and sex.
#   - Holm procedure applied to the overall MH versus MUH contrasts across
#     all index-outcome pairs (26 x 3) within each definition.
#   - Predicted probabilities between the 5th and 95th percentiles of the
#     index, at the mean age, for women.
#
# Figure 5 shows FLD for the three indices with the lowest Holm-adjusted
# p values for the MH versus MUH contrast under Definition B (FMI, WHtR,
# and BRI). Rows = indices; columns = Definition A (left) and Definition B
# (right).
#
# Output (printed to the graphics device and console; no files are written):
#   - Figure 5
#   - Console summaries: FLD indices with the lowest Holm-adjusted p values
#     under Definition B; per-SD ORs by phenotype, interaction p values,
#     Holm-adjusted p values, and overall ratios of ORs for the indices
#     shown in Figure 5; and Supplementary Table 4 (global interaction p,
#     overall ratio of ORs, and Holm-adjusted p for all index-outcome pairs
#     under both definitions)
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs; MUH = 2 or more of 5 components
# Definition B: empirically derived definition (Zembic et al.);
#               MUH = 1 or more of 3 components
#
# Analytical population: BMI >= 18.5 kg/m2, sex recorded as female or male,
# and metabolic health status resolvable under both definitions (same as
# Table 1). Each model uses participants with complete data for the index,
# outcome, phenotype, age, and sex.
#
# R version 4.6.0
# Packages: dplyr, readr, stringr, tidyr, ggplot2, ggrepel, patchwork, emmeans
###############################################################################

# -- Packages -----------------------------------------------------------------
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(emmeans)

# -- Settings -----------------------------------------------------------------
# Path to the oriGen questionnaire data file (tab-separated).
# Column names correspond to the original oriGen data release.
data_path <- "path/to/oriGen_questionnaire_data.tsv"

# Reference year used to estimate age from year of birth
REFERENCE_YEAR <- 2025L

raw_data <- read_tsv(data_path, show_col_types = FALSE)

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
  dplyr::case_when(
    y %in% c("TRUE", "VERDADERO", "SI", "1", "YES") ~ 1L,
    y %in% c("FALSE", "FALSO", "NO", "0")           ~ 0L,
    TRUE                                            ~ NA_integer_
  )
}

# Values outside the prespecified physiological range [lo, hi] are set to NA.
clip_range <- function(x, lo, hi) dplyr::if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)

###############################################################################
# 1. DATA PREPARATION: indices, definitions, and phenotypes
###############################################################################

dat <- raw_data %>%
  transmute(
    # Age sources
    AGE_REGISTRATION_RAW = suppressWarnings(as.numeric(PARTICIPANTE.EDAD_REGISTRO)),
    BIRTH_YEAR           = suppressWarnings(as.numeric(PARTICIPANTE.YYYY_NACIMIENTO)),
    AGE_INBODY_RAW       = suppressWarnings(as.numeric(INBODY.Age)),

    SEX_RAW = str_to_upper(trimws(as.character(PARTICIPANTE.SEXO))),
    Sex = dplyr::case_when(
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
    WHR        = suppressWarnings(as.numeric(INBODY.WHR_Waist_Hip_Ratio)),
    TG         = suppressWarnings(as.numeric(METABOLICOS.TRIGLICERIDOS)),
    HDL        = suppressWarnings(as.numeric(METABOLICOS.HDL)),
    LDL        = suppressWarnings(as.numeric(METABOLICOS.LDL)),
    TC         = suppressWarnings(as.numeric(METABOLICOS.COLESTEROL)),
    SBP        = suppressWarnings(as.numeric(METABOLICOS.SISTOLICA)),
    DBP        = suppressWarnings(as.numeric(METABOLICOS.DIASTOLICA)),
    FASTING_HOURS = suppressWarnings(as.numeric(ANTROPOMETRIA.HORAS_AYUNO)),
    # Glucose is considered a fasting value only after >= 8 h of fasting.
    GLUCOSE    = if_else(FASTING_HOURS >= 8,
                         suppressWarnings(as.numeric(METABOLICOS.GLUCOSA)), NA_real_),
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
  mutate(
    AGE_BIRTH_YEAR = if_else(
      !is.na(BIRTH_YEAR) &
        BIRTH_YEAR >= 1900 &
        BIRTH_YEAR <= REFERENCE_YEAR,
      as.numeric(REFERENCE_YEAR - BIRTH_YEAR),
      NA_real_
    ),

    AGE_REGISTRATION = clip_range(AGE_REGISTRATION_RAW, 18, 110),
    AGE_BIRTH_YEAR   = clip_range(AGE_BIRTH_YEAR,       18, 110),
    AGE_INBODY       = clip_range(AGE_INBODY_RAW,       18, 110),

    # Age: age at registration; if missing, age estimated from year of
    # birth; if missing, age recorded by the bioimpedance analyser
    Age = coalesce(AGE_REGISTRATION, AGE_BIRTH_YEAR, AGE_INBODY)
  ) %>%

  # Prespecified physiological ranges, applied before calculating the indices
  mutate(
    WEIGHT    = clip_range(WEIGHT,     30,  300),
    HEIGHT_CM = clip_range(HEIGHT_CM, 130,  220),
    HEIGHT_M  = clip_range(HEIGHT_M, 1.30, 2.20),
    BMI       = clip_range(BMI,        12,   70),
    WAIST     = clip_range(WAIST,      50,  200),
    WHR       = clip_range(WHR,       0.6,  1.5),
    GLUCOSE   = clip_range(GLUCOSE,    40,  600),
    TG        = clip_range(TG,         20, 2000),
    TC        = clip_range(TC,         50,  500),
    HDL       = clip_range(HDL,        10,  150),
    LDL       = clip_range(LDL,        10,  400),
    SBP       = clip_range(SBP,        70,  260),
    DBP       = clip_range(DBP,        40,  160),
    SMM       = clip_range(SMM,         5,   70),
    BFM       = clip_range(BFM,         2,  120)
  ) %>%
  mutate(
    # -- Outcomes ----------------------------------------------------------------
    # 1 if any item is positive; 0 if none is positive and at least one item
    # is evaluable; NA only if all items are missing or uninterpretable.
    CVD = dplyr::case_when(
      cvd_mi == 1L | cvd_stroke == 1L | cvd_ischaemic == 1L |
        cvd_pad == 1L | cvd_athero == 1L | cvd_hf == 1L ~ 1L,
      !is.na(cvd_mi) | !is.na(cvd_stroke) | !is.na(cvd_ischaemic) |
        !is.na(cvd_pad) | !is.na(cvd_athero) | !is.na(cvd_hf) ~ 0L,
      TRUE ~ NA_integer_
    ),

    FLD = dplyr::case_when(
      fatty_liver == 1L ~ 1L,
      fatty_liver == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),

    CKD = dplyr::case_when(
      ckd_diagnosis == 1L | ckd_renal_insf == 1L |
        ckd_dialysis == 1L | ckd_haemodial == 1L ~ 1L,
      !is.na(ckd_diagnosis) | !is.na(ckd_renal_insf) |
        !is.na(ckd_dialysis) | !is.na(ckd_haemodial) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # BMI category (NA for BMI < 18.5 kg/m2 or missing BMI)
    BMI_CATEGORY = dplyr::case_when(
      BMI >= 18.5 & BMI < 25 ~ "NW",
      BMI >= 25   & BMI < 30 ~ "OW",
      BMI >= 30              ~ "OB",
      TRUE                   ~ NA_character_
    ),

    # -- 26 non-BMI cardiometabolic indices --------------------------------------

    # Biochemical
    TG_HDL  = TG / HDL,
    TC_HDL  = TC / HDL,
    LDL_HDL = LDL / HDL,
    HDL_LDL = HDL / LDL,
    NON_HDL = TC - HDL,
    AIP     = log10(TG / HDL),
    TyG     = if_else(!is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0,
                      log((TG * GLUCOSE) / 2), NA_real_),
    LCI     = if_else(!is.na(HDL) & HDL > 0 & !is.na(TC) & !is.na(TG) & !is.na(LDL),
                      (TC * TG * LDL) / HDL, NA_real_),

    # Anthropometric and body composition
    WHtR    = if_else(!is.na(WAIST) & !is.na(HEIGHT_CM) & HEIGHT_CM >= 130 & HEIGHT_CM <= 220,
                      WAIST / HEIGHT_CM, NA_real_),
    BRI     = if_else(!is.na(WAIST) & !is.na(HEIGHT_CM) & HEIGHT_CM > 0,
                      364.2 - 365.5 * sqrt(pmax(0, 1 - ((WAIST / (2 * pi))^2 / ((HEIGHT_CM / 2)^2)))),
                      NA_real_),
    RFM     = dplyr::case_when(
      Sex == "Male"   & !is.na(HEIGHT_M) & !is.na(WAIST) & WAIST > 0 ~ 64 - 20 * (HEIGHT_M / (WAIST / 100)),
      Sex == "Female" & !is.na(HEIGHT_M) & !is.na(WAIST) & WAIST > 0 ~ 76 - 20 * (HEIGHT_M / (WAIST / 100)),
      TRUE ~ NA_real_
    ),
    FMI     = if_else(!is.na(BFM) & !is.na(HEIGHT_M) & HEIGHT_M >= 1.30 & HEIGHT_M <= 2.20,
                      BFM / HEIGHT_M^2, NA_real_),
    SMI     = if_else(!is.na(SMM) & !is.na(HEIGHT_M) & HEIGHT_M >= 1.30 & HEIGHT_M <= 2.20,
                      SMM / HEIGHT_M^2, NA_real_),
    FAT_MUSCLE = if_else(!is.na(BFM) & !is.na(SMM) & SMM > 0, BFM / SMM, NA_real_),
    WWI     = if_else(!is.na(WAIST) & !is.na(WEIGHT) & WEIGHT > 0, WAIST / sqrt(WEIGHT), NA_real_),
    ABSI    = if_else(!is.na(BMI) & BMI > 0 & !is.na(HEIGHT_M) & HEIGHT_M > 0 & !is.na(WAIST),
                      WAIST / (BMI^(2/3) * HEIGHT_M^(1/2)), NA_real_),
    C_INDEX = if_else(!is.na(WEIGHT) & WEIGHT > 0 & !is.na(HEIGHT_M) & HEIGHT_M > 0 & !is.na(WAIST),
                      (WAIST / 100) / (0.109 * sqrt(WEIGHT / HEIGHT_M)), NA_real_),
    # FFMI = (weight - body fat mass) / height^2 (same as Table 1 and Figure 4)
    FFMI    = if_else(!is.na(WEIGHT) & !is.na(BFM) & !is.na(HEIGHT_M) & HEIGHT_M > 0 & WEIGHT > BFM,
                      (WEIGHT - BFM) / HEIGHT_M^2, NA_real_),

    # Composite
    TyG_BMI  = if_else(!is.na(TyG) & !is.na(BMI), TyG * BMI, NA_real_),
    TyG_WC   = if_else(!is.na(TyG) & !is.na(WAIST), TyG * WAIST, NA_real_),
    TyG_WHtR = if_else(!is.na(TyG) & !is.na(WHtR), TyG * WHtR, NA_real_),
    CMI      = if_else(!is.na(WHtR) & !is.na(TG) & !is.na(HDL) & HDL > 0, (TG / HDL) * WHtR, NA_real_),
    LAP      = dplyr::case_when(
      Sex == "Male"   & !is.na(WAIST) & WAIST > 65 & !is.na(TG) ~ (WAIST - 65) * TG,
      Sex == "Female" & !is.na(WAIST) & WAIST > 58 & !is.na(TG) ~ (WAIST - 58) * TG,
      TRUE ~ NA_real_
    ),
    VAI      = dplyr::case_when(
      Sex == "Male"   & !is.na(HDL) & HDL > 0 & !is.na(BMI) & !is.na(TG) & !is.na(WAIST) ~
        (WAIST / (39.68 + 1.88 * BMI)) * (TG / 1.03) * (1.31 / HDL),
      Sex == "Female" & !is.na(HDL) & HDL > 0 & !is.na(BMI) & !is.na(TG) & !is.na(WAIST) ~
        (WAIST / (36.58 + 1.89 * BMI)) * (TG / 0.81) * (1.52 / HDL),
      TRUE ~ NA_real_
    ),
    MCMI     = if_else(!is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0 & !is.na(HDL) & HDL > 0 & !is.na(WHtR),
                       log((TG * GLUCOSE) / HDL) * WHtR, NA_real_),
    METS_IR  = if_else(!is.na(GLUCOSE) & !is.na(TG) & TG > 0 & GLUCOSE > 0 & !is.na(HDL) & HDL > 0 & !is.na(BMI) & BMI > 0,
                       (log((2 * GLUCOSE) + TG) * BMI) / log(HDL), NA_real_),

    # -------------------------------------------------------------------------
    # METABOLIC HEALTH DEFINITIONS (same as Table 1)
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
    defB_diabetes = dplyr::case_when(
      diabetes == 1L ~ 1L,
      diabetes == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),

    # B2) Systolic blood pressure >= 130 mm Hg or antihypertensive treatment
    defB_bp = dplyr::case_when(
      (!is.na(SBP) & SBP >= 130) |
        antihypertensive == 1L ~ 1L,
      !is.na(SBP) |
        !is.na(antihypertensive) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # B3) Waist-to-hip ratio >= 0.95 in women or >= 1.03 in men
    defB_whr = dplyr::case_when(
      Sex == "Female" & !is.na(WHR) & WHR >= 0.95 ~ 1L,
      Sex == "Female" & !is.na(WHR) & WHR <  0.95 ~ 0L,
      Sex == "Male"   & !is.na(WHR) & WHR >= 1.03 ~ 1L,
      Sex == "Male"   & !is.na(WHR) & WHR <  1.03 ~ 0L,
      TRUE ~ NA_integer_
    ),

    defB_valid_components = rowSums(!is.na(cbind(
      defB_diabetes,
      defB_bp,
      defB_whr
    ))),

    defB_n = if_else(
      defB_valid_components == 3L,
      defB_diabetes + defB_bp + defB_whr,
      NA_integer_
    ),

    defB_status = dplyr::case_when(
      !is.na(defB_n) & defB_n == 0L ~ "MH",
      !is.na(defB_n) & defB_n >= 1L ~ "MUH",
      TRUE ~ NA_character_
    ),

    # ============================ DEFINITION A ==============================
    # MH = 0-1 components; MUH = 2 or more of 5 components.

    # A1) Waist circumference >= 90 cm in men or >= 80 cm in women
    defA_waist = dplyr::case_when(
      Sex == "Male"   & !is.na(WAIST) & WAIST >= 90 ~ 1L,
      Sex == "Male"   & !is.na(WAIST) & WAIST <  90 ~ 0L,
      Sex == "Female" & !is.na(WAIST) & WAIST >= 80 ~ 1L,
      Sex == "Female" & !is.na(WAIST) & WAIST <  80 ~ 0L,
      TRUE ~ NA_integer_
    ),

    # A2) Triglycerides >= 150 mg/dL or lipid-lowering treatment
    defA_tg = dplyr::case_when(
      (!is.na(TG) & TG >= 150) |
        lipid_lowering == 1L ~ 1L,
      !is.na(TG) |
        !is.na(lipid_lowering) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # A3) HDL cholesterol < 40 mg/dL in men or < 50 mg/dL in women,
    #     or lipid-lowering treatment
    defA_hdl = dplyr::case_when(
      (Sex == "Male"   & !is.na(HDL) & HDL < 40) |
        (Sex == "Female" & !is.na(HDL) & HDL < 50) |
        lipid_lowering == 1L ~ 1L,
      !is.na(HDL) |
        !is.na(lipid_lowering) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # A4) Systolic >= 130 mm Hg, diastolic >= 85 mm Hg, or antihypertensive
    #     treatment
    defA_bp = dplyr::case_when(
      (!is.na(SBP) & SBP >= 130) |
        (!is.na(DBP) & DBP >= 85) |
        antihypertensive == 1L ~ 1L,
      !is.na(SBP) |
        !is.na(DBP) |
        !is.na(antihypertensive) ~ 0L,
      TRUE ~ NA_integer_
    ),

    # A5) Fasting glucose >= 100 mg/dL, self-reported diabetes (any type),
    #     or glucose-lowering treatment (oral agents, insulin, or any)
    defA_glucose = dplyr::case_when(
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
      defA_waist,
      defA_tg,
      defA_hdl,
      defA_bp,
      defA_glucose
    ))),

    defA_n = if_else(
      defA_valid_components == 5L,
      defA_waist + defA_tg + defA_hdl + defA_bp + defA_glucose,
      NA_integer_
    ),

    defA_status = dplyr::case_when(
      !is.na(defA_n) & defA_n <= 1L ~ "MH",
      !is.na(defA_n) & defA_n >= 2L ~ "MUH",
      TRUE ~ NA_character_
    ),

    # -- Six phenotypes per definition (metabolic status x BMI category) -------
    pheno6_defB = dplyr::case_when(
      defB_status == "MH"  & BMI_CATEGORY == "NW" ~ "MHNW",
      defB_status == "MUH" & BMI_CATEGORY == "NW" ~ "MUNW",
      defB_status == "MH"  & BMI_CATEGORY == "OW" ~ "MHOW",
      defB_status == "MUH" & BMI_CATEGORY == "OW" ~ "MUHOW",
      defB_status == "MH"  & BMI_CATEGORY == "OB" ~ "MHO",
      defB_status == "MUH" & BMI_CATEGORY == "OB" ~ "MUHO",
      TRUE ~ NA_character_
    ),

    pheno6_defA = dplyr::case_when(
      defA_status == "MH"  & BMI_CATEGORY == "NW" ~ "MHNW",
      defA_status == "MUH" & BMI_CATEGORY == "NW" ~ "MUNW",
      defA_status == "MH"  & BMI_CATEGORY == "OW" ~ "MHOW",
      defA_status == "MUH" & BMI_CATEGORY == "OW" ~ "MUHOW",
      defA_status == "MH"  & BMI_CATEGORY == "OB" ~ "MHO",
      defA_status == "MUH" & BMI_CATEGORY == "OB" ~ "MUHO",
      TRUE ~ NA_character_
    ),

    # Women are the reference category (predicted probabilities are for women)
    Sex = factor(Sex, levels = c("Female", "Male"))
  ) %>%

  # Analytical population (same as Table 1); a non-missing BMI category
  # implies valid BMI >= 18.5 kg/m2
  filter(
    !is.na(BMI_CATEGORY),
    !is.na(Sex),
    !is.na(defA_status),
    !is.na(defB_status)
  )

cat("Analytical population:", nrow(dat), "\n")

###############################################################################
# 2. SPECIFICATIONS
###############################################################################

# The 26 non-BMI indices (variable -> label)
index_def <- tibble::tribble(
  ~var,          ~label,
  "TG_HDL",      "TG/HDL ratio",
  "TC_HDL",      "TC/HDL ratio",
  "LDL_HDL",     "LDL/HDL ratio",
  "HDL_LDL",     "HDL/LDL ratio",
  "NON_HDL",     "Non-HDL cholesterol",
  "AIP",         "AIP",
  "TyG",         "TyG index",
  "LCI",         "LCI",
  "WHtR",        "WHtR",
  "BRI",         "BRI",
  "RFM",         "RFM",
  "FMI",         "Fat mass index (FMI)",
  "SMI",         "SMI",
  "FAT_MUSCLE",  "Fat-to-muscle ratio",
  "WWI",         "WWI",
  "ABSI",        "ABSI",
  "C_INDEX",     "Conicity index",
  "FFMI",        "FFMI",
  "TyG_BMI",     "TyG\u2013BMI",
  "TyG_WC",      "TyG\u2013waist circumference",
  "TyG_WHtR",    "TyG\u2013WHtR",
  "CMI",         "CMI",
  "LAP",         "LAP",
  "VAI",         "VAI",
  "MCMI",        "MCMI",
  "METS_IR",     "METS-IR"
)

outcome_def <- tibble::tribble(
  ~var,  ~label,
  "CVD", "Cardiovascular disease (CVD)",
  "FLD", "Fatty liver disease (FLD)",
  "CKD", "Chronic kidney disease (CKD)"
)

# Definitions, in the order used for the tables. In Figure 5 the columns
# are shown as Definition A (left) and Definition B (right).
definition_def <- tibble::tribble(
  ~key,   ~pheno_var,     ~status_var,    ~label,
  "defB", "pheno6_defB",  "defB_status",  "Definition B",
  "defA", "pheno6_defA",  "defA_status",  "Definition A"
)

# Fixed colours for the six phenotypes
phenotype_colours <- c(
  "MHNW"  = "#2166AC", "MUNW"  = "#92C5DE",
  "MHOW"  = "#4DAC26", "MUHOW" = "#B8E186",
  "MHO"   = "#D6604D", "MUHO"  = "#F4A582"
)
phenotype_order <- c("MHNW", "MUNW", "MHOW", "MUHOW", "MHO", "MUHO")

# MH vs MUH pairs within each BMI category
bmi_matched_pairs <- list(c("MHNW", "MUNW"), c("MHOW", "MUHOW"), c("MHO", "MUHO"))

###############################################################################
# 3. CORE: index-by-phenotype interaction model for one index-outcome pair
#    Returns predicted curves, per-SD ORs by phenotype, global interaction p,
#    and BMI-matched MH vs MUH slope contrasts
###############################################################################

analyse_pair <- function(df, idx_var, idx_label, out_var, pheno_var) {

  d <- df %>%
    transmute(
      out   = .data[[out_var]],
      idx   = .data[[idx_var]],
      pheno = .data[[pheno_var]],
      Age, Sex
    ) %>%
    filter(is.finite(idx), !is.na(out), !is.na(pheno), !is.na(Age), !is.na(Sex))

  # Keep phenotypes with >= 30 participants and >= 5 events
  kept <- d %>% group_by(pheno) %>%
    summarise(n = n(), ev = sum(out == 1), .groups = "drop") %>%
    filter(n >= 30, ev >= 5) %>% pull(pheno)
  d <- d %>% filter(pheno %in% kept)
  if (n_distinct(d$pheno) < 2) return(NULL)

  d$pheno <- factor(d$pheno, levels = phenotype_order[phenotype_order %in% kept])

  # Standardise the index -> OR per 1 SD
  mu <- mean(d$idx); sdv <- sd(d$idx)
  if (!is.finite(sdv) || sdv == 0) return(NULL)
  d$z_idx <- (d$idx - mu) / sdv

  f_full <- out ~ z_idx * pheno + Age + Sex
  f_red  <- out ~ z_idx + pheno + Age + Sex

  m_full <- tryCatch(glm(f_full, data = d, family = binomial()), error = function(e) NULL)
  m_red  <- tryCatch(glm(f_red,  data = d, family = binomial()), error = function(e) NULL)
  if (is.null(m_full) || is.null(m_red)) return(NULL)

  # Global index-by-phenotype interaction (likelihood-ratio test)
  p_int <- tryCatch(anova(m_red, m_full, test = "LRT")[["Pr(>Chi)"]][2],
                    error = function(e) NA_real_)

  # Slope of the index (log-odds per SD) within each phenotype
  emt <- tryCatch(emmeans::emtrends(m_full, ~ pheno, var = "z_idx"),
                  error = function(e) NULL)
  if (is.null(emt)) return(NULL)
  emt_s <- as.data.frame(summary(emt, infer = c(TRUE, TRUE)))
  trend_col <- grep("\\.trend$", names(emt_s), value = TRUE)[1]
  slopes <- tibble(
    index = idx_label, outcome = out_var, pheno = as.character(emt_s$pheno),
    slope  = emt_s[[trend_col]],
    OR     = exp(emt_s[[trend_col]]),
    OR_low = exp(emt_s$asymp.LCL), OR_high = exp(emt_s$asymp.UCL),
    p_int_global = p_int
  ) %>%
    # Instability flag: non-finite CI, CI spanning > 3 orders of magnitude,
    # or per-SD OR outside a plausible range (separation or collinearity)
    mutate(unstable = !is.finite(OR_low) | !is.finite(OR_high) |
             (OR_high / OR_low) > 1000 | OR > 10 | OR < 0.1)
  if (any(slopes$unstable))
    cat("   [warning] possible separation:", idx_label, "/", out_var, "->",
        paste(slopes$pheno[slopes$unstable], collapse = ", "), "\n")

  # BMI-matched MH vs MUH slope contrasts (difference in slopes)
  V   <- as.matrix(vcov(emt))
  est <- emt_s[[trend_col]]; names(est) <- as.character(emt_s$pheno)
  rownames(V) <- colnames(V) <- as.character(emt_s$pheno)
  slope_contrast <- function(a, b) {
    if (!(a %in% names(est)) || !(b %in% names(est)))
      return(tibble(pair = paste(a, "vs", b), diff = NA, p_raw = NA_real_))
    diff <- est[a] - est[b]
    se   <- sqrt(V[a, a] + V[b, b] - 2 * V[a, b])
    z    <- diff / se
    tibble(pair = paste(a, "vs", b),
           diff = as.numeric(diff),
           OR_ratio = exp(as.numeric(diff)),     # ratio of per-SD ORs (MH / MUH)
           p_raw = 2 * pnorm(-abs(z)))
  }
  pairs_df <- bind_rows(lapply(bmi_matched_pairs, function(p) slope_contrast(p[1], p[2]))) %>%
    mutate(index = idx_label, outcome = out_var, type = "BMI-matched")

  # -- Predicted probability curves (women, mean age, 5th-95th percentile) ---
  age_mean <- mean(d$Age)
  sex_ref  <- levels(d$Sex)[1]
  z_grid <- seq(quantile(d$z_idx, 0.05), quantile(d$z_idx, 0.95), length.out = 80)
  grid <- expand.grid(z_idx = z_grid, pheno = levels(d$pheno),
                      Age = age_mean, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  grid$pheno <- factor(grid$pheno, levels = levels(d$pheno))
  grid$Sex   <- factor(sex_ref, levels = levels(d$Sex))
  pr <- predict(m_full, newdata = grid, type = "link", se.fit = TRUE)
  curves <- grid %>%
    mutate(
      idx_raw = z_idx * sdv + mu,
      prob    = plogis(pr$fit),
      lo      = plogis(pr$fit - 1.96 * pr$se.fit),
      hi      = plogis(pr$fit + 1.96 * pr$se.fit),
      index   = idx_label, outcome = out_var
    )

  list(curves = curves, slopes = slopes, pairs = pairs_df, p_int = p_int)
}

# Overall MH vs MUH contrast (two-level model): ratio of per-SD ORs and p value
analyse_overall_contrast <- function(df, idx_var, idx_label, out_var, status_var) {
  d <- df %>%
    transmute(out = .data[[out_var]], idx = .data[[idx_var]],
              mh = .data[[status_var]], Age, Sex) %>%
    filter(is.finite(idx), !is.na(out), !is.na(mh), !is.na(Age), !is.na(Sex))
  d$mh <- factor(d$mh, levels = c("MH", "MUH"))
  if (n_distinct(d$mh) < 2) return(NULL)
  sdv <- sd(d$idx); if (!is.finite(sdv) || sdv == 0) return(NULL)
  d$z_idx <- (d$idx - mean(d$idx)) / sdv
  m <- tryCatch(glm(out ~ z_idx * mh + Age + Sex,
                    data = d, family = binomial()), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  emt <- tryCatch(emmeans::emtrends(m, ~ mh, var = "z_idx"), error = function(e) NULL)
  if (is.null(emt)) return(NULL)
  es <- as.data.frame(summary(emt, infer = TRUE))
  tcol <- grep("\\.trend$", names(es), value = TRUE)[1]
  V <- as.matrix(vcov(emt)); est <- es[[tcol]]
  diff <- est[1] - est[2]
  se   <- sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
  tibble(index = idx_label, outcome = out_var, pair = "MH overall vs MUH overall",
         diff = as.numeric(diff), OR_ratio = exp(as.numeric(diff)),
         p_raw = 2 * pnorm(-abs(diff / se)), type = "overall")
}

###############################################################################
# 4. ALL INDEX-OUTCOME PAIRS, WITH HOLM CORRECTION WITHIN EACH DEFINITION
###############################################################################

holm_table <- function(pheno_var, status_var, def_label) {
  cat("\n==== All index-outcome pairs |", def_label, "|", nrow(index_def), "indices x",
      nrow(outcome_def), "outcomes (this may take a while) ====\n")
  rows <- list()
  for (oi in seq_len(nrow(outcome_def))) {
    ov <- outcome_def$var[oi]
    for (ii in seq_len(nrow(index_def))) {
      iv <- index_def$var[ii]; ilb <- index_def$label[ii]
      pr  <- tryCatch(analyse_pair(dat, iv, ilb, ov, pheno_var), error = function(e) NULL)
      ovl <- tryCatch(analyse_overall_contrast(dat, iv, ilb, ov, status_var), error = function(e) NULL)
      rows[[length(rows) + 1L]] <- tibble(
        outcome = ov, outcome_lab = outcome_def$label[oi], var = iv, index = ilb,
        p_int = if (is.null(pr)) NA_real_ else pr$p_int,
        ROR   = if (is.null(ovl)) NA_real_ else ovl$OR_ratio,
        p_raw = if (is.null(ovl)) NA_real_ else ovl$p_raw,
        unstable = if (is.null(pr)) TRUE else any(pr$slopes$unstable)
      )
    }
  }
  bind_rows(rows) %>% mutate(p_adj = p.adjust(p_raw, method = "holm"))
}

tab_defB <- holm_table("pheno6_defB", "defB_status", "Definition B")
tab_defA <- holm_table("pheno6_defA", "defA_status", "Definition A")

# Lookup of Holm-adjusted p by (definition, outcome, index)
padj_lookup <- bind_rows(
  tab_defB %>% transmute(key = "defB", outcome, var, p_adj),
  tab_defA %>% transmute(key = "defA", outcome, var, p_adj)
)
get_padj <- function(k, ov, v) {
  x <- padj_lookup$p_adj[padj_lookup$key == k & padj_lookup$outcome == ov & padj_lookup$var == v]
  if (length(x) == 0) NA_real_ else x[1]
}

###############################################################################
# 5. INDICES SHOWN IN FIGURE 5
###############################################################################

# FLD indices with the lowest Holm-adjusted p values for the overall MH vs
# MUH contrast under Definition B (verification of the selection below)
cat("\n=== FLD: indices with the lowest Holm-adjusted p values (Definition B) ===\n")
print(as.data.frame(
  tab_defB %>%
    filter(outcome == "FLD", !is.na(p_adj)) %>%
    arrange(p_adj) %>%
    slice_head(n = 3) %>%
    select(Index = index, ROR, p_raw, p_adj)
), row.names = FALSE)

# Indices shown in Figure 5 (rows, top to bottom)
FIGURE_INDICES <- c("FMI", "WHtR", "BRI")

label_of <- function(v) index_def$label[index_def$var == v]

selected <- tibble(out_var = "FLD", idx_var = FIGURE_INDICES) %>%
  mutate(out_lab   = "Fatty liver disease (FLD)",
         idx_label = vapply(idx_var, label_of, character(1)))

###############################################################################
# 6. FIGURE 5
###############################################################################

# Short index labels for the figure only (tables keep the full label)
short_label_map <- c("Fat mass index (FMI)" = "FMI")
short_label <- function(lab) {
  if (!is.null(lab) && lab %in% names(short_label_map)) unname(short_label_map[lab]) else lab
}

make_panel <- function(idx_var, idx_label, out_var, def_row) {
  r <- analyse_pair(dat, idx_var, idx_label, out_var, def_row$pheno_var)
  if (is.null(r)) return(
    ggplot() + theme_void() +
      annotate("text", x = 0, y = 0, size = 3,
               label = paste0(idx_label, "\n(", out_var, " / ", def_row$label,
                              ")\ninsufficient data")))
  padj <- get_padj(def_row$key, out_var, idx_var)
  curve_labels <- r$curves %>% group_by(pheno) %>% slice_max(idx_raw, n = 1) %>% ungroup()
  p_txt <- ifelse(is.na(r$p_int), "p-int = NA",
                  ifelse(r$p_int < 0.001, "p-int < 0.001", sprintf("p-int = %.3f", r$p_int)))
  holm_txt <- if (!is.na(padj)) paste0(" | Holm p_adj=",
                                       ifelse(padj < 0.001, formatC(padj, format = "e", digits = 1), sprintf("%.3f", padj))) else ""

  ggplot(r$curves, aes(idx_raw, prob, color = pheno, fill = pheno)) +
    geom_line(linewidth = 0.9) +
    ggrepel::geom_text_repel(
      data = curve_labels, aes(label = pheno), size = 2.3, fontface = "bold",
      hjust = 0, direction = "y", nudge_x = diff(range(r$curves$idx_raw)) * 0.04,
      segment.size = 0.2, segment.alpha = 0.4, box.padding = 0.1,
      min.segment.length = 0, show.legend = FALSE) +
    # Interaction p value in the lower-right corner of the panel
    annotate("text", x = Inf, y = -Inf, label = paste0(p_txt, holm_txt),
             hjust = 1.05, vjust = -2.0, size = 2.4, color = "gray30",
             inherit.aes = FALSE) +
    scale_color_manual(values = phenotype_colours) + scale_fill_manual(values = phenotype_colours) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    scale_y_continuous(labels = function(x) round(x * 100)) +
    labs(title = NULL, subtitle = NULL,
         x = short_label(idx_label), y = "Predicted probability (%)") +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
}

# One panel per (index, definition), filled by row in the order of
# definition_def (Definition B, Definition A)
panels <- list(); k <- 0
for (i in seq_len(nrow(selected))) {
  for (j in seq_len(nrow(definition_def))) {
    k <- k + 1
    panels[[k]] <- make_panel(selected$idx_var[i], selected$idx_label[i],
                              selected$out_var[i], definition_def[j, ])
  }
}

# Column order reversed for the figure: Definition A (left), Definition B (right)
n_def   <- nrow(definition_def)
n_rows  <- nrow(selected)
panel_order <- unlist(lapply(seq_len(n_rows), function(i) (i - 1) * n_def + rev(seq_len(n_def))))
figure_panels <- panels[panel_order]

# Column labels at the bottom: "A" (Definition A) and "B" (Definition B)
column_label <- function(txt) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = txt, fontface = "bold", size = 5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}
column_labels <- lapply(LETTERS[seq_len(n_def)], column_label)

figure_5 <- patchwork::wrap_plots(c(figure_panels, column_labels),
                                  nrow = n_rows + 1, ncol = n_def, byrow = TRUE,
                                  heights = c(rep(1, n_rows), 0.12))
print(figure_5)

###############################################################################
# 7. CONSOLE SUMMARY: indices shown in Figure 5
#    Per-SD OR by phenotype, interaction p, Holm-adjusted p, and overall
#    ratio of ORs (MH vs MUH)
###############################################################################
fig5_rows <- list()
for (i in seq_len(nrow(selected))) {
  for (j in seq_len(nrow(definition_def))) {
    dr <- definition_def[j, ]
    r  <- tryCatch(analyse_pair(dat, selected$idx_var[i], selected$idx_label[i],
                                selected$out_var[i], dr$pheno_var), error = function(e) NULL)
    ov <- tryCatch(analyse_overall_contrast(dat, selected$idx_var[i], selected$idx_label[i],
                                            selected$out_var[i], dr$status_var), error = function(e) NULL)
    if (is.null(r)) next
    sl <- r$slopes %>%
      transmute(pheno, val = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high)) %>%
      tidyr::pivot_wider(names_from = pheno, values_from = val)
    pa <- get_padj(dr$key, selected$out_var[i], selected$idx_var[i])
    fig5_rows[[length(fig5_rows) + 1L]] <- bind_cols(
      tibble(
        Definition = dr$label, Outcome = selected$out_lab[i], Index = selected$idx_label[i],
        `p-int`      = ifelse(is.na(r$p_int), "NA",
                              ifelse(r$p_int < 0.001, "<0.001", sprintf("%.3f", r$p_int))),
        `Holm p_adj` = ifelse(is.na(pa), "NA",
                              ifelse(pa < 0.001, formatC(pa, format = "e", digits = 1), sprintf("%.3f", pa))),
        ROR = if (is.null(ov)) NA_real_ else round(ov$OR_ratio, 2)
      ),
      sl)
  }
}
fig5_summary <- bind_rows(fig5_rows)
cat("\n===== Indices shown in Figure 5 =====\n")
cat("Phenotype columns: OR (95% CI) per 1 SD of the index. ROR = ratio of per-SD ORs, MH vs MUH overall.\n")
print(as.data.frame(fig5_summary), row.names = FALSE)

###############################################################################
# 8. SUPPLEMENTARY TABLE 4: all index-outcome pairs under both definitions
###############################################################################
fmt_p   <- function(p) ifelse(is.na(p), "NA",
                              ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
fmt_ror <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
mark    <- function(p) ifelse(!is.na(p) & p < 0.05, "*", "")

sum_B <- tab_defB %>% transmute(outcome, outcome_lab, var, index,
                                `Definition B p-int` = fmt_p(p_int), `Definition B ROR` = fmt_ror(ROR),
                                `Definition B p_adj` = paste0(fmt_p(p_adj), mark(p_adj)), unstable_B = unstable)
sum_A <- tab_defA %>% transmute(outcome, var,
                                `Definition A p-int` = fmt_p(p_int), `Definition A ROR` = fmt_ror(ROR),
                                `Definition A p_adj` = paste0(fmt_p(p_adj), mark(p_adj)), unstable_A = unstable)

all_pairs_summary <- sum_B %>% left_join(sum_A, by = c("outcome", "var")) %>%
  arrange(match(outcome, outcome_def$var), index) %>%
  mutate(Unstable = if_else(unstable_B | unstable_A, "yes", "")) %>%
  select(Outcome = outcome_lab, Index = index,
         `Definition A p-int`, `Definition A ROR`, `Definition A p_adj`,
         `Definition B p-int`, `Definition B ROR`, `Definition B p_adj`, Unstable)

cat("\n===== Supplementary Table 4: all index-outcome pairs =====\n")
cat("* = Holm-adjusted p < 0.05 (within each definition; family =",
    nrow(index_def) * nrow(outcome_def), "tests).\n\n")
print(as.data.frame(all_pairs_summary), row.names = FALSE)
