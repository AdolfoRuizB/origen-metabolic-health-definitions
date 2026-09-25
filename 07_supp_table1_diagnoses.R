###############################################################################
# Supplementary Table 1. Prevalence of the individual diagnoses included in
# the composite outcomes, overall, by BMI category, and by sex
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Columns: Diagnosis | Overall | Normal weight | Overweight | Obesity |
#          p (BMI category) | Male | Female | p (sex)
#
# Data are n (%); percentages are calculated among participants with
# available data for each diagnosis (1 decimal, or 2 decimals when the
# percentage is < 0.1%). p values: chi-squared test; * = p < 0.05.
#
# Output (printed to the console; no files are written): Supplementary
# Table 1, valid denominators for each diagnosis, and the table footnote.
#
# Analytical population (same as Table 1): BMI >= 18.5 kg/m2, sex recorded
# as female or male, and metabolic health status resolvable under both
# Definition A and Definition B.
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs; MUH = 2 or more of 5 components
# Definition B: empirically derived definition (Zembic et al.);
#               MUH = 1 or more of 3 components
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

# == Data preparation ==========================================================
# Only the variables needed to define the analytical population and the
# outcomes are prepared here.

dat_unfiltered <- raw_data %>%
  transmute(

    SEX_RAW = str_to_upper(trimws(as.character(PARTICIPANTE.SEXO))),
    Sex = case_when(
      SEX_RAW %in% c("H", "HOMBRE", "MASCULINO", "MALE", "MAN")         ~ "Male",
      SEX_RAW %in% c("M", "MUJER", "F", "FEMENINO", "FEMALE", "WOMAN") ~ "Female",
      TRUE ~ NA_character_
    ),

    BMI   = suppressWarnings(as.numeric(ANTROPOMETRIA.IMC)),
    WAIST = suppressWarnings(as.numeric(ANTROPOMETRIA.CIRCUNFERENCIA_CINTURA)),
    WHR   = suppressWarnings(as.numeric(INBODY.WHR_Waist_Hip_Ratio)),

    FASTING_HOURS = suppressWarnings(as.numeric(ANTROPOMETRIA.HORAS_AYUNO)),
    GLUCOSE_RAW   = suppressWarnings(as.numeric(METABOLICOS.GLUCOSA)),
    TG            = suppressWarnings(as.numeric(METABOLICOS.TRIGLICERIDOS)),
    HDL           = suppressWarnings(as.numeric(METABOLICOS.HDL)),
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
  # Prespecified physiological ranges (same as Table 1)
  mutate(
    BMI         = clip_range(BMI,          12,   70),
    WAIST       = clip_range(WAIST,        50,  200),
    WHR         = clip_range(WHR,         0.6,  1.5),
    GLUCOSE_RAW = clip_range(GLUCOSE_RAW,  40,  600),
    TG          = clip_range(TG,           20, 2000),
    HDL         = clip_range(HDL,          10,  150),
    SBP         = clip_range(SBP,          70,  260),
    DBP         = clip_range(DBP,          40,  160)
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

    # -------------------------------------------------------------------------
    # METABOLIC HEALTH DEFINITIONS (used only to define the analytical
    # population)
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

cat("\nAnalytical population:", nrow(dat), "\n")

# == Groups ====================================================================

dat <- dat %>%
  mutate(
    BMI_CATEGORY = case_when(
      BMI < 25             ~ "Normal weight",
      BMI >= 25 & BMI < 30 ~ "Overweight",
      BMI >= 30            ~ "Obesity",
      TRUE                 ~ NA_character_
    ),
    BMI_CATEGORY = factor(BMI_CATEGORY, levels = c("Normal weight", "Overweight", "Obesity"))
  )

idx_nw <- !is.na(dat$BMI_CATEGORY) & dat$BMI_CATEGORY == "Normal weight"
idx_ow <- !is.na(dat$BMI_CATEGORY) & dat$BMI_CATEGORY == "Overweight"
idx_ob <- !is.na(dat$BMI_CATEGORY) & dat$BMI_CATEGORY == "Obesity"
idx_m  <- dat$Sex == "Male"
idx_f  <- dat$Sex == "Female"

# == Formatting functions ======================================================

fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

# n (%) with thousands separator; 1 decimal, or 2 decimals if % < 0.1
fmt_n_pct <- function(x) {
  n   <- sum(x == 1L, na.rm = TRUE)
  tot <- sum(!is.na(x))
  if (tot == 0) return(NA_character_)
  pct <- 100 * n / tot
  dec <- if (pct > 0 & round(pct, 1) < 0.1) 2 else 1
  paste0(fmt_n(n), " (", formatC(pct, format = "f", digits = dec), "%)")
}

p_chi <- function(x, grp) {
  df <- data.frame(x = x, g = grp) %>% filter(!is.na(x), !is.na(g))
  tt <- table(df$x, df$g)
  if (nrow(tt) < 2 || ncol(tt) < 2) return(NA_real_)
  suppressWarnings(tryCatch(chisq.test(tt)$p.value, error = function(e) NA_real_))
}

fmt_p_star <- function(p) {
  if (is.na(p)) return(NA_character_)
  txt <- if (p < 0.001) "<0.001" else formatC(p, format = "f", digits = 3)
  if (p < 0.05) paste0(txt, "*") else txt
}

# == Row specification =========================================================

table_spec <- tribble(
  ~Section,                   ~Label,                                 ~Variable,
  "CARDIOVASCULAR DISEASE",   "Cardiovascular disease (overall)",     "CVD",
  "CARDIOVASCULAR DISEASE",   "Acute myocardial infarction",          "cvd_mi",
  "CARDIOVASCULAR DISEASE",   "Cerebral infarction (stroke)",         "cvd_stroke",
  "CARDIOVASCULAR DISEASE",   "Ischaemic coronary artery disease",    "cvd_ischaemic",
  "CARDIOVASCULAR DISEASE",   "Peripheral arterial disease",          "cvd_pad",
  "CARDIOVASCULAR DISEASE",   "Atherosclerosis",                      "cvd_athero",
  "CARDIOVASCULAR DISEASE",   "Heart failure",                        "cvd_hf",
  "FATTY LIVER DISEASE",      "Fatty liver disease",                  "FLD",
  "CHRONIC KIDNEY DISEASE",   "Chronic kidney disease (overall)",     "CKD",
  "CHRONIC KIDNEY DISEASE",   "Chronic kidney disease (diagnosis)",   "ckd_diagnosis",
  "CHRONIC KIDNEY DISEASE",   "Renal insufficiency",                  "ckd_renal_insf",
  "CHRONIC KIDNEY DISEASE",   "Dialysis due to renal disease",        "ckd_dialysis",
  "CHRONIC KIDNEY DISEASE",   "Haemodialysis",                        "ckd_haemodial"
)

summarise_row <- function(var, label) {
  x <- dat[[var]]
  tibble(
    Diagnosis = label,
    Overall   = fmt_n_pct(x),
    NW        = fmt_n_pct(x[idx_nw]),
    OW        = fmt_n_pct(x[idx_ow]),
    OB        = fmt_n_pct(x[idx_ob]),
    p_bmi     = fmt_p_star(p_chi(x, dat$BMI_CATEGORY)),
    Male      = fmt_n_pct(x[idx_m]),
    Female    = fmt_n_pct(x[idx_f]),
    p_sex     = fmt_p_star(p_chi(x, dat$Sex))
  )
}

section_row <- function(section) {
  tibble(Diagnosis = section, Overall = "", NW = "", OW = "", OB = "",
         p_bmi = "", Male = "", Female = "", p_sex = "")
}

supp_table_1 <- bind_rows(
  lapply(unique(table_spec$Section), function(sec) {
    sub <- table_spec %>% filter(Section == sec)
    bind_rows(
      section_row(sec),
      bind_rows(lapply(seq_len(nrow(sub)), function(i) {
        summarise_row(sub$Variable[i], sub$Label[i])
      }))
    )
  })
)

names(supp_table_1) <- c(
  "Diagnosis",
  paste0("Overall (N=",       fmt_n(nrow(dat)),  ")"),
  paste0("Normal weight (n=", fmt_n(sum(idx_nw)), ")"),
  paste0("Overweight (n=",    fmt_n(sum(idx_ow)), ")"),
  paste0("Obesity (n=",       fmt_n(sum(idx_ob)), ")"),
  "p (BMI category)",
  paste0("Male (n=",          fmt_n(sum(idx_m)),  ")"),
  paste0("Female (n=",        fmt_n(sum(idx_f)),  ")"),
  "p (sex)"
)

# Check: valid denominators for each diagnosis (should be close to N)
cat("\nValid denominators for each diagnosis:\n")
print(sapply(table_spec$Variable, function(v) sum(!is.na(dat[[v]]))))

# == Output ====================================================================

print(supp_table_1, n = Inf, width = Inf)

supp_table_1_footnote <- paste0(
  "Data are n (%). Percentages were calculated among participants with ",
  "available data for each diagnosis. p values were obtained with the ",
  "chi-squared test; * p<0.05."
)
cat("\n", supp_table_1_footnote, "\n")
