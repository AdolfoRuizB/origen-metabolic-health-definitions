###############################################################################
# Figure 2. Relative contribution of the components of Definition A and
# Definition B to prevalent cardiovascular disease (CVD), fatty liver
# disease (FLD), and chronic kidney disease (CKD)
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Analysis: for each definition and outcome, a logistic regression model
# including all components of the definition, adjusted for age and sex.
#   - Adjusted odds ratios (95% CI) for each component
#   - Removal of each component in turn: likelihood-ratio test and paired
#     DeLong test of AUCs (full vs reduced model)
#   - General dominance analysis based on McFadden's pseudo-R2 (age and sex
#     included in every model); contributions expressed as percentages of
#     the total attributable to the components
#   - Complementary measure: share of the partial likelihood-ratio
#     chi-square attributable to each component
#
# Output (printed to the graphics device and console; no files are written):
#   - Figure 2, panels a-b:
#       (a) relative contribution of each component (general dominance, %)
#       (b) adjusted odds ratios (95% CI) of each component
#   - Console summaries: likelihood-ratio chi-square, dominance weights,
#     component-level discrimination (AUC, delta AUC, DeLong), and
#     component-block discrimination by definition and outcome
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs (5 components)
# Definition B: empirically derived definition (Zembic et al.) (3 components)
#
# Analytical population: BMI >= 18.5 kg/m2, sex recorded as female or male,
# and metabolic health status resolvable under both definitions (same as
# Table 1). Each model uses participants with complete data for the outcome,
# components, and covariates.
#
# R version 4.6.0
# Packages: dplyr, tidyr, readr, stringr, ggplot2, pROC, patchwork
###############################################################################

# ======================== 0. PACKAGES ========================================
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(pROC)
library(patchwork)

# ======================== 1. SETTINGS ========================================
# Path to the oriGen questionnaire data file (tab-separated).
# Column names correspond to the original oriGen data release.
data_path <- "path/to/oriGen_questionnaire_data.tsv"

# Reference year used to estimate age from year of birth
REFERENCE_YEAR <- 2025L

raw_data <- read_tsv(
  data_path,
  show_col_types = FALSE,
  locale = locale(encoding = "latin1")
)

# ======================== 2. HELPER FUNCTIONS ================================

# Harmonises free-text responses: upper case, no accents, single spaces.
normalise_text <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- trimws(x); x <- str_to_upper(x)
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
clip_range <- function(x, lo, hi) if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)

# ROC curve with fixed direction (controls = 0, cases = 1).
roc_fixed <- function(y, p) pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)

# ======================== 3. DATA PREPARATION (same as Table 1) ==============
dat <- raw_data %>%
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
    BMI           = suppressWarnings(as.numeric(ANTROPOMETRIA.IMC)),
    WAIST         = suppressWarnings(as.numeric(ANTROPOMETRIA.CIRCUNFERENCIA_CINTURA)),
    WHR           = suppressWarnings(as.numeric(INBODY.WHR_Waist_Hip_Ratio)),
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
    fatty_liver    = safe_bool(HISTORIA_MEDICA.HIGADO_GRASO),
    ckd_diagnosis  = safe_bool(HISTORIA_MEDICA.ENFERMEDAD_RENAL_CRONICA),
    ckd_renal_insf = safe_bool(HISTORIA_MEDICA.INSUFICIENCIA_RENAL),
    ckd_dialysis   = safe_bool(HISTORIA_MEDICA.DIALISIS_POR_RENAL),
    ckd_haemodial  = safe_bool(HISTORIA_MEDICA.HEMODIALISIS)
  ) %>%

  # Prespecified physiological ranges, including the three age sources
  mutate(
    BMI              = clip_range(BMI,          12,   70),
    WAIST            = clip_range(WAIST,        50,  200),
    WHR              = clip_range(WHR,         0.6,  1.5),
    GLUCOSE_RAW      = clip_range(GLUCOSE_RAW,  40,  600),
    TG               = clip_range(TG,           20, 2000),
    HDL              = clip_range(HDL,          10,  150),
    SBP              = clip_range(SBP,          70,  260),
    DBP              = clip_range(DBP,          40,  160),
    AGE_REGISTRATION = clip_range(AGE_REGISTRATION_RAW, 18, 110),
    BIRTH_YEAR       = clip_range(BIRTH_YEAR, 1900, REFERENCE_YEAR),
    AGE_BIRTH_YEAR   = clip_range(REFERENCE_YEAR - BIRTH_YEAR, 18, 110),
    AGE_INBODY       = clip_range(AGE_INBODY_RAW, 18, 110)
  ) %>%
  mutate(
    # Age: age at registration; if missing, age estimated from year of
    # birth; if missing, age recorded by the bioimpedance analyser
    Age = coalesce(AGE_REGISTRATION, AGE_BIRTH_YEAR, AGE_INBODY),

    # Glucose is considered a fasting value only after >= 8 h of fasting.
    GLUCOSE = if_else(FASTING_HOURS >= 8, GLUCOSE_RAW, NA_real_),
    Sex     = factor(Sex, levels = c("Female", "Male")),

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
    defB_diabetes = case_when(
      diabetes == 1L ~ 1L,
      diabetes == 0L ~ 0L,
      TRUE ~ NA_integer_
    ),

    # B2) Systolic blood pressure >= 130 mm Hg or antihypertensive treatment
    defB_bp = case_when(
      (!is.na(SBP) & SBP >= 130) |
        antihypertensive == 1L ~ 1L,

      !is.na(SBP) |
        !is.na(antihypertensive) ~ 0L,

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
      defB_diabetes,
      defB_bp,
      defB_whr
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
      (!is.na(TG) & TG >= 150) |
        lipid_lowering == 1L ~ 1L,

      !is.na(TG) |
        !is.na(lipid_lowering) ~ 0L,

      TRUE ~ NA_integer_
    ),

    # A3) HDL cholesterol < 40 mg/dL in men or < 50 mg/dL in women,
    #     or lipid-lowering treatment
    defA_hdl = case_when(
      (Sex == "Male"   & !is.na(HDL) & HDL < 40) |
        (Sex == "Female" & !is.na(HDL) & HDL < 50) |
        lipid_lowering == 1L ~ 1L,

      !is.na(HDL) |
        !is.na(lipid_lowering) ~ 0L,

      TRUE ~ NA_integer_
    ),

    # A4) Systolic >= 130 mm Hg, diastolic >= 85 mm Hg, or antihypertensive
    #     treatment
    defA_bp = case_when(
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

    defA_muh = case_when(
      !is.na(defA_n) & defA_n <= 1L ~ 0L,
      !is.na(defA_n) & defA_n >= 2L ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%

  # Analytical population (same as Table 1)
  filter(
    !is.na(BMI),
    BMI >= 18.5,
    !is.na(Sex),
    !is.na(defA_muh),
    !is.na(defB_muh)
  )

cat("Analytical population (same as Table 1):", nrow(dat), "\n")
cat("With valid age (entering the models):", sum(!is.na(dat$Age)), "\n")

# ======================== 4. SPECIFICATIONS ==================================
components_by_definition <- list(
  defA = c(
    defA_bp      = "Blood pressure",
    defA_glucose = "Glucose / diabetes",
    defA_waist   = "Waist circumference",
    defA_tg      = "Triglycerides",
    defA_hdl     = "Low HDL"
  ),
  defB = c(
    defB_bp       = "Blood pressure",
    defB_diabetes = "Glucose / diabetes",
    defB_whr      = "Central adiposity (WHR)"
  )
)
definition_label <- c(defA = "Definition A", defB = "Definition B")
outcomes <- c(CVD = "Cardiovascular disease (CVD)",
              FLD = "Fatty liver disease (FLD)",
              CKD = "Chronic kidney disease (CKD)")

# Adjustment covariates (all models)
covariates <- c("Age", "Sex")

# ======================== 5. CORE ANALYSIS: one definition x one outcome ====
all_subsets <- function(x) {
  do.call(c, lapply(0:length(x), function(s) {
    if (s == 0) list(character(0)) else utils::combn(x, s, simplify = FALSE)
  }))
}

analyse_components <- function(d, outcome, comps, covars) {
  y <- d[[outcome]]

  # Full model: all components + covariates
  m_full   <- glm(reformulate(c(comps, covars), outcome), data = d, family = binomial())
  ll_full  <- as.numeric(logLik(m_full))
  roc_full <- roc_fixed(y, predict(m_full, type = "response"))
  auc_full <- as.numeric(pROC::auc(roc_full))
  coefs    <- summary(m_full)$coefficients
  vc       <- vcov(m_full)

  # Null model and covariates-only model
  ll_null  <- as.numeric(logLik(glm(reformulate("1", outcome), data = d, family = binomial())))
  m_cov    <- glm(reformulate(covars, outcome), data = d, family = binomial())
  roc_cov  <- roc_fixed(y, predict(m_cov, type = "response"))
  auc_cov  <- as.numeric(pROC::auc(roc_cov))
  delong_block <- tryCatch(pROC::roc.test(roc_full, roc_cov, method = "delong", paired = TRUE)$p.value,
                           error = function(e) NA_real_)

  # General dominance analysis (McFadden pseudo-R2); covariates in every model
  subset_key <- function(S) if (length(S) == 0) "<none>" else paste(sort(S), collapse = "+")
  r2_map <- new.env(parent = emptyenv())
  for (S in all_subsets(comps)) {
    m <- glm(reformulate(c(covars, S), outcome), data = d, family = binomial())
    assign(subset_key(S), 1 - as.numeric(logLik(m)) / ll_null, envir = r2_map)
  }
  get_r2 <- function(S) get(subset_key(S), envir = r2_map)
  dominance <- vapply(comps, function(i) {
    others <- setdiff(comps, i)
    mean(vapply(0:length(others), function(s) {
      Ss <- if (s == 0) list(character(0)) else utils::combn(others, s, simplify = FALSE)
      mean(vapply(Ss, function(S) get_r2(c(S, i)) - get_r2(S), numeric(1)))
    }, numeric(1)))
  }, numeric(1))

  # Component-level results
  rows <- lapply(comps, function(cmp) {
    beta <- if (cmp %in% rownames(coefs)) coefs[cmp, "Estimate"] else NA_real_
    se   <- if (cmp %in% rownames(coefs)) sqrt(vc[cmp, cmp])     else NA_real_

    # Reduced model without the component: LR test and paired DeLong test
    m_red   <- glm(reformulate(c(setdiff(comps, cmp), covars), outcome), data = d, family = binomial())
    ll_red  <- as.numeric(logLik(m_red))
    lr_chi2 <- max(0, 2 * (ll_full - ll_red))
    p_lr    <- pchisq(lr_chi2, df = 1, lower.tail = FALSE)
    roc_red <- roc_fixed(y, predict(m_red, type = "response"))
    auc_red <- as.numeric(pROC::auc(roc_red))
    delong_p <- tryCatch(pROC::roc.test(roc_full, roc_red, method = "delong", paired = TRUE)$p.value,
                         error = function(e) NA_real_)

    # Component alone (unadjusted)
    m_alone   <- glm(reformulate(cmp, outcome), data = d, family = binomial())
    auc_alone <- tryCatch(as.numeric(pROC::auc(roc_fixed(y, predict(m_alone, type = "response")))),
                          error = function(e) NA_real_)

    tibble(Component_var = cmp, beta = beta, se = se,
           lr_chi2 = lr_chi2, p_LR = p_lr,
           auc_alone = auc_alone, auc_red = auc_red, dAUC = auc_full - auc_red,
           delong_p = delong_p)
  })

  out <- bind_rows(rows) %>%
    mutate(
      OR = exp(beta), OR_lo = exp(beta - 1.96 * se), OR_hi = exp(beta + 1.96 * se),
      McFadden_inc = unname(dominance[Component_var]),
      share_lr = 100 * lr_chi2 / sum(lr_chi2),
      dom_pct  = 100 * McFadden_inc / sum(McFadden_inc),
      auc_full = auc_full
    )
  attr(out, "model") <- tibble(
    N = nrow(d), Events = sum(y == 1L),
    AUC_covars = auc_cov, AUC_full = auc_full,
    dAUC_block = auc_full - auc_cov, DeLong_p_block = delong_block
  )
  out
}

# ======================== 6. RUN ALL MODELS ==================================
results <- list(); models <- list()
for (def in names(components_by_definition)) {
  comp_map <- components_by_definition[[def]]; comps <- names(comp_map)
  for (out in names(outcomes)) {
    d <- dat %>% select(all_of(c(out, comps, covariates))) %>% tidyr::drop_na()
    if (nrow(d) < 1L || length(unique(d[[out]])) < 2L) next
    r <- analyse_components(d, out, comps, covariates)
    r <- r %>% mutate(
      Definition  = unname(definition_label[[def]]),
      Outcome     = unname(outcomes[[out]]),
      Outcome_var = out,
      Criterion   = unname(comp_map[Component_var])
    )
    results[[paste(def, out)]] <- r
    models[[paste(def, out)]] <- attr(r, "model") %>%
      mutate(Definition = unname(definition_label[[def]]), Outcome = unname(outcomes[[out]]), .before = 1)
  }
}
master <- bind_rows(results)
order_outcome <- function(x) factor(x, levels = unname(outcomes))

# ======================== 7. CONSOLE SUMMARIES ===============================
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))

summary_lr <- master %>%
  arrange(Definition, order_outcome(Outcome), desc(share_lr)) %>%
  transmute(Definition, Outcome, Criterion,
            `LR chi-square` = round(lr_chi2, 1), df = 1L,
            `p (LR)` = fmt_p(p_LR), `Share, %` = round(share_lr, 1))

summary_dominance <- master %>%
  arrange(Definition, order_outcome(Outcome), desc(dom_pct)) %>%
  transmute(Definition, Outcome, Criterion,
            `McFadden R2 increment` = round(McFadden_inc, 4),
            `Dominance weight, %` = round(dom_pct, 1))

summary_auc <- master %>%
  arrange(Definition, order_outcome(Outcome), desc(dAUC)) %>%
  transmute(Definition, Outcome, Criterion,
            `Adjusted OR (95% CI)`    = sprintf("%.2f (%.2f-%.2f)", OR, OR_lo, OR_hi),
            `AUC component alone`     = round(auc_alone, 3),
            `AUC full model`          = round(auc_full, 3),
            `AUC without component`   = round(auc_red, 3),
            `Delta AUC`               = round(dAUC, 3),
            `DeLong p (full vs reduced)` = fmt_p(delong_p))

summary_models <- bind_rows(models) %>%
  arrange(Definition, order_outcome(Outcome)) %>%
  transmute(Definition, Outcome, N, Events,
            `AUC covariates only` = round(AUC_covars, 3),
            `AUC full` = round(AUC_full, 3),
            `Delta AUC (component block)` = round(dAUC_block, 3),
            `DeLong p (block)` = fmt_p(DeLong_p_block))

cat("\n### All models adjusted for age and sex ###\n")
cat("\n===== Share of partial likelihood-ratio chi-square by component =====\n")
print(as.data.frame(summary_lr), row.names = FALSE)
cat("\n===== General dominance (McFadden R2) by component (Figure 2a) =====\n")
print(as.data.frame(summary_dominance), row.names = FALSE)
cat("\n===== Adjusted ORs (Figure 2b) and component-level discrimination =====\n")
print(as.data.frame(summary_auc), row.names = FALSE)
cat("\n===== Component block by definition and outcome =====\n")
print(as.data.frame(summary_models), row.names = FALSE)

# ======================== 8. FIGURE 2 ========================================
criterion_colours <- c(
  "Blood pressure"          = "#2171B5",
  "Glucose / diabetes"      = "#6A51A3",
  "Waist circumference"     = "#D7301F",
  "Central adiposity (WHR)" = "#F16913",
  "Triglycerides"           = "#238B45",
  "Low HDL"                 = "#A1D99B"
)
criterion_levels <- c("Blood pressure", "Glucose / diabetes", "Waist circumference",
                      "Central adiposity (WHR)", "Triglycerides", "Low HDL")

# Shared labels: outcome abbreviations and definitions as A / B
outcome_abbrev <- c(
  "Cardiovascular disease (CVD)" = "CVD",
  "Fatty liver disease (FLD)"    = "FLD",
  "Chronic kidney disease (CKD)" = "CKD"
)
definition_AB <- c(
  "Definition A" = "A",
  "Definition B" = "B"
)

# -- Panel a: relative contribution (general dominance, %) --------------------
# Stacked bars by definition within each outcome (CVD, FLD, CKD, left to
# right). Outcome labels are shown once, below panel b.
make_dominance_panel <- function(value) {
  df2 <- master %>%
    mutate(Pct = .data[[value]],
           Criterion  = factor(Criterion, levels = criterion_levels),
           Outcome    = factor(Outcome, levels = unname(outcomes)),
           Definition = factor(Definition, levels = c("Definition A", "Definition B")),
           bar_label  = ifelse(Pct >= 3, sprintf("%.0f", Pct), ""))
  ggplot(df2, aes(x = Definition, y = Pct, fill = Criterion)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.4,
             position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = bar_label), position = position_stack(vjust = 0.5, reverse = TRUE),
              color = "white", fontface = "bold", size = 3.6) +
    facet_wrap(
      ~ Outcome, ncol = 3,
      strip.position = "bottom",
      labeller = labeller(Outcome = outcome_abbrev)
    ) +
    scale_fill_manual(values = criterion_colours, name = "Criterion", breaks = criterion_levels) +
    scale_x_discrete(labels = definition_AB, expand = expansion(add = 0.45)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(title = NULL, subtitle = NULL, x = NULL, y = "Relative Contribution (%)") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 11, face = "bold"),
          axis.text.y = element_text(hjust = 1, margin = margin(r = 2)),
          axis.title.y = element_text(margin = margin(r = 2)),
          axis.ticks.length.y = unit(1.5, "pt"),
          legend.position = "left",
          legend.direction = "vertical",
          legend.justification = "center",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 10),
          legend.key.size = unit(11, "pt"),
          legend.margin = margin(r = 4)) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE, title.position = "top"))
}

fig_dominance <- make_dominance_panel("dom_pct")

# -- Panel b: adjusted odds ratios (95% CI) -----------------------------------
# Rows = definition (A top, B bottom); columns = outcome.
criterion_order_y <- c("Low HDL", "Triglycerides", "Blood pressure",
                       "Glucose / diabetes", "Central adiposity (WHR)", "Waist circumference")
definition_colours <- c("Definition A" = "gray55", "Definition B" = "black")

df_forest <- master %>%
  mutate(
    Criterion  = factor(Criterion, levels = criterion_order_y),
    Outcome    = factor(Outcome, levels = unname(outcomes)),
    Definition = factor(Definition, levels = c("Definition A", "Definition B"))
  )

fig_forest <- ggplot(df_forest, aes(x = OR, y = Criterion, color = Definition)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.22, linewidth = 0.6) +
  geom_point(size = 2.6) +
  facet_grid(
    Definition ~ Outcome,
    scales = "free_y",
    space  = "free_y",
    switch = "x",
    labeller = labeller(
      Outcome    = outcome_abbrev,
      Definition = definition_AB
    )
  ) +
  scale_color_manual(values = definition_colours, guide = "none") +
  scale_x_log10(breaks = c(0.5, 1, 2, 3)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) +
  labs(title = NULL, subtitle = NULL, x = "OR (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    strip.text.x     = element_text(face = "bold", size = 8.8 * 1.5),
    strip.background = element_rect(fill = "gray92", color = NA),
    strip.placement  = "outside",
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off")

# -- Composite Figure 2: a) dominance above b) forest ------------------------
# Panels are stacked directly with patchwork so that each outcome column of
# panel b falls below the corresponding pair of bars in panel a.
tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 15, hjust = 0, vjust = 1),
  plot.tag.position = c(0.005, 0.99)
)

fig_dominance_tagged <- fig_dominance +
  labs(tag = "a)") +
  tag_theme +
  theme(plot.margin = margin(10, 5, 5, 5))

fig_forest_tagged <- fig_forest +
  labs(tag = "b)") +
  tag_theme +
  theme(plot.margin = margin(10, 5, 5, 5))

figure_2 <-
  (fig_dominance_tagged / fig_forest_tagged) +
  patchwork::plot_layout(heights = c(0.60, 0.40))

# Text sizes are absolute (points). The published figure was rendered on an
# 11 x 11 inch canvas at 600 dpi; other canvas sizes change the relative
# size of the text.
print(figure_2)
