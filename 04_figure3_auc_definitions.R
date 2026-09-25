###############################################################################
# Figure 3. Discrimination of prevalent cardiovascular disease (CVD), fatty
# liver disease (FLD), and chronic kidney disease (CKD) by Definition A and
# Definition B, in the total population and by BMI category
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Analysis: unadjusted area under the receiver-operating-characteristic
# curve (AUC) of each binary definition (MUH = 1) for each outcome, with
# 95% CIs estimated by the DeLong method; AUCs of the two definitions are
# compared within the same participants with the paired DeLong test.
#
# Output (printed to the graphics device and console; no files are written):
#   - Figure 3: AUC (95% CI) by definition (rows A and B), stratum (total,
#     normal weight, overweight, obesity), and outcome (panels CVD, FLD, CKD)
#   - Console summary: AUC (95% CI) of each definition, difference
#     (Definition B minus Definition A), and paired DeLong p value
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs; MUH = 2 or more of 5 components
# Definition B: empirically derived definition (Zembic et al.);
#               MUH = 1 or more of 3 components
#
# Analytical population: BMI >= 18.5 kg/m2, sex recorded as female or male,
# and metabolic health status resolvable under both definitions (same as
# Table 1). Each analysis uses participants with available outcome data.
#
# R version 4.6.0
# Packages: dplyr, readr, stringr, pROC, ggplot2
###############################################################################

library(dplyr)
library(readr)
library(stringr)
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
clip_range <- function(x, lo, hi) {
  if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)
}

# -- Data preparation (same as Table 1) ---------------------------------------

dat <- raw_data %>%
  transmute(
    SEX_RAW = str_to_upper(trimws(as.character(PARTICIPANTE.SEXO))),
    SEX = case_when(
      SEX_RAW %in% c("H", "HOMBRE", "MASCULINO", "MALE", "MAN")         ~ "Male",
      SEX_RAW %in% c("M", "MUJER", "F", "FEMENINO", "FEMALE", "WOMAN") ~ "Female",
      TRUE ~ NA_character_
    ),

    BMI           = suppressWarnings(as.numeric(ANTROPOMETRIA.IMC)),
    WAIST         = suppressWarnings(as.numeric(ANTROPOMETRIA.CIRCUNFERENCIA_CINTURA)),
    WHR           = suppressWarnings(as.numeric(INBODY.WHR_Waist_Hip_Ratio)),
    TG            = suppressWarnings(as.numeric(METABOLICOS.TRIGLICERIDOS)),
    HDL           = suppressWarnings(as.numeric(METABOLICOS.HDL)),
    SBP           = suppressWarnings(as.numeric(METABOLICOS.SISTOLICA)),
    DBP           = suppressWarnings(as.numeric(METABOLICOS.DIASTOLICA)),
    FASTING_HOURS = suppressWarnings(as.numeric(ANTROPOMETRIA.HORAS_AYUNO)),
    GLUCOSE_RAW   = suppressWarnings(as.numeric(METABOLICOS.GLUCOSA)),

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

  # Prespecified physiological ranges
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
    GLUCOSE = if_else(
      !is.na(FASTING_HOURS) & FASTING_HOURS >= 8,
      GLUCOSE_RAW,
      NA_real_
    ),

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
      SEX == "Female" & !is.na(WHR) & WHR >= 0.95 ~ 1L,
      SEX == "Female" & !is.na(WHR) & WHR <  0.95 ~ 0L,
      SEX == "Male"   & !is.na(WHR) & WHR >= 1.03 ~ 1L,
      SEX == "Male"   & !is.na(WHR) & WHR <  1.03 ~ 0L,
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
      SEX == "Male"   & !is.na(WAIST) & WAIST >= 90 ~ 1L,
      SEX == "Male"   & !is.na(WAIST) & WAIST <  90 ~ 0L,
      SEX == "Female" & !is.na(WAIST) & WAIST >= 80 ~ 1L,
      SEX == "Female" & !is.na(WAIST) & WAIST <  80 ~ 0L,
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
      (SEX == "Male"   & !is.na(HDL) & HDL < 40) |
        (SEX == "Female" & !is.na(HDL) & HDL < 50) |
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

  # Analytical population (same as Table 1); a non-missing BMI category
  # implies valid BMI >= 18.5 kg/m2
  filter(
    !is.na(BMI_CATEGORY),
    !is.na(SEX),
    !is.na(defA_muh),
    !is.na(defB_muh)
  )

# -- AUC with 95% CI (DeLong) for one definition, outcome, and stratum -------

calc_auc_ci <- function(df, outcome_var, def_var, def_label, outcome_label, stratum_label) {

  d <- df %>%
    filter(
      !is.na(.data[[outcome_var]]),
      !is.na(.data[[def_var]])
    ) %>%
    mutate(
      Y = as.integer(.data[[outcome_var]]),
      X = as.integer(.data[[def_var]])
    )

  if (n_distinct(d$Y) < 2 | n_distinct(d$X) < 2 | nrow(d) < 20) return(NULL)

  roc_obj <- tryCatch(
    roc(d$Y, d$X, quiet = TRUE, ci = TRUE),
    error = function(e) NULL
  )

  if (is.null(roc_obj)) return(NULL)

  ci_vals <- tryCatch(
    ci.auc(roc_obj, method = "delong"),
    error = function(e) NULL
  )

  if (is.null(ci_vals)) return(NULL)

  tibble(
    Stratum    = stratum_label,
    Outcome    = outcome_label,
    Definition = def_label,
    AUC        = round(as.numeric(auc(roc_obj)), 3),
    CI_low     = round(ci_vals[1], 3),
    CI_high    = round(ci_vals[3], 3)
  )
}

# -- AUC table ----------------------------------------------------------------

outcome_list <- list(
  list(var = "CVD", label = "CVD"),
  list(var = "FLD", label = "FLD"),
  list(var = "CKD", label = "CKD")
)

definition_list <- list(
  list(var = "defA_muh", label = "A"),
  list(var = "defB_muh", label = "B")
)

strata <- list(
  list(label = "Total",         df = dat),
  list(label = "Normal weight", df = dat %>% filter(BMI_CATEGORY == "Normal weight")),
  list(label = "Overweight",    df = dat %>% filter(BMI_CATEGORY == "Overweight")),
  list(label = "Obesity",       df = dat %>% filter(BMI_CATEGORY == "Obesity"))
)

auc_table <- bind_rows(lapply(strata, function(s) {
  bind_rows(lapply(outcome_list, function(o) {
    bind_rows(lapply(definition_list, function(def) {
      calc_auc_ci(s$df, o$var, def$var, def$label, o$label, s$label)
    }))
  }))
})) %>%
  mutate(
    Stratum    = factor(Stratum, levels = c("Total", "Normal weight", "Overweight", "Obesity")),
    Outcome    = factor(Outcome, levels = sapply(outcome_list, `[[`, "label")),
    Definition = factor(Definition, levels = c("A", "B"))
  )

# -- Paired DeLong comparison: Definition A vs Definition B -------------------

compare_auc <- function(df, outcome_var, outcome_label, stratum_label) {

  d <- df %>%
    filter(
      !is.na(.data[[outcome_var]]),
      !is.na(defA_muh),
      !is.na(defB_muh)
    ) %>%
    mutate(Y = as.integer(.data[[outcome_var]]))

  if (n_distinct(d$Y) < 2 || nrow(d) < 20) return(NULL)

  roc_a <- tryCatch(
    roc(d$Y, as.numeric(d$defA_muh), quiet = TRUE, ci = TRUE),
    error = function(e) NULL
  )

  roc_b <- tryCatch(
    roc(d$Y, as.numeric(d$defB_muh), quiet = TRUE, ci = TRUE),
    error = function(e) NULL
  )

  if (is.null(roc_a) || is.null(roc_b)) return(NULL)

  ci_a <- as.numeric(ci.auc(roc_a, method = "delong"))
  ci_b <- as.numeric(ci.auc(roc_b, method = "delong"))

  auc_a <- as.numeric(auc(roc_a))
  auc_b <- as.numeric(auc(roc_b))

  test <- tryCatch(
    roc.test(roc_a, roc_b, method = "delong", paired = TRUE),
    error = function(e) NULL
  )

  tibble(
    Stratum         = stratum_label,
    Outcome         = outcome_label,
    N               = nrow(d),
    Events          = sum(d$Y == 1L),
    AUC_A           = sprintf("%.3f (%.3f-%.3f)", auc_a, ci_a[1], ci_a[3]),
    AUC_B           = sprintf("%.3f (%.3f-%.3f)", auc_b, ci_b[1], ci_b[3]),
    Diff_B_minus_A  = round(auc_b - auc_a, 3),
    p_DeLong        = if (is.null(test)) NA_real_ else signif(test$p.value, 3)
  )
}

comparison_table <- bind_rows(lapply(strata, function(s) {
  bind_rows(lapply(outcome_list, function(o) {
    compare_auc(s$df, o$var, o$label, s$label)
  }))
}))

cat("\n=== AUC by definition and paired DeLong comparison (Definition A vs Definition B) ===\n")
cat("Diff = AUC(Definition B) - AUC(Definition A)\n\n")
print(as.data.frame(comparison_table), row.names = FALSE)

# -- Figure 3 -----------------------------------------------------------------
# y axis: one row per definition (A top, B bottom); within each row, one
# point per stratum. x axis: AUC with 95% CI. Panels: CVD, FLD, CKD.

stratum_offset <- c(
  "Total"         =  0.30,
  "Normal weight" =  0.10,
  "Overweight"    = -0.10,
  "Obesity"       = -0.30
)

df_plot <- auc_table %>%
  mutate(
    y_base = if_else(Definition == "A", 2, 1),
    y_pos  = y_base + stratum_offset[as.character(Stratum)]
  )

y_definition_labels <- df_plot %>%
  group_by(Definition) %>%
  summarise(y_mid = mean(y_pos), .groups = "drop")

stratum_colours <- c(
  "Total"         = "black",
  "Normal weight" = "#2171B5",
  "Overweight"    = "#238B45",
  "Obesity"       = "#CB181D"
)

stratum_shapes <- c(
  "Total"         = 18,
  "Normal weight" = 16,
  "Overweight"    = 17,
  "Obesity"       = 15
)

figure_3 <- ggplot(
  df_plot,
  aes(
    x = AUC,
    y = y_pos,
    color = Stratum,
    shape = Stratum
  )
) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.64
  ) +
  geom_hline(
    yintercept = 1.5,
    color = "gray80",
    linewidth = 0.64,
    inherit.aes = FALSE
  ) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.10,
    linewidth = 0.96
  ) +
  geom_point(aes(size = Stratum)) +
  facet_wrap(~ Outcome, nrow = 1, scales = "free_x") +
  scale_color_manual(values = stratum_colours, name = "BMI category") +
  scale_shape_manual(values = stratum_shapes, name = "BMI category") +
  scale_size_manual(
    values = c(
      "Total"         = 5.6,
      "Normal weight" = 4.0,
      "Overweight"    = 4.0,
      "Obesity"       = 4.0
    ),
    guide = "none"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 7, linewidth = 1.4)),
    shape = guide_legend(override.aes = list(size = 7, linewidth = 1.4))
  ) +
  scale_x_continuous(
    breaks = seq(0.45, 0.80, by = 0.05),
    labels = function(x) sprintf("%.2f", x)
  ) +
  scale_y_continuous(
    breaks = y_definition_labels$y_mid,
    labels = y_definition_labels$Definition,
    limits = c(0.5, 2.7)
  ) +
  labs(
    title = NULL,
    x     = "AUC (95% CI)",
    y     = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 13,
      color = "gray20"
    ),
    strip.background = element_rect(
      fill = "gray95",
      color = NA
    ),
    axis.text.x = element_text(
      size = 13,
      color = "gray20",
      angle = 30,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 13,
      color = "gray20",
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 14,
      margin = margin(t = 6)
    ),
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.key.size = grid::unit(1.1, "cm"),
    panel.grid.major.x = element_line(
      color = "gray92",
      linewidth = 0.3
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.margin = margin(24, 32, 24, 32)
  )

print(figure_3)
