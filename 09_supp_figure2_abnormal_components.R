###############################################################################
# Supplementary Figure 2. Combinations of abnormal components under
# Definition A and Definition B
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# A component is considered abnormal when the participant meets its
# criterion (eg, waist circumference >= 90 cm in men or >= 80 cm in women).
# Percentages are calculated over the analytical sample.
#
# Output (printed to the graphics device and console; no files are written):
#   - Supplementary Figure 2, panels a-b (UpSet-type plots):
#       (a) combinations of abnormal components under Definition A
#       (b) combinations of abnormal components under Definition B
#     In each panel, the upper bars show the percentage of participants
#     with each combination (participants with no abnormal component are
#     not shown), the dot matrix indicates the components included in each
#     combination, and the left bars show the percentage of participants
#     with each abnormal component.
#   - Console summaries: analytical sample size, prevalence of each
#     abnormal component by definition (ranked), and frequency of each
#     combination of abnormal components by definition
#
# Data cleaning, physiological ranges, and definitions are identical to
# those in 01_table1_characteristics.R.
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs; MUH = 2 or more of 5 components
# Definition B: empirically derived definition (Zembic and colleagues);
#               MUH = 1 or more of 3 components
#
# Analytical population: BMI >= 18.5 kg/m2, sex recorded as female or male,
# and metabolic health status resolvable under both definitions.
#
# R version 4.6.0
# Packages: dplyr, readr, stringr, tibble, ggplot2, patchwork
###############################################################################

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

# Path to the oriGen questionnaire data file (tab-separated).
# Column names correspond to the original oriGen data release.
data_path <- "path/to/oriGen_questionnaire_data.tsv"
raw_data  <- read_tsv(data_path, show_col_types = FALSE)

# Wide console output so that the summary tables print on one line per row
options(width = 250)

# -- Figure settings ----------------------------------------------------------

INCLUDE_ZERO   <- FALSE      # TRUE adds a bar for participants with no abnormal component
SHOW_SET_SIZE  <- TRUE       # TRUE adds the left panel (% with each abnormal component)
BAR_COLOUR     <- "#9E9AC8"
DOT_PRESENT    <- "gray10"
DOT_ABSENT     <- "gray85"

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
clip_range <- function(x, lo, hi) if_else(!is.na(x) & x >= lo & x <= hi, x, NA_real_)

fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

# -- Data preparation (same as 01_table1_characteristics.R) -------------------

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
    lipid_lowering   = safe_bool(HISTORIA_MEDICA.TRATA_MEDICAMENTO)
  ) %>%

  # Prespecified physiological ranges (same as in the Table 1 script)
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

    # -------------------------------------------------------------------------
    # METABOLIC HEALTH DEFINITIONS (same as Table 1)
    #
    # Component coding: 1 = present (abnormal); 0 = absent.
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
    )
  ) %>%

  # Analytical population (same as Table 1)
  filter(
    !is.na(BMI),
    BMI >= 18.5,
    !is.na(SEX),
    !is.na(defA_muh),
    !is.na(defB_muh)
  )

n_total <- nrow(dat)

cat("\n=== Analytical sample ===\n")
cat("N =", fmt_n(n_total), "\n")

# -- Components of each definition (variable -> label) ------------------------

components_A <- c(
  defA_waist   = "Elevated waist circumference",
  defA_tg      = "Elevated triglycerides",
  defA_hdl     = "Low HDL cholesterol",
  defA_bp      = "Elevated blood pressure",
  defA_glucose = "Elevated glucose or diabetes"
)
components_B <- c(
  defB_whr      = "Elevated waist-to-hip ratio",
  defB_bp       = "Elevated systolic blood pressure",
  defB_diabetes = "Diabetes"
)

# Internal check: all components are non-missing in the analytical sample.
stopifnot(
  !anyNA(dat[, names(components_A)]),
  !anyNA(dat[, names(components_B)])
)

# -- Combinations of abnormal components --------------------------------------
# One row per observed combination, with its frequency and a 0/1 column for
# each component.

count_combinations <- function(comps) {
  m   <- as.data.frame(dat[, names(comps)])
  key <- do.call(paste0, m)
  out <- tibble(key = key) %>%
    count(key, name = "n") %>%
    mutate(pct = 100 * n / n_total)
  for (j in seq_along(comps)) {
    out[[names(comps)[j]]] <- as.integer(substr(out$key, j, j))
  }
  out$n_components <- rowSums(out[, names(comps), drop = FALSE])
  out$Combination <- apply(out[, names(comps), drop = FALSE], 1, function(r) {
    s <- unname(comps[names(comps)[r == 1L]])
    if (length(s) == 0) "None" else paste(s, collapse = " + ")
  })
  out %>% arrange(desc(n))
}

combinations_A <- count_combinations(components_A)
combinations_B <- count_combinations(components_B)

# -- Console summaries --------------------------------------------------------

# Prevalence of each abnormal component, ranked within each definition
component_prevalence <- bind_rows(
  tibble(Definition = "A", Component = unname(components_A),
         n = vapply(names(components_A), function(v) sum(dat[[v]] == 1L), numeric(1))),
  tibble(Definition = "B", Component = unname(components_B),
         n = vapply(names(components_B), function(v) sum(dat[[v]] == 1L), numeric(1)))
) %>%
  group_by(Definition) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  mutate(Rank = row_number(),
         `n (%)` = sprintf("%s (%.1f%%)", fmt_n(n), 100 * n / n_total)) %>%
  ungroup() %>%
  select(Definition, Rank, Component, `n (%)`)

cat("\n=== Prevalence of each abnormal component, by definition ===\n")
print(as.data.frame(component_prevalence), row.names = FALSE, right = FALSE)

print_combinations <- function(comb, title) {
  cat("\n===", title, "===\n")
  print(as.data.frame(
    comb %>% transmute(
      Rank = row_number(),
      Combination,
      `Number of components` = n_components,
      `n (%)` = sprintf("%s (%.1f%%)", fmt_n(n), pct)
    )
  ), row.names = FALSE, right = FALSE)
}

print_combinations(combinations_A,
                   "Combinations of abnormal components, Definition A (by frequency)")
print_combinations(combinations_B,
                   "Combinations of abnormal components, Definition B (by frequency)")

# -- UpSet-type plot (one panel per definition) -------------------------------

make_upset <- function(comb, comps, include_zero = INCLUDE_ZERO,
                       show_set_size = SHOW_SET_SIZE, label_size = 2.6) {

  df <- comb
  if (!include_zero) df <- df %>% filter(n_components > 0)
  df <- df %>% arrange(desc(n)) %>% mutate(x = row_number())
  x_limits <- c(0.5, nrow(df) + 0.5)

  # Upper panel: percentage of participants with each combination
  p_bars <- ggplot(df, aes(x = x, y = pct)) +
    geom_col(fill = BAR_COLOUR, color = "gray30", linewidth = 0.25, width = 0.7) +
    geom_text(aes(label = paste0(fmt_n(n), "\n", sprintf("%.1f%%", pct))),
              vjust = -0.25, size = label_size, lineheight = 0.85, color = "gray20") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(labels = function(v) paste0(v, "%"),
                       expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = "Participants (%)") +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x  = element_blank(),
          axis.title.y = element_text(size = 11, margin = margin(r = 6)),
          plot.margin  = margin(5, 5, 0, 5))

  # Dot matrix: components included in each combination
  y_levels <- rev(unname(comps))
  matrix_df <- bind_rows(lapply(names(comps), function(v) {
    tibble(x = df$x, Component = comps[[v]], present = df[[v]] == 1L)
  })) %>%
    mutate(Component = factor(Component, levels = y_levels))

  segments <- matrix_df %>%
    filter(present) %>%
    mutate(yy = as.integer(Component)) %>%
    group_by(x) %>%
    summarise(y_min = min(yy), y_max = max(yy), .groups = "drop") %>%
    filter(y_max > y_min)

  stripes <- tibble(yy = seq_along(y_levels)) %>% filter(yy %% 2 == 0)

  p_matrix <- ggplot(matrix_df, aes(x = x, y = Component)) +
    geom_rect(data = stripes, inherit.aes = FALSE,
              aes(xmin = x_limits[1], xmax = x_limits[2],
                  ymin = yy - 0.5, ymax = yy + 0.5),
              fill = "gray96") +
    geom_segment(data = segments, inherit.aes = FALSE,
                 aes(x = x, xend = x, y = y_min, yend = y_max),
                 linewidth = 0.8, color = DOT_PRESENT) +
    geom_point(aes(color = present), size = 3) +
    scale_color_manual(values = c(`TRUE` = DOT_PRESENT, `FALSE` = DOT_ABSENT),
                       guide = "none") +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(),
          panel.grid  = element_blank(),
          axis.text.y = element_text(size = 10, color = "gray15"),
          plot.margin = margin(0, 5, 5, 5))

  matrix_height <- 0.28 * length(comps) + 0.2

  if (!show_set_size) {
    return(p_bars / p_matrix + plot_layout(heights = c(2.2, matrix_height)))
  }

  # Left panel: percentage of participants with each abnormal component
  set_size <- tibble(
    Component = factor(unname(comps), levels = y_levels),
    pct = vapply(names(comps), function(v) 100 * sum(dat[[v]] == 1L) / n_total, numeric(1))
  )

  p_set_size <- ggplot(set_size, aes(x = pct, y = Component)) +
    geom_col(fill = "gray55", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              hjust = 1.1, size = 2.8, color = "gray20") +
    scale_x_reverse(expand = expansion(mult = c(0.35, 0)),
                    labels = function(v) paste0(v, "%")) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Participants with\nabnormal component (%)", y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y  = element_blank(),
          panel.grid   = element_blank(),
          axis.title.x = element_text(size = 8.5),
          axis.text.x  = element_text(size = 7.5),
          plot.margin  = margin(0, 0, 5, 5))

  layout <- "
  AB
  CD
  "
  plot_spacer() + p_bars + p_set_size + p_matrix +
    plot_layout(design = layout, widths = c(1.1, 5),
                heights = c(2.2, matrix_height))
}

upset_A <- make_upset(combinations_A, components_A)
upset_B <- make_upset(combinations_B, components_B, label_size = 3.2)

# =============================================================================
# COMPOSITE SUPPLEMENTARY FIGURE 2: a) Definition A above b) Definition B
# =============================================================================

add_panel_tag <- function(p, tag) {
  wrap_elements(full = p) +
    labs(tag = tag) +
    theme(plot.tag          = element_text(size = 14, face = "bold"),
          plot.tag.position = c(0.005, 0.985))
}

supp_figure_2 <- add_panel_tag(upset_A, "a") / add_panel_tag(upset_B, "b") +
  plot_layout(heights = c(1.15, 1))

# Text sizes are absolute (points). The published figure was rendered on a
# 14 x 12 inch canvas; other canvas sizes change the relative size of the
# text.
print(supp_figure_2)
