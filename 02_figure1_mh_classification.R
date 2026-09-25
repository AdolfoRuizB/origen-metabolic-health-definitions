###############################################################################
# Figure 1. Metabolic health classification under Definition A and
# Definition B, overall and by BMI category
#
# Data: oriGen Project, data release of July 2026. Individual participant
# data are not included in this repository; they are available upon
# specific request to the oriGen Project.
#
# Output (printed to the graphics device and console; no files are written):
#   - Figure 1, panels a-e:
#       (a) prevalence of metabolically healthy (MH) and metabolically
#           unhealthy (MUH) participants under Definition A and
#           Definition B, by BMI category and in the total population
#       (b-e) reclassification of participants between definitions
#           (alluvial plots): normal weight (b), overweight (c),
#           obesity (d), and total population (e)
#   - Console summaries: analytical sample sizes, prevalence of MH/MUH,
#     and reclassification between definitions by BMI category
#
# Definition A: harmonised metabolic syndrome criteria with Latin American
#               waist-circumference cut-offs; MUH = 2 or more of 5 components
# Definition B: empirically derived definition (Zembic et al.);
#               MUH = 1 or more of 3 components
#
# Analytical population: BMI >= 18.5 kg/m2, sex recorded as female or male,
# and metabolic health status resolvable under both definitions.
#
# R version 4.6.0
# Packages: dplyr, readr, stringr, tidyr, ggplot2 (>= 3.5.0), ggalluvial,
#           patchwork
###############################################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(patchwork)

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

# -- Data preparation ---------------------------------------------------------

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
    diabetes                = safe_bool(HISTORIA_MEDICA.DIABETES),
    diabetes_type1          = safe_bool(HISTORIA_MEDICA.DIABETES_TIPO1),
    diabetes_type2          = safe_bool(HISTORIA_MEDICA.DIABETES_TIPO2),
    glucose_lowering_oral   = safe_bool(HISTORIA_MEDICA.PASTILL_CONTROL_AZUCAR),
    glucose_lowering_insulin = safe_bool(HISTORIA_MEDICA.INSULINA_CONTROL_AZUCAR),
    glucose_lowering_any    = safe_bool(HISTORIA_MEDICA.TRATAMIENTO_AZUCAR),

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

    BMI_CATEGORY = case_when(
      BMI >= 18.5 & BMI < 25 ~ "Normal weight",
      BMI >= 25   & BMI < 30 ~ "Overweight",
      BMI >= 30              ~ "Obesity",
      TRUE                   ~ NA_character_
    ),

    # -------------------------------------------------------------------------
    # METABOLIC HEALTH DEFINITIONS
    #
    # Component coding: 1 = present; 0 = absent.
    # Handling of partially missing data within a component:
    #   - 1 if any of its evaluable criteria is met (measured value meets the
    #     cut-off, or the corresponding diagnosis or treatment is reported);
    #   - 0 if none is met and at least one criterion is evaluable;
    #   - NA only if all criteria are missing, uninterpretable, or out of range.
    # Single-variable components are NA when that variable is missing.
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

  # Analytical population (same as Table 1)
  filter(
    !is.na(BMI),
    BMI >= 18.5,
    !is.na(SEX),
    !is.na(defA_muh),
    !is.na(defB_muh)
  ) %>%

  mutate(
    BMI_CATEGORY = factor(
      BMI_CATEGORY,
      levels = c("Normal weight", "Overweight", "Obesity")
    )
  )

# -- Analytical sample sizes (should match Table 1) --------------------------

sample_sizes <- tibble(
  Group = c(
    "Total",
    "Definition A, MH",
    "Definition A, MUH",
    "Definition B, MH",
    "Definition B, MUH",
    "Normal weight",
    "Overweight",
    "Obesity"
  ),
  N = c(
    nrow(dat),
    sum(dat$defA_muh == 0L),
    sum(dat$defA_muh == 1L),
    sum(dat$defB_muh == 0L),
    sum(dat$defB_muh == 1L),
    sum(dat$BMI_CATEGORY == "Normal weight"),
    sum(dat$BMI_CATEGORY == "Overweight"),
    sum(dat$BMI_CATEGORY == "Obesity")
  )
)

cat("\n=== Analytical sample sizes ===\n")
print(sample_sizes, n = Inf)

# Internal check: categories must add up to the total.
stopifnot(
  sum(dat$defA_muh == 0L) + sum(dat$defA_muh == 1L) == nrow(dat),
  sum(dat$defB_muh == 0L) + sum(dat$defB_muh == 1L) == nrow(dat),
  sum(!is.na(dat$BMI_CATEGORY)) == nrow(dat)
)

# -- Strata ------------------------------------------------------------------

strata <- list(
  list(label = "Total",         df = dat),
  list(label = "Normal weight", df = dat %>% filter(BMI_CATEGORY == "Normal weight")),
  list(label = "Overweight",    df = dat %>% filter(BMI_CATEGORY == "Overweight")),
  list(label = "Obesity",       df = dat %>% filter(BMI_CATEGORY == "Obesity"))
)

# -- Colours and axis labels ------------------------------------------------

status_colours    <- c("MH" = "#6ABF69", "MUH" = "#F4A460")
definition_labels <- c("A", "B")

# -- Alluvial plot (Figure 1b-e) ----------------------------------------------
# y axis: percentage of participants in the stratum (0-100).

make_alluvial <- function(df) {

  n_total <- nrow(df)

  df_alluvial <- df %>%
    mutate(
      defA_status = if_else(defA_muh == 0L, "MH", "MUH"),
      defB_status = if_else(defB_muh == 0L, "MH", "MUH")
    ) %>%
    count(defA_status, defB_status, name = "n") %>%
    mutate(pct = 100 * n / n_total)   # unrounded, so the stack sums to 100

  ggplot(df_alluvial,
         aes(axis1 = defA_status, axis2 = defB_status, y = pct)) +
    geom_alluvium(aes(fill = defA_status), alpha = 0.6, width = 1/3) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/3,
                 color = "white", linewidth = 0.6) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),
              size = 4.2, fontface = "bold", color = "white") +
    scale_x_discrete(
      limits = c("defA_status", "defB_status"),
      labels = definition_labels,
      expand = c(0.25, 0.1)
    ) +
    scale_fill_manual(values = status_colours, guide = "none") +
    scale_y_continuous(
      breaks = seq(0, 100, by = 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = NULL, subtitle = NULL, x = NULL, y = "Participants (%)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x        = element_text(size = 12, face = "bold", color = "gray20"),
      axis.text.y        = element_text(size = 9,  color = "gray20"),
      axis.title.y       = element_text(size = 10, angle = 90, margin = margin(r = 8)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.margin        = margin(15, 25, 15, 25)
    )
}

alluvial_plots <- lapply(strata, function(s) make_alluvial(s$df))
names(alluvial_plots) <- vapply(strata, function(s) s$label, character(1))

# -- Prevalence bar chart (Figure 1a) ----------------------------------------

panel_order <- c("Normal weight", "Overweight", "Obesity", "Total")

prevalence_data <- bind_rows(lapply(strata, function(s) {
  df <- s$df
  n  <- nrow(df)
  bind_rows(
    tibble(Stratum = s$label, Definition = "A", Status = "MH",  Pct = round(100 * sum(df$defA_muh == 0L) / n, 1)),
    tibble(Stratum = s$label, Definition = "A", Status = "MUH", Pct = round(100 * sum(df$defA_muh == 1L) / n, 1)),
    tibble(Stratum = s$label, Definition = "B", Status = "MH",  Pct = round(100 * sum(df$defB_muh == 0L) / n, 1)),
    tibble(Stratum = s$label, Definition = "B", Status = "MUH", Pct = round(100 * sum(df$defB_muh == 1L) / n, 1))
  )
})) %>%
  mutate(
    Stratum    = factor(Stratum,    levels = panel_order),
    Definition = factor(Definition, levels = c("A", "B")),
    Status     = factor(Status,     levels = c("MH", "MUH"))
  )

p_prevalence <- ggplot(prevalence_data,
                       aes(x = Definition, y = Pct, fill = Status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "white") +
  geom_text(
    aes(label = format(Pct, trim = TRUE)),
    position = position_dodge(width = 0.7),
    vjust = -0.5, size = 4.2, fontface = "bold", color = "gray15"
  ) +
  facet_wrap(~ Stratum, nrow = 1, strip.position = "bottom") +
  scale_x_discrete(labels = definition_labels) +
  scale_fill_manual(values = status_colours, name = "Metabolic status") +
  scale_y_continuous(
    limits = c(0, 110),                      # headroom for bar labels
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(title = NULL, x = NULL, y = "Participants (%)") +
  theme_minimal(base_size = 11) +
  theme(
    strip.placement    = "outside",
    strip.text         = element_text(face = "bold", size = 10, color = "gray20",
                                      margin = margin(t = 5, b = 5)),
    strip.background   = element_rect(fill = "gray95", color = NA),
    axis.text.x        = element_text(size = 12, face = "bold", color = "gray20",
                                      margin = margin(t = 4)),
    axis.text.y        = element_text(size = 9,  color = "gray20"),
    axis.title.y       = element_text(size = 10, angle = 90, margin = margin(r = 8)),
    legend.title       = element_text(face = "bold", size = 10),
    legend.text        = element_text(size = 9),
    legend.position    = "bottom",
    panel.spacing.x    = grid::unit(0.6, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(15, 20, 15, 20)
  )

# =============================================================================
# COMPOSITE FIGURE 1 (panels a-e)
# =============================================================================
# Layout: panel a spans the top row; each alluvial plot (b-e) sits below its
# corresponding facet of panel a. Stratum labels are drawn once, below the
# bottom row, and label each full column.
#
#   +-------------+-------------+-------------+-------------+
#   | a)                                                    |
#   +-------------+-------------+-------------+-------------+
#   |     b)      |     c)      |     d)      |     e)      |
#   +-------------+-------------+-------------+-------------+
#   |Normal weight| Overweight  |   Obesity   |    Total    |
#   +-------------+-------------+-------------+-------------+
#   |            Metabolic status: MH / MUH                 |
#   +-------------------------------------------------------+

# Adds the panel tag (a-e).
add_panel_tag <- function(p, tag) {
  p +
    labs(tag = tag) +
    theme(
      plot.tag          = element_text(size = 14, face = "bold", color = "gray15"),
      plot.tag.position = "topleft",
      plot.margin       = margin(t = 12, r = 10, b = 8, l = 10)
    )
}

# Removes the y axis (shared scale; drawn only in the leftmost panel, which
# also keeps the bottom panels aligned with the facets of panel a).
drop_y_axis <- function(p) {
  p + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
}

# Removes the facet strips from panel a (labels are drawn below the bottom row).
drop_strips <- function(p) {
  p + theme(
    strip.text       = element_blank(),
    strip.background = element_blank()
  )
}

# Adds a grey stratum label below an alluvial plot, using the same strip
# style as panel a.
add_bottom_strip <- function(p, label) {
  p_data <- p$data
  p_data[["Panel_label"]] <- label
  (p %+% p_data) +
    facet_wrap(~ Panel_label, strip.position = "bottom") +
    theme(
      strip.placement  = "outside",
      strip.text       = element_text(face = "bold", size = 15, color = "gray20",
                                      margin = margin(t = 7, b = 7)),
      strip.background = element_rect(fill = "gray95", color = NA),
      panel.border     = element_blank()
    )
}

p_a <- add_panel_tag(drop_strips(p_prevalence), "a)")

p_b <- add_panel_tag(
  add_bottom_strip(alluvial_plots[["Normal weight"]], "Normal weight"), "b)")

p_c <- drop_y_axis(add_panel_tag(
  add_bottom_strip(alluvial_plots[["Overweight"]],    "Overweight"),    "c)"))

p_d <- drop_y_axis(add_panel_tag(
  add_bottom_strip(alluvial_plots[["Obesity"]],       "Obesity"),       "d)"))

p_e <- drop_y_axis(add_panel_tag(
  add_bottom_strip(alluvial_plots[["Total"]],         "Total"),         "e)"))

figure_layout <- "
AAAA
BCDE
"

figure_1 <- p_a + p_b + p_c + p_d + p_e +
  plot_layout(design = figure_layout, heights = c(1, 1), guides = "collect") &
  theme(
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.title       = element_text(face = "bold", size = 14),
    legend.text        = element_text(size = 12),
    legend.key.size    = grid::unit(1.6, "lines"),
    legend.key.spacing = grid::unit(0.6, "lines"),   # requires ggplot2 >= 3.5.0
    legend.margin      = margin(t = 10, b = 4),
    plot.background    = element_rect(fill = "white", color = NA)
  )

# Text sizes are absolute (points). The published figure was rendered on a
# 17 x 9.5 inch canvas at 300 dpi; other canvas sizes change the relative
# size of the text.
print(figure_1)

# -- Console summaries (numbers reported in the Results) ---------------------

fmt_np <- function(n, N) paste0(formatC(n, format = "d", big.mark = ","),
                                " (", round(100 * n / N, 1), "%)")

# Prevalence of MH/MUH by definition and BMI category
prevalence_summary <- bind_rows(lapply(strata, function(s) {
  df <- s$df; N <- nrow(df)
  tibble(
    Stratum  = s$label,
    N        = N,
    defA_MH  = fmt_np(sum(df$defA_muh == 0L), N),
    defA_MUH = fmt_np(sum(df$defA_muh == 1L), N),
    defB_MH  = fmt_np(sum(df$defB_muh == 0L), N),
    defB_MUH = fmt_np(sum(df$defB_muh == 1L), N)
  )
}))

# Reclassification between Definition A and Definition B by BMI category
reclassification_summary <- bind_rows(lapply(strata, function(s) {
  df <- s$df; N <- nrow(df)
  both_mh      <- sum(df$defA_muh == 0L & df$defB_muh == 0L)
  both_muh     <- sum(df$defA_muh == 1L & df$defB_muh == 1L)
  a_mh_b_muh   <- sum(df$defA_muh == 0L & df$defB_muh == 1L)
  a_muh_b_mh   <- sum(df$defA_muh == 1L & df$defB_muh == 0L)
  tibble(
    Stratum         = s$label,
    N               = N,
    Both_MH         = fmt_np(both_mh, N),
    A_MH_to_B_MUH   = fmt_np(a_mh_b_muh, N),
    A_MUH_to_B_MH   = fmt_np(a_muh_b_mh, N),
    Both_MUH        = fmt_np(both_muh, N),
    Concordance_pct = round(100 * (both_mh + both_muh) / N, 1)
  )
}))

cat("\n=== Prevalence of MH/MUH by definition and BMI category ===\n")
print(as.data.frame(prevalence_summary), row.names = FALSE)

cat("\n=== Reclassification between Definition A and Definition B by BMI category ===\n")
print(as.data.frame(reclassification_summary), row.names = FALSE)
