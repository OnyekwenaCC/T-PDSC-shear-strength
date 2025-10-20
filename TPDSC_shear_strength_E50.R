# TPDSC: Shear strength + E50 

# ---- Packages ----
library(readxl)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# ---- User switches ----
excel_path        <- "strain_data.xlsx"   # input
export_dir        <- "TPDSC_out_from_excel"
export_inspection <- FALSE                # Set TRUE for detailed per-sample files

dir.create(export_dir, showWarnings = FALSE)

# ---- Model parameters ----
# qf_ref are experimental peak deviator stresses per (t, s3)
df_params <- data.frame(
  t      = c(3,3,3, 9,9,9, 11,11,11, 14,14,14),
  s3     = c(100,200,300, 100,200,300, 100,200,300, 100,200,300),
  A      = c(1.49E+01, 3.45E+01, 3.50E+01, 4.01E-02, 1.33E-02, 2.80E-01,
             3.30E-02, 1.26E-02, 1.40E+01, 1.75E-01, 6.76E-01, 6.94E-01),
  B      = c(4.16E-02, 7.39E-02, 1.99E-02, 1.75E-06, 6.09E-07, 1.05E-05,
             1.35E-06, 4.76E-07, 3.72E-04, 1.16E-06, 4.35E-06, 6.25E-06),
  alpha  = c(7.62E-01, 4.12E-01, 1.11E+00, 4.37E+00, 4.94E+00, 4.06E+00,
             4.67E+00, 5.13E+00, 2.76E+00, 3.54E+00, 3.28E+00, 3.48E+00),
  F      = c(2.05E-01, 5.23E-01, 4.70E-02, 1.19E-02, 1.30E-02, 1.28E-02,
             3.97E-03, 8.11E-03, 1.06E-02, 5.49E-03, 6.53E-03, 9.86E-03),
  i      = c(9.92E-01, 1.22E+00, 7.61E-01, 3.73E-01, 4.24E-01, 4.73E-01,
             2.26E-01, 3.63E-01, 4.20E-01, 1.82E-01, 3.78E-01, 5.26E-01),
  eps0   = c(1.12E-04, 2.86E-04, -3.18E-05, 4.45E-05, -3.15E-04, -3.06E-04,
             -4.20E-04, -4.03E-04, -3.17E-04, 1.40E-03, 9.46E-05, -9.30E-05),
  qf_ref = c(54.85884, 86.64505, 118.14998,
             671.5017, 873.50757, 1203.4559,
             879.0721, 1162.472, 1505.1271,
             1464.8936, 2107.8298, 2513.0015)
)

# ---- TPDSC Eq. (13) ----
sigma_tpds <- function(eps, A, B, alpha, F, i, eps0, qf_ref) {
  g <- (A * eps) / (B + eps^alpha)
  h <- 1 - (F / (eps + eps0))
  g * h + i * qf_ref
}

# ---- Read and normalize strains ----
strain_df <- read_excel(excel_path, sheet = 1)
names(strain_df) <- tolower(gsub("\\s+", "_", names(strain_df)))

map_col <- function(df, candidates, target){
  hit <- intersect(candidates, names(df))
  if (length(hit) > 0 && !(target %in% names(df))) {
    names(df)[match(hit[1], names(df))] <- target
  }
  df
}
strain_df <- strain_df |>
  map_col(c("sample_id","sample","id","series","specimen"), "sample_id") |>
  map_col(c("t","time","day","days","curing","curing_time"), "t") |>
  map_col(c("s3","sigma3","sigma_3","conf","confining","cell","s3_kpa"), "s3") |>
  map_col(c("eps","strain","epsilon","e","axial_strain","eps_decimal","eps_percent","eps_%"), "eps")

need <- c("t","s3","eps")
missing <- setdiff(need, names(strain_df))
if (length(missing) > 0) {
  stop(paste0("Excel missing columns: ", paste(missing, collapse = ", "),
              "\nPresent: ", paste(names(strain_df), collapse = ", ")))
}

if (!"sample_id" %in% names(strain_df)) {
  strain_df <- strain_df |>
    group_by(t, s3) |>
    mutate(sample_id = cur_group_id()) |>
    ungroup()
}

strain_df <- strain_df |>
  mutate(
    t = as.numeric(t),
    s3 = as.numeric(s3),
    eps = as.numeric(eps)
  ) |>
  arrange(sample_id, t, s3, eps)

if (is.finite(max(strain_df$eps, na.rm = TRUE)) &&
    max(strain_df$eps, na.rm = TRUE) > 1 && max(strain_df$eps, na.rm = TRUE) <= 100) {
  strain_df$eps <- strain_df$eps / 100
}

# ---- Predict σ–ε at provided strains (TPDSC) ----
pred_df <- strain_df |>
  left_join(df_params, by = c("t","s3")) |>
  mutate(q = sigma_tpds(eps, A, B, alpha, F, i, eps0, qf_ref)) |>
  select(sample_id, t, s3, eps, q)

# ---- Peak detection: qf, ε_peak, and pf (TPDSC) ----
peak_tbl <- pred_df |>
  group_by(sample_id, t, s3) |>
  slice_max(q, with_ties = FALSE) |>
  transmute(sample_id, t, s3, qf = q, eps_peak = eps, pf = (qf + 3*s3)/3) |>
  ungroup()

# ---- Secant modulus E50 from pre-peak branch (TPDSC) ----
E50_tbl <- pred_df |>
  left_join(peak_tbl, by = c("sample_id","t","s3")) |>
  filter(eps <= eps_peak) |>
  group_by(sample_id, t, s3) |>
  arrange(eps, .by_group = TRUE) |>
  do({
    qvec <- .$q; epsvec <- .$eps; qf <- unique(.$qf)
    q50 <- 0.5 * qf
    if (min(qvec) <= q50 && max(qvec) >= q50) {
      eps50 <- approx(x = qvec, y = epsvec, xout = q50, ties = "ordered")$y
      E50 <- q50 / eps50
    } else {
      eps50 <- NA_real_; E50 <- NA_real_
    }
    tibble(q50 = q50, eps50 = eps50, E50 = E50)
  }) |>
  ungroup()

peak_E50_tbl <- peak_tbl |>
  left_join(E50_tbl, by = c("sample_id","t","s3")) |>
  relocate(qf, eps_peak, q50, eps50, E50)

# ---- Aggregate pf–qf by (t, s3) for MC fit (TPDSC) ----
pfqf_agg_model <- peak_tbl |>
  group_by(t, s3) |>
  summarise(qf = mean(qf), pf = mean(pf), eps_peak = mean(eps_peak),
            n_samples = n(), .groups = "drop")

# ---- Mohr–Coulomb from p–q regression ----
fit_mc_pq_by_time <- function(pfqf) {
  groups <- split(pfqf, pfqf$t)
  lst <- lapply(groups, function(g) {
    fit <- lm(qf ~ pf, data = g)
    co  <- coef(summary(fit))
    M   <- unname(co["pf","Estimate"])
    k   <- unname(co["(Intercept)","Estimate"])
    sinphi <- (3*M) / (6 + M); sinphi <- max(min(sinphi, 0.999999), -0.999999)
    phi <- asin(sinphi) * 180/pi
    c   <- k * (3 - sinphi) / (6 * cos(phi*pi/180))
    data.frame(
      t = g$t[1], M = M, k = k, phi_deg = phi, c_kPa = c,
      R2 = summary(fit)$r.squared, n_points = nrow(g), check.names = FALSE
    )
  })
  out <- do.call(rbind, lst); rownames(out) <- NULL; out[order(out$t), ]
}

# ---- Mohr–Coulomb (p–q) - TPDSC model ----
mc_fit_model <- fit_mc_pq_by_time(pfqf_agg_model)

# ---- Mohr–Coulomb (p–q) - Experimental (from qf,ref) ----
pfqf_exp <- df_params |>
  transmute(t, s3, qf = qf_ref, pf = (qf_ref + 3*s3)/3)
mc_fit_exp <- fit_mc_pq_by_time(pfqf_exp)

# ---- Exports ----

# 0) Optional inspection: per-sample and combined (σ–ε)
if (export_inspection) {
  # combined predictions
  write.csv(pred_df,
            file.path(export_dir, "TPDSC_sigma_predictions_ALL.csv"),
            row.names = FALSE)
  
  # per (sample_id, t, s3): split and write CSVs
  splits <- split(pred_df, list(pred_df$sample_id, pred_df$t, pred_df$s3), drop = TRUE)
  for (nm in names(splits)) {
    parts <- strsplit(nm, "\\.")[[1]]
    sid <- parts[1]; tt <- parts[2]; cc <- parts[3]
    fn  <- sprintf("sigma_pred_sid%s_t%sd_s3_%skPa.csv", sid, tt, cc)
    write.csv(splits[[nm]], file.path(export_dir, fn), row.names = FALSE)
  }
}

# 1) Peaks and E50 per sample (TPDSC)
write.csv(peak_E50_tbl |>
            arrange(sample_id, t, s3),
          file.path(export_dir, "TPDSC_peak_E50_per_sample.csv"),
          row.names = FALSE)

# 2) Aggregated pf–qf (TPDSC) for MC fits
write.csv(pfqf_agg_model |>
            arrange(t, s3),
          file.path(export_dir, "TPDSC_peak_pf_qf_aggregated.csv"),
          row.names = FALSE)

# 3) Mohr–Coulomb (p–q) summaries - model and experimental
write.csv(mc_fit_model,
          file.path(export_dir, "TPDSC_MC_fit_summary_MODEL_pq.csv"),
          row.names = FALSE)

write.csv(mc_fit_exp,
          file.path(export_dir, "TPDSC_MC_fit_summary_EXPERIMENTAL_pq.csv"),
          row.names = FALSE)

# ---- Console summary ----
cat("\n--- Peak & E50 (per sample, TPDSC) ---\n")
print(peak_E50_tbl |>
        arrange(sample_id, t, s3))

cat("\n--- Aggregated pf–qf (TPDSC, for MC) ---\n")
print(pfqf_agg_model |>
        arrange(t, s3))

cat("\n--- Mohr–Coulomb (p–q) — TPDSC model ---\n")
print(mc_fit_model)

cat("\n--- Mohr–Coulomb (p–q) — EXPERIMENTAL ---\n")
print(mc_fit_exp)