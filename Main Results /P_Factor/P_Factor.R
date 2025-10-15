setwd("/users/2/march341/multiple_imputation")

# Packages ---------------------------------------------------------------
packages <- c("tidyverse", "lavaan", "semTools", "lavaan.mi")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Read stacked data ------------------------------------------------------
data_raw <- read.csv("completed_data_for_analysis_final.csv")

# CFA model -----------------------------------------------
#Updated to remove items laoding below .40
cfa_model <- '
  Internalizing =~ Worry_Frequency + Nervous_Freq_30_Day + Restless_Freq_30_Day + 
         Afraid_Freq_30_Day + Frustrated_Freq_30_Day + 
         Upset_Freq_30_Day + Hopeless_Freq + Everything_Effort + 
         Worthless_Freq + Irritable_Freq + Ashamed_Freq + 
         Lonely_Frequency + Angry_Freq

  Anatagonistic =~ People_Mean + Upset_Think_Day + Minor_Setback_Irritate + 
         Mood_Up_Down + Angry_Ready_Hit + Enjoy_Hurting_Say_Mean + 
         When_Insult_Get_Even + Sometimes_Just_Hit_Someone + More_Sucessful

  Detachment =~ Others_More_Life + Disapoint_Acheivment + Gave_Up + 
         No_Good_Sense + Difficult_Arrange + 
         Maintaining_Relationships_Hard + Dont_Fit_Community + 
         Few_Close_Friends + No_Warm_Relationship + 
         Activities_Trivial + Self_Attitude

  Somaticizing =~ Backache_Frequency + Swea_Frequency + Hot_Flashes + 
         Falling_Asleep + Extremities_Aches

  Disinhibited_Externalizing =~ Actively_Carry_Plans +  Like_Difficult_Things + 
         Like_Hard_Work + High_Standards + 
         Plan_Future + Will + Change_For_Better + Too_Much_Get_Done + 
         Dont_Give_Up + Know_Want_Life + Goal_Keep_Benefits +   
         Weigh_Possibilities 
         
  P_factor =~ Internalizing + Anatagonistic + Detachment + Somaticizing + Disinhibited_Externalizing
'


# --- Pull only observed indicators (exclude second-order RHS latents) ----
pt <- lavaan::lavaanify(cfa_model, fixed.x = FALSE)
latent_names  <- unique(pt$lhs[pt$op == "=~"])         # any latent on LHS
rhs_all       <- pt$rhs[pt$op == "=~"]                  # all RHS of loadings
observed_items <- unique(rhs_all[!rhs_all %in% latent_names])  # keep observed-only

# Warn if any observed indicators aren't in file
missing_obs <- setdiff(observed_items, names(data_raw))
if (length(missing_obs)) {
  warning("Missing observed indicators in data: ", paste(missing_obs, collapse = ", "))
}

# Keep ONLY .imp + observed items (preserve model order) ------------------
cols_to_keep <- c(".imp", observed_items)
data <- data_raw %>%
  dplyr::select(dplyr::any_of(cols_to_keep)) %>%
  dplyr::select(c(intersect(".imp", names(.)), intersect(observed_items, names(.))))

# Prep ordinals (be cautious: if already integers 1..k, this is fine) -----
ordinal_items <- setdiff(names(data), ".imp")
data_num <- data %>%
  dplyr::mutate(across(dplyr::all_of(ordinal_items),
                       ~ as.integer(round(as.numeric(.))))) %>%
  dplyr::mutate(across(dplyr::all_of(ordinal_items),
                       ~ ordered(., levels = sort(unique(.))))) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(".imp"), as.integer))

data_cfa <- data_num %>%
  dplyr::select(dplyr::any_of(".imp"), dplyr::all_of(intersect(observed_items, names(data_num))))

rev_ord <- function(x) {
  #x is ordered/factor/character with numeric categories like 1..k
  xi <- as.numeric(as.character(x))
  lo <- min(xi, na.rm=TRUE); hi <- max(xi, na.rm=TRUE)
  xr <- lo + hi - xi
  ordered(xr, levels = sort(unique(xr)))
}

rev_items <- c("Angry_Freq", "More_Sucessful")
data_cfa[rev_items] <- lapply(data_cfa[rev_items], rev_ord)

# Split stacked data to list ---------------------------------------------
if (!".imp" %in% names(data_cfa)) stop(".imp column not found; stacked data expected.")
imp_list <- split(data_cfa, data_cfa$.imp)
imp_list <- lapply(imp_list, function(df) dplyr::select(df, - .imp))

# Fit CFA across ALL imputations (lavaan.mi) ------------------------------
fit_mi <- lavaan.mi::cfa.mi(
  model            = cfa_model,
  data             = imp_list,
  estimator        = "WLSMV",
  parameterization = "delta",     # match test run for stability
  std.lv           = TRUE,
  ordered          = intersect(observed_items, names(imp_list[[1]]))
)

# Outputs and summary ----------------------------------------------------
{
  sum_txt <- capture.output({
    cat("=== lavaan.mi summary (pooled across imputations) ===\n")
    print(summary(fit_mi))
    cat("\n=== Scaled/robust fit measures for WLSMV ===\n")
    print(fitMeasures(fit_mi, c(
      "chisq.scaled","df.scaled","pvalue.scaled",
      "cfi.scaled","tli.scaled","rmsea.scaled","srmr"
    )))
  })
  writeLines(sum_txt, "Pooled_P_Factor_Model_CFA_summary.txt")
}

# Parameter estimates (pooled) -------------------------------------------
pe <- parameterEstimates.mi(
  fit_mi, se = TRUE, ci = TRUE, level = 0.95,
  rsquare = TRUE, standardized = TRUE, output = "data.frame"
)
write.csv(pe, "P_Factor_parameters_standardized.csv", row.names = FALSE)

# Reliability (pooled) ---------------------------------------------------
rel <- semTools::reliability(fit_mi)
rel_df <- as.data.frame(rel)
rel_df$factor <- rownames(rel_df)
rel_df <- rel_df[, c("factor", setdiff(names(rel_df), "factor"))]
write.csv(rel_df, "_P_Factor_reliability.csv", row.names = FALSE)

# Modification Indices (pooled) ------------------------------------------
mi <- lavaan.mi::modindices.mi(fit_mi, sort. = TRUE)
write.csv(mi, "_P_Factor_modindices.csv", row.names = FALSE)

# Φ (latent covariances) and latent correlations (pooled) -----------------
# --- Build pooled Φ (covariances) and Φ as correlations from parameterEstimates.mi ----
pt <- lavaan::lavaanify(cfa_model, fixed.x = FALSE)

# Identify latent names & observed items
all_latents <- unique(pt$lhs[pt$op == "=~"])
rhs_all     <- pt$rhs[pt$op == "=~"]
observed_items <- unique(rhs_all[!rhs_all %in% all_latents])

# First-order latents = those that load on observed items
fo_latents <- unique(pt$lhs[pt$op == "=~" & !(pt$rhs %in% all_latents)])

# Pull pooled parameter table
pe_pool <- parameterEstimates.mi(
  fit_mi, se = TRUE, ci = TRUE, level = 0.95,
  standardized = TRUE, output = "data.frame"
)

# Filter latent-latent (~~) rows for first-order latents
psi_tab <- subset(pe_pool, op == "~~" & lhs %in% fo_latents & rhs %in% fo_latents)

# Build covariance matrix (using raw 'est')
lv_order <- sort(unique(c(psi_tab$lhs, psi_tab$rhs)))
Phi_cov <- matrix(NA_real_, nrow = length(lv_order), ncol = length(lv_order),
                  dimnames = list(lv_order, lv_order))
for (i in seq_len(nrow(psi_tab))) {
  a <- psi_tab$lhs[i]; b <- psi_tab$rhs[i]; v <- psi_tab$est[i]
  Phi_cov[a, b] <- v; Phi_cov[b, a] <- v
}
# Ensure diagonals filled (variances)
for (a in lv_order) {
  v <- psi_tab$est[psi_tab$lhs == a & psi_tab$rhs == a]
  if (length(v) == 1L) Phi_cov[a, a] <- v
}

# Correlation matrix: either compute from cov or use standardized ('std.all')
# (std.lv=TRUE typically fixes latent variances to 1, but we compute safely)
sd_vec <- sqrt(diag(Phi_cov))
Phi_cor <- sweep(sweep(Phi_cov, 1, sd_vec, "/"), 2, sd_vec, "/")

# Write to CSV
write.csv(Phi_cov, "P_Factor_phi_covariances.csv", row.names = TRUE)
write.csv(Phi_cor, "P_Factor_phi_correlations.csv", row.names = TRUE)

cat("\nDone. All imputations fit.\n")

# =========================================================================
# FACTOR SCORE INDETERMINACY TESTING
# =========================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("FACTOR SCORE INDETERMINACY ANALYSIS (ESSENTIAL METHODS)\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# =========================================================================
# FACTOR DETERMINACY FOR REGRESSION FACTOR SCORES
# Source: Grice (2001); Brown (2015, Chapter 5)
# =========================================================================

calculate_factor_determinacy <- function(fit_mi) {
  pe <- parameterEstimates.mi(fit_mi, standardized = TRUE, output = "data.frame")
  
  # --- Lambda (items x factors) ---
  L_df <- pe %>%
    dplyr::filter(op == "=~") %>%
    dplyr::select(item = rhs, factor = lhs, val = std.all) %>%
    tidyr::pivot_wider(names_from = factor, values_from = val, values_fill = 0)
  items   <- L_df$item
  L       <- as.matrix(L_df[,-1, drop = FALSE]); rownames(L) <- items
  factors <- colnames(L)
  
  # --- Psi (uniquenesses of items; standardized residual variances) ---
  uniq <- pe %>%
    dplyr::filter(op == "~~", lhs == rhs, lhs %in% items) %>% # Extracts item variances
    dplyr::arrange(match(lhs, items)) %>%
    dplyr::pull(std.all)
  Psi <- diag(uniq); rownames(Psi) <- items; colnames(Psi) <- items
  
  # --- Phi (latent correlations; std.all on latent-latent '~~') ---
  Phi <- diag(length(factors)); dimnames(Phi) <- list(factors, factors)
  cors <- pe %>% dplyr::filter(op == "~~", lhs %in% factors, rhs %in% factors, lhs != rhs)
  for (i in seq_len(nrow(cors))) {
    a <- cors$lhs[i]; b <- cors$rhs[i]; v <- cors$std.all[i]
    Phi[a, b] <- v; Phi[b, a] <- v
  }
  
  # --- Determinacy: p_j^2 = diag(Phi * L' * Sigma^{-1} * L * Phi) ---
  Sigma <- L %*% Phi %*% t(L) + Psi
  if (rcond(Sigma) < 1e-8) Sigma <- Sigma + diag(1e-8, nrow(Sigma))
  precision_weighted <- t(L) %*% solve(Sigma) %*% L
  p2  <- diag(Phi %*% precision_weighted %*% Phi)          # diag(Phi)=1 under std.lv=TRUE
  p   <- sqrt(p2)
  n_items <- colSums(abs(L) > 0)
  
  determinacy_results <- data.frame(
    Factor              = factors,
    N_Items             = as.integer(n_items),
    Determinacy_Precise = round(as.numeric(p), 3),
    Min_Correlation     = round(2 * p2 - 1, 3),
    Interpretation      = dplyr::case_when(
      p >= .90 ~ "Excellent",
      p >= .80 ~ "Good",
      p >= .70 ~ "Acceptable",
      TRUE     ~ "Poor - Use SEM instead"
    ),
    row.names = NULL
  )
  
  
  print(determinacy_results)
  cat("\nInterpretation (Grice, 2001):\n")
  cat("≥0.90 = Excellent (factor scores can replace latent variables)\n")
  cat("≥0.80 = Good (factor scores acceptable for most purposes)\n")
  cat("<0.80 = Problematic (substantial indeterminacy)\n")
  
  return(determinacy_results)
}

# =========================================================================
# METHOD 2: FACTOR SCORE RELIABILITY 
# Source: Hancock & Mueller (2001); Raykov (2004)
# =========================================================================

calculate_construct_reliability <- function(fit_mi) {
  pe <- parameterEstimates.mi(fit_mi, standardized = TRUE, output = "data.frame")
  
  # --- Lambda (items x factors) ---
  L_df <- pe |>
    dplyr::filter(op == "=~") |>
    dplyr::select(item = rhs, factor = lhs, val = std.all) |>
    tidyr::pivot_wider(names_from = factor, values_from = val, values_fill = 0)
  items   <- L_df$item
  L       <- as.matrix(L_df[ , -1, drop = FALSE]); rownames(L) <- items
  factors <- colnames(L)
  
  # --- Phi (latent correlations) for fallback communalities ---
  Phi <- diag(length(factors)); dimnames(Phi) <- list(factors, factors)
  cors <- pe |>
    dplyr::filter(op == "~~", lhs %in% factors, rhs %in% factors, lhs != rhs)
  for (i in seq_len(nrow(cors))) {
    a <- cors$lhs[i]; b <- cors$rhs[i]; v <- cors$std.all[i]
    Phi[a, b] <- v; Phi[b, a] <- v
  }
  
  # --- Theta (residual var/cov among items), standardized scale ---
  Theta <- matrix(0, nrow = length(items), ncol = length(items),
                  dimnames = list(items, items))
  th <- pe |> dplyr::filter(op == "~~", lhs %in% items, rhs %in% items)
  if (nrow(th) > 0) {
    for (k in seq_len(nrow(th))) {
      Theta[th$lhs[k], th$rhs[k]] <- th$std.all[k]
      Theta[th$rhs[k], th$lhs[k]] <- th$std.all[k]
    }
  }
  # Fill any missing/zero diagonals with 1 - communality (fallback)
  diag_fallback <- 1 - diag(L %*% Phi %*% t(L))
  need_fill <- which(is.na(diag(Theta)) | diag(Theta) <= 0)
  if (length(need_fill)) diag(Theta)[need_fill] <- diag_fallback[need_fill]
  
  out <- lapply(factors, function(f) {
    idx <- which(abs(L[, f]) > 0)
    lam_v <- matrix(L[idx, f], ncol = 1)
    Theta_f <- Theta[idx, idx, drop = FALSE]
    if (rcond(Theta_f) < 1e-8) Theta_f <- Theta_f + diag(1e-8, nrow(Theta_f))
    
    # Hancock & Mueller's H: H = (λ' Θ^{-1} λ) / (1 + λ' Θ^{-1} λ)
    q <- drop(t(lam_v) %*% solve(Theta_f) %*% lam_v)
    H <- q / (1 + q)
    
    # Fornell–Larcker Composite Reliability (CR, ρc):
    # ρc = (Σλ)^2 / ((Σλ)^2 + Σθ_i), using Θ for these indicators
    sum_lambda <- sum(lam_v)
    sum_theta  <- sum(diag(Theta_f))
    rho_c <- (sum_lambda^2) / (sum_lambda^2 + sum_theta)
    
    data.frame(
      Factor = f,
      Hancock_H = round(as.numeric(H), 3),
      Fornell_Larcker = round(as.numeric(rho_c), 3),
      H_Interpretation = dplyr::case_when(
        H >= .90 ~ "Excellent",
        H >= .80 ~ "Good",
        H >= .70 ~ "Acceptable",
        TRUE     ~ "Poor"
      ),
      FL_Interpretation = dplyr::case_when(   # CR (ρc) conventions
        rho_c >= .90 ~ "Excellent",
        rho_c >= .80 ~ "Good",
        rho_c >= .70 ~ "Acceptable",
        TRUE         ~ "Poor"
      ),
      row.names = NULL
    )
  })
  
  res <- do.call(rbind, out)
  print(res)
  cat("\nInterpretation notes:\n")
  cat("- H (MaxR(H)) ≥ 0.80 suggests reliable factor definition (Hancock & Mueller, 2001).\n")
  cat("- Fornell–Larcker ρc (composite reliability) ≥ 0.70 is commonly acceptable; ≥ 0.80 good; ≥ 0.90 excellent.\n")
  cat("- These are NOT the same as factor-score determinacy (report p and 2p^2-1 separately).\n")
  res
}

# =========================================================================
# METHOD 3: FACTOR SCORE UNCERTAINTY (ESSENTIAL FOR INTERPRETATION)
# Source: Beauducel (2011); Rodriguez et al. (2016)
# =========================================================================

calculate_practical_uncertainty <- function(fit_mi) {
  pe <- parameterEstimates.mi(fit_mi, standardized = TRUE, output = "data.frame")
  
  # --- Lambda (items x factors) ---
  L_df <- pe |>
    dplyr::filter(op == "=~") |>
    dplyr::select(item = rhs, factor = lhs, val = std.all) |>
    tidyr::pivot_wider(names_from = factor, values_from = val, values_fill = 0)
  items   <- L_df$item
  L       <- as.matrix(L_df[, -1, drop = FALSE]); rownames(L) <- items
  factors <- colnames(L)
  
  # --- Phi (latent correlations; standardized) ---
  Phi <- diag(length(factors)); dimnames(Phi) <- list(factors, factors)
  cors <- pe |> dplyr::filter(op == "~~", lhs %in% factors, rhs %in% factors, lhs != rhs)
  for (i in seq_len(nrow(cors))) {
    a <- cors$lhs[i]; b <- cors$rhs[i]; v <- cors$std.all[i]
    Phi[a, b] <- v; Phi[b, a] <- v
  }
  diag(Phi) <- 1
  
  # --- Theta (residual var/cov among items), standardized scale ---
  Theta <- matrix(0, nrow = length(items), ncol = length(items),
                  dimnames = list(items, items))
  th <- pe |> dplyr::filter(op == "~~", lhs %in% items, rhs %in% items)
  if (nrow(th) > 0) {
    for (k in seq_len(nrow(th))) {
      Theta[th$lhs[k], th$rhs[k]] <- th$std.all[k]
      Theta[th$rhs[k], th$lhs[k]] <- th$std.all[k]
    }
  }
  # Fill missing/zero diagonals with 1 - communality fallback
  diag_fallback <- 1 - diag(L %*% Phi %*% t(L))
  d <- diag(Theta); idx <- which(is.na(d) | d <= 0)
  if (length(idx)) d[idx] <- diag_fallback[idx]
  diag(Theta) <- d
  
  # --- Determinacy p from Sigma = L Phi L' + Theta ---
  Sigma <- L %*% Phi %*% t(L) + Theta
  if (rcond(Sigma) < 1e-8) Sigma <- Sigma + diag(1e-8, nrow(Sigma))
  mid <- t(L) %*% solve(Sigma) %*% L
  p2  <- diag(Phi %*% mid %*% Phi)                 # std.lv -> diag(Phi)=1
  p   <- sqrt(p2)
  names(p) <- factors
  
  # --- Uncertainty metrics from p ---
  sem_fs      <- sqrt(1 - p^2)                     # SEE for standardized true scores
  ci_width_sd <- 2 * 1.96 * sem_fs                 # full 95% width (in SD units)
  
  # Exact median misclassification: P(sign(score) != sign(true))
  p_median <- 0.5 - asin(p) / pi                   # bivariate normal identity
  
  # Extreme-group error (top tercile FDR): P(true <= t | score > t)
  tthr <- qnorm(2/3)                               # ~0.4307
  Phi2 <- function(a, b, rho) {
    if (requireNamespace("pbivnorm", quietly = TRUE)) {
      pbivnorm::pbivnorm(a, b, rho)
    } else if (requireNamespace("mvtnorm", quietly = TRUE)) {
      mvtnorm::pmvnorm(upper = c(a, b),
                       corr = matrix(c(1, rho, rho, 1), 2, 2))[1]
    } else {
      NA_real_  # if neither package available
    }
  }
  phi_t  <- pnorm(tthr)
  bvn_tt <- vapply(as.numeric(p), function(rho) Phi2(tthr, tthr, rho), numeric(1))
  p_ext  <- (phi_t - bvn_tt) / (1 - phi_t)         # symmetric for bottom tercile
  
  # --- Assemble with your interpretations ---
  out <- data.frame(
    Factor              = factors,
    Determinacy         = round(p, 3),
    SEM_FactorScore     = round(sem_fs, 3),
    CI95_Width_SD       = round(ci_width_sd, 2),
    P_Misclass_Median   = round(p_median, 3),
    P_Misclass_Extreme  = round(p_ext, 3),
    Classification_Reliability = dplyr::case_when(
      p_median < 0.10 ~ "Excellent",
      p_median < 0.20 ~ "Good",
      p_median < 0.30 ~ "Fair",
      TRUE            ~ "Poor"
    ),
    stringsAsFactors = FALSE
  )
  
  print(out)
  cat("\nInterpretation:\n")
  cat("- P_Misclass_Median: exact bivariate-normal error for a median split (top vs bottom half).\n")
  cat("- P_Misclass_Extreme: P(true ≤ t | score > t) with t = Φ^{-1}(2/3)≈0.431 (top-tercile false discovery rate; bottom is symmetric).\n")
  cat("- SEM_FactorScore = sqrt(1 - p^2); 95% width ≈ 3.92*SEM (both assume standardized scale).\n")
  out
}


# =========================================================================
# INTEGRATED DECISION MATRIX
# Source: Reise & Waller (2009); Beauducel (2011)
# =========================================================================

create_decision_matrix <- function(determinacy, reliability, uncertainty) {
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("FACTOR SCORE USAGE DECISION MATRIX\n")
  cat("(Based on Reise & Waller, 2009; Beauducel, 2011)\n")
  cat(paste(rep("=", 60), collapse=""), "\n\n")
  
  decision_matrix <- determinacy %>%
    select(Factor, Determinacy = Determinacy_Precise, Min_Correlation) %>%
    left_join(reliability[, c("Factor", "Hancock_H")], by = "Factor") %>%
    left_join(uncertainty[, c("Factor", "P_Misclass_Median")], by = "Factor")
  
  # Create decision rules based on literature
  decision_matrix$Use_Factor_Scores <- apply(decision_matrix, 1, function(x) {
    det <- as.numeric(x["Determinacy"])
    h <- as.numeric(x["Hancock_H"])
    misclass <- as.numeric(x["P_Misclass_Median"])
    
    if (det >= 0.90 && h >= 0.90 && misclass < 0.10) {
      "YES - Highly Recommended"
    } else if (det >= 0.80 && h >= 0.80 && misclass < 0.20) {
      "YES - Acceptable"
    } else if (det >= 0.70 && misclass < 0.30) {
      "CAUTION - Use with limitations noted"
    } else {
      "NO - Use latent variable modeling"
    }
  })
  
  decision_matrix$Recommended_Use <- apply(decision_matrix, 1, function(x) {
    det <- as.numeric(x["Determinacy"])
    
    if (det >= 0.90) {
      "Individual assessment, group comparisons, regression"
    } else if (det >= 0.80) {
      "Group comparisons, correlation/regression"
    } else if (det >= 0.70) {
      "Group-level analysis only"
    } else {
      "Avoid factor scores - use SEM/latent variables"
    }
  })
  
  print(decision_matrix)
  
  return(decision_matrix)
}

# =========================================================================
# RUN ANALYSIS
# =========================================================================

# Run essential analyses
determinacy_results <- calculate_factor_determinacy(fit_mi)
reliability_results <- calculate_construct_reliability(fit_mi)
uncertainty_results <- calculate_practical_uncertainty(fit_mi)

# Create decision matrix
decision_matrix <- create_decision_matrix(
  determinacy_results, 
  reliability_results, 
  uncertainty_results
)

# Save results
write.csv(decision_matrix, "_P_Factor_decision_matrix.csv", row.names = FALSE)

# =========================================================================
# FINAL RECOMMENDATIONS
# =========================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RECOMMENDATIONS BASED ON LITERATURE\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

cat("DECISION CRITERIA (from psychometric literature):\n\n")

cat("1. INDIVIDUAL ASSESSMENT (e.g., clinical diagnosis):\n")
cat("   - Require determinacy ≥ 0.90 (Grice, 2001)\n")
cat("   - Hancock's H ≥ 0.90\n")
cat("   - P(misclassification) < 0.10\n\n")

cat("2. GROUP COMPARISONS (e.g., t-tests, ANOVA):\n")
cat("   - Require determinacy ≥ 0.80 (Brown, 2015)\n")
cat("   - Hancock's H ≥ 0.80\n")
cat("   - P(misclassification) < 0.20\n\n")

cat("3. CORRELATION/REGRESSION STUDIES:\n")
cat("   - Minimum determinacy ≥ 0.70 (Beauducel, 2011)\n")
cat("   - Note: Correlations will be attenuated\n\n")

cat("4. IF CRITERIA NOT MET:\n")
cat("   - Use SEM with latent variables\n")
cat("   - Or use simple sum/mean scores (may be more robust)\n")
cat("   - Report limitations clearly\n\n")

cat(paste(rep("=", 60), collapse=""), "\n")

