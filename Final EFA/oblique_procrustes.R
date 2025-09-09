# ===================================================================
# ENHANCED PROCRUSTES EFA WITH MAP-BASED FACTOR SELECTION
# ===================================================================

# Load necessary libraries
library(dplyr)
library(psych)
library(tidyr)
library(mice)
library(tidyverse)

cat("=== ITEM REMOVAL ANALYSIS ===\n")

# Analyze which items were consistently removed across imputations
removal_log <- read.csv("/Users/march341/Desktop/MIDUS/Final_EFA/EFA_removal_log_all_imputations_Velicer_MAP_Final.csv")

percent_removed <- removal_log %>%
  distinct(imputation, item_removed) %>%
  count(item_removed, name = "n_imputations") %>%
  mutate(total_imps = n_distinct(removal_log$imputation),
         percent = 100 * n_imputations / total_imps)

cat("Items removed across imputations (≥50% removal threshold):\n")
print(percent_removed)

# Function to analyze factor suggestions 
analyze_factor_suggestions <- function(csv_file_path) {
  data <- readr::read_csv(csv_file_path, show_col_types = FALSE)
  
  final_suggestions <- data %>%
    group_by(imputation) %>%
    slice_max(iteration, n = 1) %>%
    select(imputation, suggested_factors) %>%
    ungroup()
  
  results <- list(
    individual_suggestions = final_suggestions$suggested_factors,
    mean_factors = mean(final_suggestions$suggested_factors),
    median_factors = median(final_suggestions$suggested_factors),
    sd_factors = sd(final_suggestions$suggested_factors),
    min_factors = min(final_suggestions$suggested_factors),
    max_factors = max(final_suggestions$suggested_factors),
    mode_factors = as.numeric(names(sort(table(final_suggestions$suggested_factors), decreasing = TRUE)[1])),
    n_imputations = nrow(final_suggestions),
    final_suggestions_table = table(final_suggestions$suggested_factors)
  )
  
  return(results)
}

# Analyze factor suggestions from iterative process
factor_suggestions <- analyze_factor_suggestions("/Users/march341/Desktop/MIDUS/Final_EFA/EFA_removal_log_all_imputations_Velicer_MAP_Final.csv")

cat("\nFactor suggestions from iterative removal process:\n")
print(factor_suggestions)

cat("\n=== DATA PREPARATION ===\n")

data_raw <- read.csv("/Users/march341/Desktop/MIDUS/completed_data_for_analysis_final.csv")

data <- data_raw %>%
  dplyr::select(-c(ID, Sex, .id, Income, Age, Education, Sexual_Orientation, Status_SAQ, Race, Generalized_Anxiety, Emotional_Problems_Drinking, Desire_Drinking, Month_or_more_Drink, Need_Drink_Same, Alcohol_Affected_Work, Anhedonia, Depressed_Affect, Alcohol_More_Than_Intended, Approach_Health, Asking_Help_Naturally, Avoid_Distraction, Aware_Body, Bad_Happens_Prevented, Careless, Choose_Goals, Dissapointment_Low_Goal, Dont_Like_Help, Exciting_Look_Forward, Hate_Temp, Headache_Frequency, Intercourse_Pain, Leaking_Urine, Loud_Noise, Low_Pain_Tolerance, Lower_Expectations, Mental_Emotional, More_Or_Less_Positive, Nothing_Cheer, People_Take_Advantage, Positive_About_Self, Quick_Sense_Hunger,Releived_No_Responsibility, Should_Obey_Law, Some_Wander_Aimlessly, Sympathy_Limits, To_Reach_Goals, Trust_Friends, When_Difficult, Without_skill, Worry_Compare, Live_One_Day, Asking_Help_Naturally, Enjoy_Plans_Future, Effective_Talking_People,  Good_Influence_People, Others_Turn_To_Me, Life_Day_By_Day, Irritabiltiy_Frequency,Irritable_Freq_30_Day, Rarely_Give_Up, Frustrated_Freq, Feel_Close_Others,  Felt_Close_Others_Freq,  Enthusiastic_Freq, Full_Of_Life_Freq, Extreme_Happy_Freq, Good_Spirits_Freq, Cheerful_Freq, Important_Help,           Not_Happy_If_Friend_Trouble, Moved_By_Other_Hardship, Important_Sympathetic, Moody, Worrying, Responsible, Nervous, Helpful_Set_Future_Goal)) #Dissapointment_Low_Goal, Asking_Help_Naturally, Enjoy_Plans_Future, Live_One_Day, responsible, Life_Day_By_Day, Irritabiltiy_Frequency, Irritable_Freq_30_Day, Rarely_Give_Up, Frustrated_Freq, and Feel_Close_Others were removed to reduce doublets. Felt_Close_Others_Freq, Felt_Belong_Freq, Enthusiastic_Freq, Full_Of_Life_Freq, Extreme_Happy_Freq, Good_Spirits_Freq, and Cheerful_Freq were removed because they are not within the purview of HiTOP.  Effective_Talking_People, Good_Influence_People, and Others_Turn_To_Me were removed because they are just assertiveness and do not group with any other factors. Important_Help, Not_Happy_If_Friend_Trouble, Moved_By_Other_Hardship, and Important_Sympathetic were removed because they are outside of the purview of HiTOP. Like_Time_Alone, Seek_Friends, I_Am_Warm_Person, and Prefer_No_Others were removed because they are just introversion, they are not pathological, and don't group with other factors. Overwhelmed_Responsibility, Not_Afraid_Opinion, Confidence_If_Contrary, Difficult_Voice_Opinion, Influenced_By_People, Most_See_Love, Others_Describe_Giving were removed because they do not group with anything and they are outside of the purview of HiTOP.

#Pull out ordinal items
ordinal_items <- setdiff(names(data), ".imp")
cat("Number of items for analysis:", length(ordinal_items), "\n")

data_clean <- data %>%
  mutate(across(any_of(ordinal_items),
                ~ as.integer(round(as.numeric(.))))) %>%
  mutate(across(any_of(ordinal_items),
                ~ ordered(., levels = sort(unique(.))))) %>%
  dplyr::select(any_of(c(".imp", ordinal_items))) %>%
  mutate(.imp = as.integer(.imp))

n_total <- max(table(data_clean$.imp))

cat("Data cleaning complete.\n")

cat("\n=== MAP-BASED FACTOR SELECTION ===\n")

get_map_suggestions <- function(imputed_data) {
  imputation_ids <- unique(imputed_data$.imp)
  map_results <- data.frame(imputation = integer(), map_suggested_factors = integer())
  
  for (i in imputation_ids) {
    message("Computing MAP for imputation ", i)
    
    this_data <- imputed_data %>%
      filter(.imp == i) %>%
      select(-.imp)
    
    pc <- psych::polychoric(as.data.frame(this_data), correct = 0)
    R <- pc$rho
    n_obs <- nrow(this_data)
    
    vss_out <- VSS(R, n.obs = n_obs, n = min(30, ncol(this_data)), rotate = "none", plot = FALSE)
    map_suggestion <- which.min(vss_out$map)
    
    map_results <- rbind(map_results, 
                         data.frame(imputation = i, map_suggested_factors = map_suggestion))
  }
  
  return(map_results)
}

analyze_map_results <- function(map_results) {
  summary_stats <- list(
    individual_suggestions = map_results$map_suggested_factors,
    mean_factors = mean(map_results$map_suggested_factors),
    median_factors = median(map_results$map_suggested_factors),
    sd_factors = sd(map_results$map_suggested_factors),
    min_factors = min(map_results$map_suggested_factors),
    max_factors = max(map_results$map_suggested_factors),
    mode_factors = as.numeric(names(sort(table(map_results$map_suggested_factors), decreasing = TRUE)[1])),
    n_imputations = nrow(map_results),
    frequency_table = table(map_results$map_suggested_factors)
  )
  return(summary_stats)
}

# Run your MAP analysis
map_results <- get_map_suggestions(data_clean)
map_analysis <- analyze_map_results(map_results)

cat("MAP analysis results:\n")
print(map_analysis)

## =========================
## OBLIQUE PROCRUSTES PIPELINE (1st + 2nd order)
## =========================

set.seed(123)

## =========================
## Utilities
## =========================

# Ensure row/col names exist
ensure_names <- function(L, item_names = NULL) {
  if (is.null(rownames(L))) {
    if (!is.null(item_names)) rownames(L) <- item_names else rownames(L) <- paste0("Item_", seq_len(nrow(L)))
  }
  if (is.null(colnames(L))) colnames(L) <- paste0("F", seq_len(ncol(L)))
  L
}

# Print + save loadings BEFORE any cutoff
print_loadings_pre_cutoff <- function(L, file_prefix = "first_order",
                                      digits = 2, order_rows = TRUE) {
  stopifnot(is.matrix(L))
  L <- ensure_names(L)  # make sure row/col names exist
  
  if (order_rows) {
    primary <- apply(abs(L), 1, which.max)
    mag     <- apply(abs(L), 1, max)
    ord     <- order(primary, -mag)
    L_ord   <- L[ord, , drop = FALSE]
  } else {
    L_ord <- L
  }
  
  Rnum <- round(L_ord, digits)
  pretty_tab <- cbind(Item = rownames(L_ord), as.data.frame(Rnum, check.names = FALSE))
  
  cat("\n=== LOADINGS (pre-cutoff) —", file_prefix, "===\n")
  print(pretty_tab, row.names = FALSE)
  
  write.csv(pretty_tab, paste0(file_prefix, "_loadings_pre_cutoff_full.csv"), row.names = FALSE)
  write.csv(cbind(Item = rownames(L_ord), as.data.frame(L_ord, check.names = FALSE)),
            paste0(file_prefix, "_loadings_pre_cutoff_unrounded.csv"),
            row.names = FALSE)
}

# Oblique Procrustes (unconstrained linear transform)
oblique_procrustes <- function(X, Y, Phi = NULL, ridge = 1e-6, standardize_phi = TRUE) {
  stopifnot(ncol(X) == ncol(Y))
  XtX <- t(X) %*% X
  if (ridge > 0) XtX <- XtX + diag(ridge, ncol(X))
  T <- solve(XtX, t(X) %*% Y)
  
  if (is.null(Phi)) return(list(T = T, Phi_star = NULL))
  
  Phi_star <- t(T) %*% Phi %*% T
  if (standardize_phi) {
    D <- diag(1 / sqrt(pmax(diag(Phi_star), .Machine$double.eps)))
    T <- T %*% D
    Phi_star <- t(D) %*% Phi_star %*% D
  }
  list(T = T, Phi_star = Phi_star)
}

# Align a set of loading matrices (and Φ’s) to a target via oblique Procrustes.
align_to_target_oblique <- function(loadings_list, phis_list, target, ridge = 1e-6, standardize_phi = TRUE) {
  stopifnot(length(loadings_list) == length(phis_list))
  aligned_L <- vector("list", length(loadings_list))
  rotated_Phi <- vector("list", length(phis_list))
  T_list <- vector("list", length(loadings_list))
  
  for (i in seq_along(loadings_list)) {
    X <- loadings_list[[i]]
    Y <- target
    Phi_i <- phis_list[[i]]
    if (!all(dim(X) == dim(Y))) stop("Loading matrices must have same dimensions")
    
    fit <- oblique_procrustes(X, Y, Phi = Phi_i, ridge = ridge, standardize_phi = standardize_phi)
    T_i <- fit$T
    aligned_L[[i]] <- X %*% T_i
    rotated_Phi[[i]] <- fit$Phi_star
    T_list[[i]] <- T_i
    
    dimnames(aligned_L[[i]]) <- dimnames(X)
    dimnames(rotated_Phi[[i]]) <- list(colnames(X), colnames(X))
  }
  list(aligned = aligned_L, rotated_phis = rotated_Phi, Ts = T_list)
}

iterative_procrustes_oblique <- function(loadings_list, phis_list, max_iter = 25, tol = 1e-7,
                                         ridge = 1e-6, standardize_phi = TRUE, verbose = TRUE) {
  if (verbose) cat("Starting iterative oblique Procrustes alignment...\n")
  target <- Reduce("+", loadings_list) / length(loadings_list)
  aligned <- loadings_list
  rotated_phis <- phis_list
  Ts <- replicate(length(loadings_list), diag(ncol(loadings_list[[1]])), simplify = FALSE)
  
  for (iter in 1:max_iter) {
    out <- align_to_target_oblique(loadings_list, phis_list, target, ridge, standardize_phi)
    aligned_new <- out$aligned
    new_target <- Reduce("+", aligned_new) / length(aligned_new)
    change <- sqrt(sum((new_target - target)^2))
    if (verbose) cat("Iteration", iter, "- Change:", round(change, 6), "\n")
    
    target <- new_target
    aligned <- aligned_new
    rotated_phis <- out$rotated_phis
    Ts <- out$Ts
    if (change < tol) {
      if (verbose) cat("Converged after", iter, "iterations\n")
      break
    }
  }
  list(final_target = target,
       aligned_loadings = aligned,
       rotated_phis = rotated_phis,
       Ts = Ts,
       n_iterations = iter)
}


# Correct Tucker’s congruence (per factor/column)
tucker_congruence_vec <- function(A, B) {
  stopifnot(all(dim(A) == dim(B)))
  sapply(seq_len(ncol(A)), function(j) {
    a <- A[, j]; b <- B[, j]
    as.numeric(sum(a * b) / sqrt(sum(a^2) * sum(b^2)))
  })
}
assess_alignment <- function(aligned_list, target) {
  tc_mat <- sapply(aligned_list, function(Li) tucker_congruence_vec(Li, target))
  rownames(tc_mat) <- colnames(target)
  list(tucker_by_matrix = tc_mat,
       mean_congruence = rowMeans(tc_mat))
}

## =========================
## FIRST-ORDER: EFA per imputation
## =========================
cat("\n=== STEP 1: FIRST-ORDER EFA PER IMPUTATION ===\n")
imputation_ids <- sort(unique(data_clean$.imp))
all_loadings_raw <- list()
all_phis <- list()

for (imp in imputation_ids) {
  cat("Processing imputation", imp, "\n")
  temp <- data_clean %>%
    dplyr::filter(.imp == imp) %>%
    dplyr::select(all_of(ordinal_items))
  
  # polychoric correlations
  cor_result <- psych::polychoric(as.data.frame(temp), correct = 0)
  cor_matrix <- cor_result$rho
  
  efa_result <- psych::fa(cor_matrix, nfactors = 10, fm = "minres", rotate = "oblimin")
  L <- as.matrix(efa_result$loadings)  # pattern loadings
  dimnames(L) <- dimnames(efa_result$loadings)
  
  all_loadings_raw[[as.character(imp)]] <- L
  all_phis[[as.character(imp)]] <- efa_result$Phi
}

## =========================
## FIRST-ORDER: Oblique Procrustes across imputations
## =========================
cat("\n=== STEP 2: OBLIQUE PROCRUSTES ALIGNMENT (1st order) ===\n")
iter_res <- iterative_procrustes_oblique(all_loadings_raw, all_phis, max_iter = 50, tol = 1e-7)

pooled_loadings_oblique <- iter_res$final_target
rotated_phis <- iter_res$rotated_phis
pooled_phi <- Reduce("+", rotated_phis) / length(rotated_phis)
pooled_phi <- symm_unit(pooled_phi)

## Print + save PRE-CUTOFF 1st-order loadings
L <- pooled_loadings_oblique
L <- ensure_names(L, item_names = ordinal_items)
print_loadings_pre_cutoff(L, file_prefix = "first_order")

## =========================
## FIRST-ORDER: Communalities (oblique) + assignment
## =========================
cat("\n=== STEP 3: COMMUNALITIES + ASSIGNMENT (1st order) ===\n")
communalities <- diag(L %*% pooled_phi %*% t(L))

which_fac <- apply(abs(L), 1, which.max)
max_abs   <- apply(abs(L), 1, max)
fac_name  <- colnames(L)[which_fac]

primary_factor <- ifelse(max_abs >= 0.30 & communalities >= 0.25, fac_name, NA_character_)
keep <- !is.na(primary_factor)
items_by_factor <- split(rownames(L)[keep], primary_factor[keep])
unassigned <- rownames(L)[!keep]

cat("First-order factor structure:\n")
for (factor_name in names(items_by_factor)) {
  cat(" ", factor_name, "(", length(items_by_factor[[factor_name]]), "items):\n")
  cat("    ", paste(items_by_factor[[factor_name]], collapse = ", "), "\n\n")
}
cat("Unassigned items (", length(unassigned), "): ", paste(unassigned, collapse = ", "), "\n", sep = "")

# Save assignment lists
write.csv(
  data.frame(Item = unlist(items_by_factor, use.names = FALSE),
             AssignedFactor = rep(names(items_by_factor), lengths(items_by_factor))),
  "first_order_item_assignments.csv",
  row.names = FALSE
)
write.csv(data.frame(Unassigned = unassigned), "first_order_unassigned_items.csv", row.names = FALSE)


## =========================
## Exports (final pooled objects)
## =========================
cat("\n=== STEP 7: EXPORTS ===\n")
write.csv(L,                          "first_order_pooled_loadings.csv")
write.csv(pooled_phi,                 "first_order_factor_correlations_phi.csv")
write.csv(L %*% pooled_phi %*% t(L),  "first_order_reproduced_cov_from_pattern.csv")


cat("\nAll done. Files written to working directory.\n")
