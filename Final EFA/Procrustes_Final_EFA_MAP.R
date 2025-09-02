# ===================================================================
# ENHANCED PROCRUSTES EFA WITH MAP-BASED FACTOR SELECTION
# ===================================================================

# Load necessary libraries
library(dplyr)
library(psych)
library(tidyr)
library(mice)

cat("=== STEP 1: ITEM REMOVAL ANALYSIS ===\n")
# Analyze which items were consistently removed across imputations

removal_log <- read.csv("/Users/march341/Desktop/Final EFA/EFA_removal_log_all_imputations_Velicer_MAP_Final.csv")

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
factor_suggestions <- analyze_factor_suggestions("/Users/march341/Desktop/Final EFA/EFA_removal_log_all_imputations_Velicer_MAP_Final.csv")

cat("\nFactor suggestions from iterative removal process:\n")
print(factor_suggestions)

cat("\n=== STEP 2: DATA PREPARATION ===\n")

data_raw <- read.csv("/Users/march341/Desktop/MIDUS/completed_data_for_analysis_final.csv")

data <- data_raw %>%
  dplyr::select(-c(ID, Sex, .id, Income, Age, Education, Sexual_Orientation, Status_SAQ, Race, Generalized_Anxiety, Emotional_Problems_Drinking, Desire_Drinking, Month_or_more_Drink, Need_Drink_Same, Alcohol_Affected_Work, Anhedonia, Depressed_Affect, Alcohol_More_Than_Intended, Approach_Health, Asking_Help_Naturally, Avoid_Distraction, Aware_Body, Bad_Happens_Prevented, Careless, Choose_Goals, Dissapointment_Low_Goal, Dont_Like_Help, Exciting_Look_Forward, Hate_Temp, Headache_Frequency, Intercourse_Pain, Leaking_Urine, Loud_Noise, Low_Pain_Tolerance, Lower_Expectations, Mental_Emotional, More_Or_Less_Positive, Nothing_Cheer, People_Take_Advantage, Positive_About_Self, Quick_Sense_Hunger,Releived_No_Responsibility, Should_Obey_Law, Some_Wander_Aimlessly, Sympathy_Limits, To_Reach_Goals, Trust_Friends, When_Difficult, Without_skill, Worry_Compare, Live_One_Day, Asking_Help_Naturally, Enjoy_Plans_Future, Effective_Talking_People,  Good_Influence_People, Others_Turn_To_Me, Life_Day_By_Day, Irritabiltiy_Frequency,Irritable_Freq_30_Day, Rarely_Give_Up, Frustrated_Freq, Feel_Close_Others,  Felt_Close_Others_Freq,  Enthusiastic_Freq, Full_Of_Life_Freq, Extreme_Happy_Freq, Good_Spirits_Freq, Cheerful_Freq, Important_Help,           Not_Happy_If_Friend_Trouble, Moved_By_Other_Hardship, Important_Sympathetic, Moody, Worrying, Responsible, Nervous, Helpful_Set_Future_Goal)) #Dissapointment_Low_Goal, Asking_Help_Naturally, Enjoy_Plans_Future, Live_One_Day, responsible, Life_Day_By_Day, Irritabiltiy_Frequency, Irritable_Freq_30_Day, Rarely_Give_Up, Frustrated_Freq, and Feel_Close_Others were removed to reduce doublets. Felt_Close_Others_Freq, Felt_Belong_Freq, Enthusiastic_Freq, Full_Of_Life_Freq, Extreme_Happy_Freq, Good_Spirits_Freq, and Cheerful_Freq were removed because they are not within the purview of HiTOP.  Effective_Talking_People, Good_Influence_People, and Others_Turn_To_Me were removed because they are just assertiveness and do not group with any other factors. Important_Help, Not_Happy_If_Friend_Trouble, Moved_By_Other_Hardship, and Important_Sympathetic were removed because they are outside of the purview of HiTOP. Like_Time_Alone, Seek_Friends, I_Am_Warm_Person, and Prefer_No_Others were removed because they are just introversion, they are not pathological, and don't group with other factors. Overwhelmed_Responsibility, Not_Afraid_Opinion, Confidence_If_Contrary, Difficult_Voice_Opinion, Influenced_By_People, Most_See_Love, Others_Describe_Giving were removed because they do not group with anything and they are outside of the purview of HiTOP.

ordinal_items <- setdiff(names(data), ".imp")
cat("Number of items for analysis:", length(ordinal_items), "\n")

data_clean <- data %>%
  mutate(across(any_of(ordinal_items),
                ~ as.integer(round(as.numeric(.))))) %>%
  mutate(across(any_of(ordinal_items),
                ~ ordered(., levels = sort(unique(.))))) %>%
  dplyr::select(any_of(c(".imp", ordinal_items))) %>%
  mutate(.imp = as.integer(.imp))

cat("Data cleaning complete.\n")

cat("\n=== STEP 3: MAP-BASED FACTOR SELECTION ===\n")

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

cat("\n=== STEP 4: PROCRUSTES ROTATION IMPLEMENTATION ===\n")

procrustes_rotation <- function(X, Y) {
  # Find optimal orthogonal rotation T such that XT ≈ Y
  svd_result <- svd(t(Y) %*% X) #Singluar value decomposition
  T <- svd_result$v %*% t(svd_result$u)
  return(T)
}

align_to_target <- function(loadings_list, target) {
  aligned_list <- list()
  
  for (i in seq_along(loadings_list)) {
    current_loadings <- loadings_list[[i]]
    
    if (!all(dim(current_loadings) == dim(target))) {
      stop("Loading matrices must have same dimensions")
    }
    
    T <- procrustes_rotation(current_loadings, target)
    aligned_loadings <- current_loadings %*% T
    dimnames(aligned_loadings) <- dimnames(current_loadings)
    aligned_list[[i]] <- aligned_loadings
  }
  
  return(aligned_list)
}

iterative_procrustes <- function(loadings_list, max_iter = 10, tol = 1e-6) {
  cat("Starting iterative Procrustes alignment...\n")
  target <- Reduce("+", loadings_list) / length(loadings_list)
  
  for (iter in 1:max_iter) {
    old_target <- target
    aligned <- align_to_target(loadings_list, target)
    new_target <- Reduce("+", aligned) / length(aligned)
    
    change <- sqrt(sum((new_target - old_target)^2))
    cat("Iteration", iter, "- Change:", round(change, 6), "\n")
    
    if (change < tol) {
      cat("Converged after", iter, "iterations\n")
      break
    }
    target <- new_target
  }
  
  return(list(final_target = target, aligned_loadings = aligned, n_iterations = iter))
}

cat("\n=== STEP 5: ENHANCED FIRST-ORDER EFA WITH PROCRUSTES ===\n")

imputation_ids <- unique(data_clean$.imp)
all_loadings_raw <- list()
all_communalities <- list()
all_phis <- list()

# Run EFA on each imputation 
n_factors_to_use <- 11
cat("Running EFA with", n_factors_to_use, "factors on each imputation...\n")

for (imp in imputation_ids) {
  cat("Processing imputation", imp, "\n")
  
  temp <- data_clean %>%
    filter(.imp == imp) %>%
    dplyr::select(-.imp)
  
  cor_result <- psych::polychoric(as.data.frame(temp[, ordinal_items, drop = FALSE]), correct = 0)
  cor_matrix <- cor_result$rho
  
  # Use MAP-determined factors 
  efa_result <- fa(cor_matrix, nfactors = n_factors_to_use, fm = "minres", rotate = "oblimin")
  
  all_loadings_raw[[imp]] <- as.matrix(efa_result$loadings)
  all_communalities[[imp]] <- setNames(efa_result$communality, rownames(efa_result$loadings))
  all_phis[[imp]] <- efa_result$Phi
}

# Apply Procrustes alignment
cat("Applying Procrustes alignment...\n")
iterative_result <- iterative_procrustes(all_loadings_raw)
pooled_loadings_procrustes <- iterative_result$final_target

# Pool other results 
pooled_communality <- Reduce("+", all_communalities) / length(all_communalities)
pooled_phi <- Reduce("+", all_phis) / length(all_phis)

cat("First-order analysis complete.\n")

cat("\n=== STEP 7: FACTOR ASSIGNMENT  ===\n")

L <- as.matrix(pooled_loadings_procrustes)  # Use Procrustes-aligned loadings
rownames(L) <- colnames(data_clean)[!colnames(data_clean) %in% c(".imp")]

which_fac <- apply(abs(L), 1, which.max)
max_abs <- apply(abs(L), 1, max)
fac_name <- colnames(L)[which_fac]
communalities <- rowSums(L^2)

# Selection Criteria: |loading| ≥ .30 AND communality ≥ .25
primary_factor <- ifelse(max_abs >= 0.30 & communalities >= 0.25, fac_name, NA_character_)

keep <- !is.na(primary_factor)
items_by_factor <- split(rownames(L)[keep], primary_factor[keep])

cat("First-order factor structure:\n")
for(factor_name in names(items_by_factor)) {
  cat("", factor_name, "(", length(items_by_factor[[factor_name]]), "items):\n")
  cat("   ", paste(items_by_factor[[factor_name]], collapse = ", "), "\n\n")
}

unassigned <- rownames(L)[!keep]
cat("Unassigned items (", length(unassigned), "):", paste(unassigned, collapse = ", "), "\n")

cat("\n=== STEP 8: HIERARCHICAL ANALYSIS  ===\n")
# Apply hierarchical analysis with Procrustes rotation

# MAP on pooled factor correlations
vss_out <- VSS(pooled_phi, n.obs = 7350, n = min(30, ncol(pooled_phi)), rotate = "none", plot = FALSE)
map_suggestion_hierarchical <- which.min(vss_out$map)
cat("MAP suggests", map_suggestion_hierarchical, "second-order factors\n")

# MAP only suggests one factor so check parallel analysis as well.
cat("Running parallel analysis on factor correlations...\n")
fa.parallel(pooled_phi, n.obs = 7350, fa = "fa", fm = "minres", n.iter = 10, 
            error.bars = TRUE, main = "Parallel Analysis on Pooled Factor Correlations")

# Run second-order EFA with Procrustes
all_second_order_loadings <- list()
all_second_order_phis <- list()

cat("Running second-order EFA with Procrustes alignment...\n")
for (imp in imputation_ids) {
  phi_matrix <- all_phis[[imp]]
  
  # Use parallel suggestion
  n_second_order <- 5  # or use map_suggestion_hierarchical, but that is only one factor
  
  second_order_efa <- fa(phi_matrix, nfactors = n_second_order, fm = "minres", 
                         rotate = "oblimin", scores = "tenBerge")
  
  all_second_order_loadings[[imp]] <- as.matrix(second_order_efa$loadings)
  all_second_order_phis[[imp]] <- second_order_efa$Phi
}

# Align second-order factors
second_order_iterative <- iterative_procrustes(all_second_order_loadings)
pooled_loadings_second_order <- second_order_iterative$final_target

# Pool other results 
second_order_pooled_phi <- Reduce("+", all_second_order_phis) / length(all_second_order_phis)


cat("\n=== STEP 9: SECOND-ORDER FACTOR ASSIGNMENT ===\n")

L2 <- as.matrix(pooled_loadings_second_order)
if(is.null(rownames(L2))) {
  rownames(L2) <- paste0("MR", 1:nrow(L2))
}

which_fac_second <- apply(abs(L2), 1, which.max)
max_abs_second <- apply(abs(L2), 1, max)
fac_name_second <- colnames(L2)[which_fac_second]
communalities_second <- rowSums(L2^2)

primary_factor_second <- ifelse(max_abs_second >= 0.30 & communalities_second >= 0.25,
                                fac_name_second, NA_character_)

keep_second <- !is.na(primary_factor_second)
factors_by_second_order <- split(rownames(L2)[keep_second], primary_factor_second[keep_second])

cat("Hierarchical structure:\n")
for(second_factor in names(factors_by_second_order)) {
  cat("\n", toupper(second_factor), ":\n")
  for(first_factor in factors_by_second_order[[second_factor]]) {
    if(first_factor %in% names(items_by_factor)) {
      cat("  ", first_factor, " (", length(items_by_factor[[first_factor]]), " items)\n")
    }
  }
}

cat("\n=== STEP 10: ALIGNMENT QUALITY ASSESSMENT ===\n")

assess_alignment <- function(original_list, aligned_list, target) {
  tucker_congruence <- sapply(aligned_list, function(L) {
    diag(L %*% t(target) / sqrt(diag(L %*% t(L)) * diag(target %*% t(target))))
  })
  
  list(tucker_congruence = tucker_congruence,
       mean_congruence = apply(tucker_congruence, 1, mean))
}

first_order_quality <- assess_alignment(all_loadings_raw, 
                                        iterative_result$aligned_loadings, 
                                        pooled_loadings_procrustes)

cat("First-order alignment quality (Tucker's congruence by factor):\n")
print(round(first_order_quality$mean_congruence, 3))

second_order_quality <- assess_alignment(all_second_order_loadings,
                                         second_order_iterative$aligned_loadings,
                                         pooled_loadings_second_order)

cat("\nSecond-order alignment quality (Tucker's congruence by factor):\n")
print(round(second_order_quality$mean_congruence, 3))

cat("\n=== STEP 11: FINAL RESULTS AND EXPORT ===\n")

cat("ANALYSIS SUMMARY:\n")
cat("- Items analyzed:", length(ordinal_items), "\n")
cat("- Imputations:", length(imputation_ids), "\n")
cat("- MAP-suggested factors:", optimal_factors, "\n")
cat("- Items assigned to factors:", sum(sapply(items_by_factor, length)), "\n")
cat("- Items unassigned:", length(unassigned), "\n")
cat("- Second-order factors:", length(names(factors_by_second_order)), "\n")

# Export results
write.csv(pooled_loadings_procrustes, "first_order_loadings_enhanced_procrustes.csv")
write.csv(pooled_loadings_second_order, "second_order_loadings_enhanced_procrustes.csv")
write.csv(pooled_phi, "factor_correlations_enhanced.csv")
write.csv(second_order_pooled_phi, "second_order_factor_correlations.csv")


cat("Procrustes analysis complete! Results exported to CSV files.\n")