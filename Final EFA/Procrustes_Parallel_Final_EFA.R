# Packages 
library(dplyr)
library(psych)
library(tidyr)
library(mice)

removal_log <- read.csv("/Users/march341/Desktop/Final EFA/EFA_removal_log_all_imputations_final.csv")

# Count % of imputations where each item was removed
percent_removed <- removal_log %>%
  distinct(imputation, item_removed) %>%         # one row per imputation-item
  count(item_removed, name = "n_imputations") %>%
  mutate(total_imps = n_distinct(removal_log$imputation),
         percent = 100 * n_imputations / total_imps)

percent_removed

# We will choose to remove items that were removed across 50% or more imputations. 
data_raw <- read.csv("/Users/march341/Desktop/MIDUS/completed_data_for_analysis_final.csv")

data <- data_raw %>%
  dplyr::select(-c(ID, Sex, .id, Income, Age, Education, Sexual_Orientation, Status_SAQ, Race, Generalized_Anxiety, Emotional_Problems_Drinking, Desire_Drinking, Month_or_more_Drink, Need_Drink_Same, Alcohol_Affected_Work, Anhedonia, Depressed_Affect, Life_Day_By_Day, Frustrated_Freq, Irritable_Freq, Actively_Carry_Plans, Alcohol_More_Than_Intended, Approach_Health, Avoid_Distraction, Aware_Body, Bad_Happens_Prevented, Careless, Choose_Goals, Dont_Fit_Community, Everything_Effort, Exciting_Look_Forward, Expect_The_Best, Hate_Temp, Headache_Frequency, High_Standards, Intercourse_Pain, Leaking_Urine, Loud_Noise, Low_Pain_Tolerance, Mental_Emotional, Minor_Setback_Irritate, More_Or_Less_Positive,  More_Sucessful, Nothing_Cheer, Optimistic_About_Future, People_Take_Advantage, Positive_About_Self, Quick_Sense_Hunger, Releived_No_Responsibility, Sense_Of_Direction, Should_Obey_Law, Social_Others_Lead, Some_Wander_Aimlessly, Sympathy_Limits, To_Reach_Goals, Trust_Friends, Upset_Think_Day, When_Difficult, Without_skill, Worry_Compare, Worry_Frequency, Worthless_Freq, Dissapointment_Low_Goal, Responsible, Live_One_Day, Asking_Help_Naturally, Enjoy_Plans_Future, Felt_Close_Others_Freq, Felt_Belong_Freq, Enthusiastic_Freq, Full_Of_Life_Freq, Extreme_Happy_Freq, Good_Spirits_Freq, Cheerful_Freq, Effective_Talking_People,  Good_Influence_People, Others_Turn_To_Me)) #Dissapointment_Low_Goal, Asking_Help_Naturally, Enjoy_Plans_Future, Live_One_Day, and responsible were removed to reduce doublets. Felt_Close_Others_Freq, Felt_Belong_Freq, Enthusiastic_Freq, Full_Of_Life_Freq, Extreme_Happy_Freq, Good_Spirits_Freq, and Cheerful_Freq were removed because they are outside of the scope of HiTOP. Effective_Talking_People, Good_Influence_People, and Others_Turn_To_Me were removed because they are just assertiveness and do not group meaningfully with any other factor.  

# Setting up the items that are to be treated as ordinal.
ordinal_items <- setdiff(names(data), ".imp")

# Convert to numeric, then set up the factors as ordinal. Ensure that .imp is still numeric. 
data_clean <- data %>%
  mutate(across(any_of(ordinal_items),
                ~ as.integer(round(as.numeric(.))))) %>%
  mutate(across(any_of(ordinal_items),
                ~ ordered(., levels = sort(unique(.))))) %>%
  dplyr::select(any_of(c(".imp", ordinal_items))) %>%
  mutate(.imp = as.integer(.imp))

# Split by imputation
imp_slices <- split(data_clean, data_clean$.imp)

# ----------------- Parallel analysis ----------------------

# Polychoric per imputation (exclude .imp column)
cor_list <- lapply(imp_slices, function(df) {
  pc <- psych::polychoric(as.data.frame(df[, ordinal_items, drop = FALSE]), correct = 0)
  pc$rho
})

# Simplified Rubin's rules for pooling the correlation matrices.  
pool_corrs <- function(cor_list) {
  A    <- simplify2array(cor_list)      
  Z    <- atanh(A)                      # Fisher z
  Zbar <- apply(Z, c(1,2), mean)        # Average
  Rbar <- tanh(Zbar)                    # Back-transform
  diag(Rbar) <- 1                       # Reset diagonals
  dimnames(Rbar) <- dimnames(cor_list[[1]])
  Rbar
}

# We will run parallel analysis on the pooled matrix, but ONLY for parallel analysis. It is bad practice to perform analyses on the pooled matrix, rather than perform the computation on each imputation and pool across the imputations. We are only doing this because there is no way to run by imputation and pool with parallel analysis. 
R_pooled <- pool_corrs(cor_list)


# Parallel on R_pooled 
fa.parallel(R_pooled,   
            n.obs    = 7350,  #By imputation sample size.        
            fa       = "fa",       
            fm       = "minres",   
            n.iter   = 50,        
            error.bars = TRUE,     
            main     = "Parallel Analysis on Pooled R"
)

# --------------------------- EFA using procrustes ----------------------

# ===================================================================
# PROCRUSTES EFA WITH HIERARCHICAL ANALYSIS AND PARALLEL ANALYSIS
# ===================================================================


cat("=== STEP 1: FIRST-ORDER EFA ON EACH IMPUTATION ===\n")
# This step runs separate EFAs on each imputed dataset

# Get unique imputation IDs
imputation_ids <- unique(data_clean$.imp)
cat("Number of imputations found:", length(imputation_ids), "\n")

# Set up lists to store results from each imputation
all_loadings_raw <- list()        # Raw loadings before alignment
all_loadings_aligned <- list()    # Loadings after Procrustes alignment
all_communalities <- list()       # Item communalities (R²)
all_phis <- list()               # Factor correlation matrices
all_efa_fits <- list()           # Full EFA objects (for diagnostics)

# Run EFA on each imputation separately
for (imp in imputation_ids) {
  cat("Processing imputation", imp, "of", length(imputation_ids), "\n")
  
  # Extract data for this imputation only
  temp <- data_clean %>%
    filter(.imp == imp) %>%
    dplyr::select(-.imp)  # Remove imputation indicator
  
  # Calculate polychoric correlations
  cor_result <- psych::polychoric(as.data.frame(temp[, ordinal_items, drop = FALSE]), correct = 0)
  cor_matrix <- cor_result$rho
  
  # Run EFA on this imputation
  efa_result <- fa(cor_matrix, nfactors = 10, fm = "minres", rotate = "oblimin")
  
  # Store all results for later pooling
  all_loadings_raw[[imp]] <- as.matrix(efa_result$loadings)
  all_communalities[[imp]] <- setNames(efa_result$communality, rownames(efa_result$loadings))
  all_phis[[imp]] <- efa_result$Phi 
  all_efa_fits[[imp]] <- efa_result
}

cat("\n=== STEP 2: PROCRUSTES ROTATION IMPLEMENTATION ===\n")
# Factors might emerge in different orders across imputations
# Factor 1 in imputation 1 might be Factor 3 in imputation 2
# Procrustes finds optimal rotation to align them

procrustes_rotation <- function(X, Y) {
  # X = matrix to be rotated (current imputation)
  # Y = target matrix (what we want to align to)
  # Mathematical solution: Find orthogonal matrix T such that XT ≈ Y
  # Uses Singular Value Decomposition (SVD) for optimal solution
  
  # SVD of Y'X gives us the components we need
  svd_result <- svd(t(Y) %*% X)
  
  # Optimal rotation matrix T = V * U'
  # This minimizes ||XT - Y||² (squared differences)
  T <- svd_result$v %*% t(svd_result$u)
  
  return(T)
}

# Function to align all loading matrices to a target
align_to_target <- function(loadings_list, target) {
  aligned_list <- list()
  
  for (i in seq_along(loadings_list)) {
    current_loadings <- loadings_list[[i]]
    
    # Safety check: matrices must be same size
    if (!all(dim(current_loadings) == dim(target))) {
      stop("Loading matrices must have same dimensions")
    }
    
    # Find optimal rotation to align current to target
    T <- procrustes_rotation(current_loadings, target)
    
    # Apply the rotation: rotated = original * rotation_matrix
    aligned_loadings <- current_loadings %*% T
    
    # Preserve row/column names (items and factors)
    dimnames(aligned_loadings) <- dimnames(current_loadings)
    aligned_list[[i]] <- aligned_loadings
  }
  
  return(aligned_list)
}

cat("\n=== STEP 3: CHOOSING TARGET AND ALIGNMENT ===\n")

iterative_procrustes <- function(loadings_list, max_iter = 10, tol = 1e-6) {
  cat("Starting iterative Procrustes alignment...\n")
  
  # Start with simple average as initial target
  target <- Reduce("+", loadings_list) / length(loadings_list)
  
  for (iter in 1:max_iter) {
    old_target <- target
    
    # Align all matrices to current target
    aligned <- align_to_target(loadings_list, target)
    
    # Update target to average of aligned matrices
    new_target <- Reduce("+", aligned) / length(aligned)
    
    # Check if we've converged (target isn't changing much)
    change <- sqrt(sum((new_target - old_target)^2))
    cat("Iteration", iter, "- Change:", round(change, 6), "\n")
    
    if (change < tol) {
      cat("Converged after", iter, "iterations\n")
      break
    }
    
    target <- new_target
  }
  
  return(list(
    final_target = target,
    aligned_loadings = aligned,
    n_iterations = iter
  ))
}

# Run iterative alignment 
iterative_result <- iterative_procrustes(all_loadings_raw)
pooled_loadings_first_order <- iterative_result$final_target

cat("\n=== STEP 4: POOLING OTHER FIRST-ORDER RESULTS ===\n")
# Pool communalities (simple average across imputations)
pooled_communality <- Reduce("+", all_communalities) / length(all_communalities)

# Pool factor correlation matrices (Phi matrices)
# These will be used for second-order analysis
pooled_phi <- Reduce("+", all_phis) / length(all_phis)

cat("First-order pooling complete.\n")
cat("Pooled factor correlation matrix dimensions:", dim(pooled_phi), "\n")

cat("\n=== STEP 5: FIRST-ORDER FACTOR ASSIGNMENT ===\n")
# Assign items to factors based on highest loading + criteria

L1 <- as.matrix(pooled_loadings_first_order)
rownames(L1) <- colnames(data_clean)[!colnames(data_clean) %in% c(".imp")]

# Find each item's primary factor
which_fac <- apply(abs(L1), 1, which.max)      # Which factor has highest loading
max_abs <- apply(abs(L1), 1, max)              # What is that loading value
fac_name <- colnames(L1)[which_fac]            # Factor names
communalities_first <- rowSums(L1^2)           # Sum of squared loadings per item

# Apply criteria: |loading| ≥ .30 AND communality ≥ .25
primary_factor_first <- ifelse(max_abs >= 0.30 & communalities_first >= 0.25,
                               fac_name, NA_character_)

# Split items by their assigned factor
keep_first <- !is.na(primary_factor_first)
items_by_first_factor <- split(rownames(L1)[keep_first], primary_factor_first[keep_first])

cat("First-order factors and their items:\n")
print(items_by_first_factor)

unassigned_first <- rownames(L1)[!keep_first]
cat("\nUnassigned items (first-order):", paste(unassigned_first, collapse = ", "), "\n")

cat("\n=== STEP 6: SECOND-ORDER ANALYSIS ===\n")

cat("Running parallel analysis on factor correlations...\n")
fa.parallel(pooled_phi, 
            n.obs = nrow(data_clean)/length(unique(data_clean$.imp)),  # Sample size per imputation
            fa = "fa", 
            fm = "minres", 
            plot = TRUE,
            main = "Parallel Analysis: Second-Order Factors")

# Based on parallel analysis, choose number of second-order factors
n_second_order <- 5

cat("Extracting", n_second_order, "second-order factors...\n")

# Set up storage for second-order analysis across imputations
all_second_order_loadings <- list()
all_second_order_phis <- list()

# Run second-order EFA on each imputation's Phi matrix
for (imp in imputation_ids) {
  phi_matrix <- all_phis[[imp]]
  
  # EFA on factor correlations
  second_order_efa <- fa(phi_matrix, 
                         nfactors = n_second_order, 
                         fm = "minres", 
                         rotate = "oblimin",
                         n.obs = nrow(data_clean)/length(unique(data_clean$.imp)))
  
  # Store results
  all_second_order_loadings[[imp]] <- as.matrix(second_order_efa$loadings)
  all_second_order_phis[[imp]] <- second_order_efa$Phi
}

cat("\n=== STEP 7: PROCRUSTES ALIGNMENT FOR SECOND-ORDER ===\n")

second_order_iterative <- iterative_procrustes(all_second_order_loadings)
pooled_loadings_second_order <- second_order_iterative$final_target

# Pool second-order factor correlations
pooled_phi_second <- Reduce("+", all_second_order_phis) / length(all_second_order_phis)

cat("Second-order pooling complete.\n")

cat("\n=== STEP 8: SECOND-ORDER FACTOR ASSIGNMENT ===\n")
# Assign first-order factors to second-order factors

L2 <- as.matrix(pooled_loadings_second_order)
# Row names should be first-order factor names (MR1, MR2, etc.)
if(is.null(rownames(L2))) {
  rownames(L2) <- paste0("MR", 1:nrow(L2))
}

# Find primary second-order factor for each first-order factor
which_fac_second <- apply(abs(L2), 1, which.max)
max_abs_second <- apply(abs(L2), 1, max)
fac_name_second <- colnames(L2)[which_fac_second]
communalities_second <- rowSums(L2^2)

# Apply criteria for second-order assignment 
primary_factor_second <- ifelse(max_abs_second >= 0.30 & communalities_second >= 0.25,
                                fac_name_second, NA_character_)

# Group first-order factors by second-order factors
keep_second <- !is.na(primary_factor_second)
factors_by_second_order <- split(rownames(L2)[keep_second], primary_factor_second[keep_second])

cat("Second-order factor structure:\n")
print(factors_by_second_order)

unassigned_second <- rownames(L2)[!keep_second]
cat("\nUnassigned first-order factors:", paste(unassigned_second, collapse = ", "), "\n")

cat("\n=== STEP 9: ASSESSMENT AND OUTPUT ===\n")

# Tucker's congruence (how well did factors align?)
assess_alignment <- function(original_list, aligned_list, target) {
  tucker_congruence <- sapply(aligned_list, function(L) {
    diag(L %*% t(target) / sqrt(diag(L %*% t(L)) * diag(target %*% t(target))))
  })
  
  list(
    tucker_congruence = tucker_congruence,
    mean_congruence = apply(tucker_congruence, 1, mean)
  )
}

# Assess first-order alignment
first_order_quality <- assess_alignment(all_loadings_raw, 
                                        iterative_result$aligned_loadings, 
                                        pooled_loadings_first_order)

cat("First-order alignment quality (Tucker's congruence by factor):\n")
print(round(first_order_quality$mean_congruence, 3))

# Assess second-order alignment
second_order_quality <- assess_alignment(all_second_order_loadings,
                                         second_order_iterative$aligned_loadings,
                                         pooled_loadings_second_order)

cat("\nSecond-order alignment quality (Tucker's congruence by factor):\n")
print(round(second_order_quality$mean_congruence, 3))

cat("\n=== FINAL HIERARCHICAL STRUCTURE ===\n")
cat("First-order factors (items):\n")
for(factor_name in names(items_by_first_factor)) {
  cat("", factor_name, ":", paste(items_by_first_factor[[factor_name]], collapse = ", "), "\n")
}

cat("\nSecond-order factors (first-order factors):\n")
for(second_factor in names(factors_by_second_order)) {
  cat("", second_factor, ":", paste(factors_by_second_order[[second_factor]], collapse = ", "), "\n")
  
  # Show items under each second-order factor
  for(first_factor in factors_by_second_order[[second_factor]]) {
    if(first_factor %in% names(items_by_first_factor)) {
      cat("    ->", first_factor, "items:", paste(items_by_first_factor[[first_factor]], collapse = ", "), "\n")
    }
  }
}

# Export results
write.csv(pooled_loadings_first_order, "first_order_loadings_procrustes.csv")
write.csv(pooled_loadings_second_order, "second_order_loadings_procrustes.csv")
write.csv(as.data.frame(items_by_first_factor), "first_order_factor_assignments.csv")
write.csv(as.data.frame(factors_by_second_order), "second_order_factor_assignments.csv")

cat("\nAnalysis complete! Results exported to CSV files.\n")