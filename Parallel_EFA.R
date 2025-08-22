#Packages 
library(dplyr)
library(mice)
library(psych)
library(lavaan)
library(semTools)

# Load in the data.
  data <- read.csv("/Users/march341/Desktop/MIDUS/completed_data_for_analysis_final.csv")

# Remove the items that failed data prep.
  data <- data %>%
    dplyr::select(-any_of(c("ID", "Sex","Race","Age","Education", "Income", "Sexual_Orientation", 
                            "Status_SAQ", "Generalized_Anxiety", "Emotional_Problems_Drinking", "Desire_Drinking", 
                            "Month_or_more_Drink", "Need_Drink_Same", "Alcohol_Affected_Work", "Anhedonia", 
                            "Depressed_Affect", "Life_Day_By_Day", "Frustrated_Freq", "Irritable_Freq", ".id")))
  

# Separate the ordinal data.
  ordinal_items <- setdiff(names(data), ".imp")

# Convert to numeric then to ordered factors for polychoric correlations.
  data_num <- data %>%
    mutate(across(any_of(ordinal_items),
                  ~ as.integer(round(as.numeric(.))))) %>%
    mutate(across(any_of(ordinal_items),
                  ~ ordered(., levels = sort(unique(.))))) %>%
    dplyr::select(any_of(c(".imp", ordinal_items))) %>%
    mutate(.imp = as.integer(.imp))


# Split by imputation
  imp_slices <- split(data_num, data_num$.imp)

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

# ----------------- EFA --------------------

# Get unique imputation IDs
  imputation_ids <- unique(data_num$.imp)

# Set up lists to store results. 
  all_loadings <- list()
  all_communalities <- list()
  all_phis <- list()

# Loop through each imputed data set.
  for (imp in imputation_ids) {
    temp <- data_num %>%
      filter(.imp == imp) %>%
      dplyr::select(-.imp)
  
  # Calculate polychoric correlation matrix for each imputed data set.
    cor_result <- psych::polychoric(as.data.frame(temp[, ordinal_items, drop = FALSE]), correct = 0)
    cor_matrix <- cor_result$rho
  
  # Run EFA
    efa_result <- fa(cor_matrix, nfactors = 5, fm = "minres", rotate = "oblimin", scores = "tenBerge")
  
  # Store loading's and factor correlations
    all_loadings[[imp]] <- as.matrix(efa_result$loadings)
    all_communalities[[imp]] <- setNames(efa_result$communality, rownames(efa_result$loadings))
    all_phis[[imp]] <- efa_result$Phi
}

# Pool loading's
  pooled_loadings <- Reduce("+", all_loadings) / length(all_loadings)

# Pool communality 
  pooled_communality <- Reduce("+", all_communalities) / length(all_communalities)

# Pool factor correlations
  pooled_phi <- Reduce("+", all_phis) / length(all_phis)

# ------------ Now to sort every item into it's factor based on the highest loading ---------

# Determine majority-rule item assignment to factors
  primary_factor <- apply(abs(pooled_loadings), 1, function(x) {
    colnames(pooled_loadings)[which.max(x)]
  })

# Assume the items are in the same order as your original data columns
  rownames(pooled_loadings) <- colnames(data_num)[!colnames(data_num) %in% c(".imp")]


# L = your pooled loading matrix
  L <- as.matrix(pooled_loadings)

# Find each item's primary factor and the size of that loading
  which_fac <- apply(abs(L), 1, which.max)           # index of primary factor
  max_abs   <- apply(abs(L), 1, max)                  # magnitude of primary loading
  fac_name  <- colnames(L)[which_fac]                 # factor names

# Communalities = sum of squared loadings per item
  communalities <- rowSums(L^2)

# Enforce BOTH criteria: |loading| ≥ .30 AND communality ≥ .25
  primary_factor <- ifelse(max_abs >= 0.30 & communalities >= 0.25,
                           fac_name, NA_character_)

# Split only items that meet both criteria
  keep <- !is.na(primary_factor)
  items_by_factor <- split(rownames(L)[keep], primary_factor[keep])

  print(items_by_factor)

# See which items were dropped (by either cutoff)
  unassigned <- rownames(L)[!keep]
  unassigned

