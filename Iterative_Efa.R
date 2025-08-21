# Packages

# This will automatically load your packages if they aren't loaded yet. 
packages <- c("tidyverse", "psych")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Reads the stacked imputed dataset. 
data_raw <- read.csv("completed_data_for_analysis_final.csv")

data <- data_raw %>%
  dplyr::select(-any_of(c("ID", "Sex","Race","Age","Education", "Income", "Sexual_Orientation", 
  "Status_SAQ", "Generalized_Anxiety", "Emotional_Problems_Drinking", "Desire_Drinking", 
  "Month_or_more_Drink", "Need_Drink_Same", "Alcohol_Affected_Work", "Anhedonia", 
  "Depressed_Affect", "Life_Day_By_Day", "Frustrated_Freq", "Irritable_Freq", ".id")))

# Seperate the ordinal Data
ordinal_items <- setdiff(names(data), ".imp")

# Convert to numeric then to ordered factors for polychoric correlations.
Imputed_Data <- data %>%
  mutate(across(any_of(ordinal_items),
                ~ as.integer(round(as.numeric(.))))) %>%
  mutate(across(any_of(ordinal_items),
                ~ ordered(., levels = sort(unique(.))))) %>%
  select(any_of(c(".imp", ordinal_items))) %>%
  mutate(.imp = as.integer(.imp))
  
# Get unique imputation IDs
imputation_ids <- unique(Imputed_Data$.imp)

# This function runs parallel analysis (fa.parallel()) to determine number of factors. 
# Runs factor analysis (fa() from psych) with that number of factors. Identifies the worst 
# item (lowest communality or loading). Removes one item per loop until all items meet thresholds: 
# loading_threshold = |0.30|, communal_threshold = 0.25. Logs everything per iteration: Number of 
# suggested factors, item removed, reason for removal, communality, loading, and complexity. 
# Returns a log as a tibble for each imputation.

run_iterative_efa <- function(data_subset, id) {
  
  items <- colnames(data_subset)
  loading_threshold <- 0.30 #Checks against absolute value later in script
  communal_threshold <- 0.25
  iteration <- 1
  log_list <- list()
  
  # Compute polychoric ONCE here, reuse via subsetting
  pc <- psych::polychoric(as.data.frame(data_subset), correct = 0)
  rho_full <- pc$rho
  if (is.null(colnames(rho_full))) {
    colnames(rho_full) <- rownames(rho_full) <- colnames(data_subset)
  }
  n_obs <- nrow(data_subset)
  
  repeat {
    # Subset the precomputed correlation matrix for current items
    Rsub <- rho_full[items, items, drop = FALSE]
    
    # Run parallel analysis on the correlation matrix (needs n.obs)
    pa <- fa.parallel(Rsub,
                      fa = "fa", fm = "minres",
                      n.obs = n_obs,
                      n.iter = 50, error.bars = FALSE, plot = FALSE)
    nf <- pa$nfact
    
    # Fit FA from the correlation matrix
    fm <- fa(r = Rsub,
             nfactors = nf,
             fm = "minres",
             rotate = "oblimin",
             n.obs = n_obs)
    
    lambda <- unclass(fm$loadings) %>%
      as.data.frame() %>%
      rownames_to_column(var = "rowname")
    
    lambda$max_loading <- apply(abs(lambda[,-1]), 1, max)
    
    communalities <- fm$communality
    idx_min_comm <- which.min(communalities)
    item_min_comm <- names(communalities)[idx_min_comm]
    min_comm <- communalities[idx_min_comm]
    
    idx_min_load <- which.min(lambda$max_loading)
    item_min_loading <- lambda$rowname[idx_min_load]
    min_loading <- lambda$max_loading[idx_min_load]
    
    if (min_comm < communal_threshold) {
      remove_var <- item_min_comm
      reason <- "communality"
      removal_value <- min_comm
    } else if (min_loading < loading_threshold) {
      remove_var <- item_min_loading
      reason <- "max_loading"
      removal_value <- min_loading
    } else {
      break
    }
    
    this_max_loading <- lambda$max_loading[match(remove_var, lambda$rowname)]
    log_list[[iteration]] <- tibble(
      imputation = id,
      iteration = iteration,
      suggested_factors = nf,
      item_removed = remove_var,
      removal_reason = reason,
      removal_value = round(removal_value, 2),
      communality = round(communalities[remove_var], 2),
      complexity = round(fm$complexity[remove_var], 2),
      highest_loading = round(this_max_loading, 2)
    )
        
    items <- setdiff(items, remove_var)
    iteration <- iteration + 1
  }
  
  bind_rows(log_list)
}

efa_logs <- list()

for (i in imputation_ids) {
  message("Running on imputation ", i)
  
  this_data <- Imputed_Data %>%
    filter(.imp == i) %>%
    select(-.imp)  # drop tracking cols
  
  efa_log <- run_iterative_efa(this_data, id = i)
  efa_logs[[i]] <- efa_log
}

# Combine all logs
combined_log <- bind_rows(efa_logs)

# Save 
write.csv(combined_log, "EFA_removal_log_all_imputations_final.csv", row.names = FALSE)
