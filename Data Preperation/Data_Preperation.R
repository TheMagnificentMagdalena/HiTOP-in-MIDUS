#Packages 
library(psych)
library(dplyr)

# ------------- Load in the data --------------
Data_Analysis <- read.csv("/Users/march341/Desktop/MIDUS/completed_data_for_analysis_final_2.csv")

# ------------- Check the Imputed Data --------------
#Some visualization
  str(Data_Analysis)
  summary(Data_Analysis)

#Verify imputation count and equal rows in each imputation
  table(Data_Analysis$.imp)

#Check for missing values left over
  Data_Analysis %>%
    group_by(.imp) %>%
    summarise(across(everything(), ~ mean(is.na(.))))

### Data set is too big for visual diagnostics ###

#Inspect ID's to ensure no doubles
  Data_Analysis %>%
    group_by(.imp, ID) %>%
    tally(name = "n") %>%
    filter(n > 1) %>%
    ungroup()
  # Expect 0 rows

#Check distributions don't wildly change between imputations
summ_stats <- Data_Analysis %>%
  group_by(.imp) %>%
  summarise(across(where(is.numeric), list(mean = ~mean(.), sd = ~sd(.))), .groups = "drop")

# Flag variables with large spread of means across imputations
summ_long <- tidyr::pivot_longer(summ_stats, - .imp,
                                 names_to = c("var",".value"),
                                 names_sep = "_")
by_var <- summ_long %>%
  group_by(var) %>%
  summarise(mean_of_means = mean(mean, na.rm=TRUE),
            sd_of_means   = sd(mean,  na.rm=TRUE),
            mean_of_sds   = mean(sd,   na.rm=TRUE),
            sd_of_sds     = sd(sd,     na.rm=TRUE),
            .groups = "drop") %>%
  arrange(desc(sd_of_means))
head(by_var, 15)  # inspect top few variables

# -------------  Check Out the Covariance Matrix ---------------------

# Make the data numeric because most of these checks require it. 
# Keep only imputation 1
Data_Analysis_1 <- Data_Analysis %>%
  dplyr::filter(.imp == 1)

# Remove identifier columns (ID, .imp, .id if present), updated to remove variables that failed the tests in this script.
num_data <- Data_Analysis_1 %>%
  dplyr::select(-ID, -.imp, -.id, -Income, -Age, -Education, -Sexual_Orientation, -Status_SAQ, -Race, -Generalized_Anxiety, -Emotional_Problems_Drinking, -Desire_Drinking, -Month_or_more_Drink, -Need_Drink_Same, -Alcohol_Affected_Work, -Anhedonia, -Depressed_Affect, -Life_Day_By_Day, -Frustrated_Freq, -Irritable_Freq) %>%
  dplyr::mutate(across(everything(), ~ as.numeric(as.character(.))))

# Compute covariance matrix. 
  cor_matrix <- cor(num_data)

# Eigenvalues, make sure none are 0 or smaller! 
  eigenvalues <- eigen(cor_matrix)$values
  print(eigenvalues, digits = 2)

# Determinant, also shouldn't be 0 or negative, we are playing it dangerously close here. 
  det_value <- det(cor_matrix)
  print(det_value)


# Make sure rank is larger than number of variables to ensure it is full rank.
qr(cor(num_data))$rank

# Number of variables 
ncol(num_data)

## We must investigate further and re-run diagnostics, things are amiss! Update: Things are still amiss but we proceed anyway!

# -------------  Check for 0 Variance -------------
# Compute near 0 variance 95/5
  nzv_cols <- nearZeroVar(num_data) #Prints columns with near 0 variance 

# Get column names instead of numbers
  nzv_colnames <- names(num_data)[nzv_cols]

# Print the names
  print(nzv_colnames) 

#Items that need to be removed from this step: Generalized_Anxiety, Emotional_Problems_Drinking, Desire_Drinking, Month_or_more_Drink, Need_Drink_Same, Alcohol_Affected_Work, Anhedonia, Depressed_Affect 

# -------------- Check Multicollinearity and Linear Dependence --------------

# Continuous
  c_vars <- c(
    "Generalized_Anxiety",  # 1–10
    "Depressed_Affect"     # 0–7 composite
  )

# Dichotomous
  d_vars <- c(
    "Emotional_Problems_Drinking",
    "Desire_Drinking",
    "Month_or_more_Drink",
    "Need_Drink_Same"
  )

#Make Binary 0/1
num_data[d_vars] <- lapply(num_data[d_vars], function(x) ifelse(x == min(x, na.rm=TRUE), 0L, 1L))

# Polytomous (everything else that remains is Likert-type)
p_vars <- setdiff(names(num_data), c(c_vars, d_vars))

# Calcualte polychoric correlations matrix
Cor <- mixedCor(data=num_data,c=NULL,p=p_vars,d=NULL,smooth=TRUE,correct=.5,global=TRUE,ncat=8,
                       weight=NULL)

  Cor_Matrix <- Cor$rho
  smc <- smc(Cor_Matrix)  

#SMCs > .90
  high_smc <- which(smc > .90)
  high_smc 
   
  ### Variables to remove from this step: Life_Day_By_Day, Frustrated_Freq, Month_or_more_Drink, Irritable_Freq ###
  
#Now we will look at the variance inflation factor and remove any variables with VIF > 10.
#Invert correlation matrix
  invR <- solve(Cor_Matrix)

#Extract VIFs
  vif_vals <- diag(invR) 

# The diagonal elements of the inverse of the correlation matrix is the VIF. When you invert that matrix, the diagonal elements  turn out to be 1/1 - R squared which is exactly the definition of the variance-inflation factor.

  print(vif_vals) # Nothing to be removed from this stage

#Finding linear combinations
  findLinearCombos(num_data) #Nothing to be removed from this stage,

# --------------- Check for too Low of an Endorsement Rate -----------------
#Set your threshold  
  threshold <- 0.95
  
  # Compute, for each column, the maximum category frequency
  max_props <- sapply(num_data, function(x) 
    max(prop.table(table(x))))
  
  # Pick the names where that max is greater than threshold
  high_freq_items <- names(max_props)[max_props >= threshold]
  
  high_freq_items #Items to remove from this round: Generalized_Anxiety, Emotional_Problems_Drinking, Desire_Drinking, Month_or_more_Drink, Need_Drink_Same, Anhedonia
  
# -------------- Check Variance Ratios -----------------
  
#Covariance matrices in which the ratio of the largest to the smallest variance is greater than 10 are ill scaled.
  
  # Column variances
  variances <- apply(num_data, 2, stats::var)
  
  # Guard against zeros/NA
  pos_min <- min(variances[variances > 0])
  max_var <- max(variances, na.rm = TRUE)
  variance_ratio <- max_var / pos_min
  
  variance_ratio
  sort(variances, decreasing = TRUE)[1:5]   # top 5 largest
  sort(variances, decreasing = FALSE)[1:5]  # top 5 smallest
  
  # Variance ratio is a little bit higher than we want, but we are going to roll with it. 

# -------------- Checking Multivariate Normality ----------------
  
  # Perform Mardia's Test
  mardia <- mvn(data = num_data, mvnTest = "mardia")
  mardia #Fails Mardias test, use robust statistics 
  
  # Per-variable skew 
  skews <- sapply(num_data, function(col) psych::skew(col, na.rm = TRUE))
  
  # Flag extreme skew
  bad_skew <- names(skews[abs(skews) > 3])
  bad_skew
  
  # Compute excess kurtosis for each variable. 
  excess_kurtosis <- sapply(num_data, function(x) {
    k <- kurtosis(x, na.rm=TRUE)
    k - 3
  })
  
  # Print those little nuggets
  excess_kurtosis
  
  # Flag any with excess kurtosis > 8 for diagnosis of really bad kurtosis. 
  bad_kurtosis <- names(excess_kurtosis[abs(excess_kurtosis) > 8])
  bad_kurtosis
  
   #Luckily these tests tell us to remove the same items we removed previously

 # ------------- If You Need to Smooth the Matrix ----------------
  
  #Compute the raw covariance matrix
  S <- cov(data_subset)
  
  #Identify any zero‐variance items
  zero_vars <- which(apply(data_subset, 2, stats::var, na.rm = TRUE) == 0)
  
  #“patch” those rows/columns so constant items have zero covariance with everything
  S[ zero_vars, ] <- 0
  S[ , zero_vars] <- 0
  #Ensure their variance on the diagonal is zero
  diag(S)[zero_vars] <- 0
  
  #Nudge to the nearest positive‐definite covariance matrix
  S_pd <- as.matrix(nearPD(S)$mat )
  
  #Quantify how much you’ve changed it
  max_abs_change <- max(abs(S_pd - S))
  sum_sq_change  <- sum((S_pd - S)^2)
  
# --------------- Check if the Data Should be Factor Analyzed ---------------
  
  Bartlette <- cortest.bartlett(num_data)
  Bartlette # Most large datastes fail this one
  
  KMO <- KMO(Cor_Matrix)
  KMO # Technically all items pass this test 
  