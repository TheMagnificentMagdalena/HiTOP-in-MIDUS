# Packages ---------------------------------------------------------------

packages <- c("tidyverse", "lavaan", "semTools", "lavaan.mi")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Read stacked data & drop items that didn't pass the data preperation tests ---------------------------------

data_raw <- read.csv("completed_data_for_analysis_final.csv")
    
    data <- data_raw %>%
      dplyr::select(-any_of(c("ID", "Sex","Race","Age","Education", "Income", "Sexual_Orientation", 
                              "Status_SAQ", "Generalized_Anxiety", "Emotional_Problems_Drinking", "Desire_Drinking", 
                              "Month_or_more_Drink", "Need_Drink_Same", "Alcohol_Affected_Work", "Anhedonia", 
                              "Depressed_Affect", "Life_Day_By_Day", "Frustrated_Freq", "Irritable_Freq", ".id")))
                          
# Separate the ordinal Data

    ordinal_items <- setdiff(names(data), ".imp")

# Convert to numeric then to ordered factors for polychoric correlations.

data_num <- data %>%
  mutate(across(any_of(ordinal_items),
                ~ as.integer(round(as.numeric(.))))) %>%
  mutate(across(any_of(ordinal_items),
                ~ ordered(., levels = sort(unique(.))))) %>%
  dplyr::select(any_of(c(".imp", ordinal_items))) %>%
  mutate(.imp = as.integer(.imp))

# CFA model --------------------------------------------------------------

cfa_model <- '
  Internalizing =~ 
    Moody + Worrying + Nervous + Irritabiltiy_Frequency + Hot_Flashes +
    Falling_Asleep + Worry_Frequency + Nervous_Freq_30_Day + Restless_Freq_30_Day +
    Afraid_Freq_30_Day + Jittery_Freq_30_Day + Irritable_Freq_30_Day +
    Upset_Think_Day + Minor_Setback_Irritate + Frustrated_Freq_30_Day +
    Upset_Freq_30_Day + Mood_Up_Down + Nothing_Cheer + Hopeless_Freq +
    Everything_Effort + Worthless_Freq + Ashamed_Freq + Felt_Belong_Freq +
    Enthusiastic_Freq + Full_Of_Life_Freq + Extreme_Happy_Freq +
    Good_Spirits_Freq + Cheerful_Freq + Lonely_Frequency + Angry_Freq

  Disinhibited_Externalizing =~ 
    Responsible + Exciting_Look_Forward + Sense_Of_Direction +
    Good_Manage_Responsibility + Actively_Carry_Plans + Weigh_Possibilities +
    Like_Think_Over + Keep_Work_Problems + Like_Difficult_Things +
    High_Standards + Plan_Future + Enjoy_Plans_Future + Positive_About_Self +
    Will + Change_For_Better + Too_Much_Get_Done + Dont_Give_Up +
    Rarely_Give_Up + Know_Want_Life + Helpful_Set_Future_Goal +
    Goal_Keep_Benefits + Avoid_Distraction + Others_Turn_To_Me

  Antagonistic_Externalizing =~ 
    Angry_Ready_Hit + Effective_Talking_People + Good_Influence_People +
    Enjoy_Hurting_Say_Mean + When_Insult_Get_Even + Sometimes_Just_Hit_Someone

  Detachment =~ 
    People_Mean + Others_More_Life + Disapoint_Acheivment + Gave_Up +
    If_Go_Wrong_It_Will + Hardly_Expect_Go_My_Way + Rarely_Good_Happen +
    No_Good_Sense + Overwhelmed_Responsibility + Cautious_Person +
    Difficult_Arrange + Maintaining_Relationships_Hard + Dont_Fit_Community +
    Few_Close_Friends + No_Warm_Relationship + Difficult_Voice_Opinion +
    Self_Attitude + Dissapointment_Low_Goal + Too_Many_Things +
    No_Sense_Planning + More_Sucessful + Dont_Like_Help +
    Dont_Feel_Community + People_Take_Advantage

  Warmth =~ 
    Enjoy_Convo + I_Am_Warm_Person + Prefer_No_Others + Felt_Close_Others_Freq +
    Most_See_Love + Others_Describe_Giving + Important_Sympathetic +
    Feel_Close_Others
'

# Exact indicators used in the model

model_items <- c(
  "Moody","Worrying","Nervous","Irritabiltiy_Frequency","Hot_Flashes",
  "Falling_Asleep","Worry_Frequency","Nervous_Freq_30_Day","Restless_Freq_30_Day",
  "Afraid_Freq_30_Day","Jittery_Freq_30_Day","Irritable_Freq_30_Day",
  "Upset_Think_Day","Minor_Setback_Irritate","Frustrated_Freq_30_Day",
  "Upset_Freq_30_Day","Mood_Up_Down","Nothing_Cheer","Hopeless_Freq",
  "Everything_Effort","Worthless_Freq","Ashamed_Freq","Felt_Belong_Freq",
  "Enthusiastic_Freq","Full_Of_Life_Freq","Extreme_Happy_Freq",
  "Good_Spirits_Freq","Cheerful_Freq","Lonely_Frequency","Angry_Freq",
  "Responsible","Exciting_Look_Forward","Sense_Of_Direction",
  "Good_Manage_Responsibility","Actively_Carry_Plans","Weigh_Possibilities",
  "Like_Think_Over","Keep_Work_Problems","Like_Difficult_Things",
  "High_Standards","Plan_Future","Enjoy_Plans_Future","Positive_About_Self",
  "Will","Change_For_Better","Too_Much_Get_Done","Dont_Give_Up",
  "Rarely_Give_Up","Know_Want_Life","Helpful_Set_Future_Goal",
  "Goal_Keep_Benefits","Avoid_Distraction","Others_Turn_To_Me",
  "Angry_Ready_Hit","Effective_Talking_People","Good_Influence_People",
  "Enjoy_Hurting_Say_Mean","When_Insult_Get_Even","Sometimes_Just_Hit_Someone",
  "People_Mean","Others_More_Life","Disapoint_Acheivment","Gave_Up",
  "If_Go_Wrong_It_Will","Hardly_Expect_Go_My_Way","Rarely_Good_Happen",
  "No_Good_Sense","Overwhelmed_Responsibility","Cautious_Person",
  "Difficult_Arrange","Maintaining_Relationships_Hard","Dont_Fit_Community",
  "Few_Close_Friends","No_Warm_Relationship","Difficult_Voice_Opinion",
  "Self_Attitude","Dissapointment_Low_Goal","Too_Many_Things",
  "No_Sense_Planning","More_Sucessful","Dont_Like_Help",
  "Dont_Feel_Community","People_Take_Advantage",
  "Enjoy_Convo","I_Am_Warm_Person","Prefer_No_Others","Felt_Close_Others_Freq",
  "Most_See_Love","Others_Describe_Giving","Important_Sympathetic",
  "Feel_Close_Others"
)

# Keep only model items
    ordered_items_in_model <- intersect(ordinal_items, model_items)

# Rebuild data_cfa using ONLY the model items

    data_cfa <- data_num %>%
      mutate(across(all_of(ordered_items_in_model),
                    ~ ordered(., levels = sort(unique(.))))) %>%
      dplyr::select(any_of(".imp"), all_of(ordered_items_in_model)) %>%
      mutate(.imp = as.integer(.imp))

# Split stacked data to list ---------------------------------------------

imp_list <- split(data_cfa, data_cfa$.imp)
imp_list <- lapply(imp_list, function(df) df %>% select(-.imp))

# Fit CFA with lavaan.mi -------------------------------------------------

    fit_mi <- lavaan.mi::cfa.mi(
      model     = cfa_model,
      data      = imp_list,
      estimator = "WLSMV",
      parameterization = "theta",
      std.lv    = TRUE,
      ordered   = ordered_items_in_model
    )

# Outputs and summary ----------------------------------------------------

{
  sum_txt <- capture.output({
    cat("=== lavaan.mi summary (no extra args) ===\n")
    print(summary(fit_mi))
    cat("\n=== Robust/scaled fit measures (WLSMV) ===\n")
    print(fitMeasures(fit_mi, c(
      "chisq","df","pvalue",                
      "chisq.scaled","df.scaled","pvalue.scaled",  
      "cfi","tli","rmsea", "srmr",
      "cfi.scaled","tli.scaled","rmsea.scaled"
    )))
  })
  writeLines(sum_txt, "pooled_5_CFA_summary.txt")
}

# Extract parameters 

    pe <- parameterEstimates.mi(fit_mi, se = TRUE, ci = TRUE, level = 0.95,
    rsquare = TRUE, standardized = TRUE, output = "data.frame")

write.csv(pe, "pooled_CFA_parameters_standardized.csv", row.names = FALSE)

# --- Reliability & convergent validity ---------------------------------

# Factor-level reliability (alpha, omega, etc.)
    rel <- semTools::reliability(fit_mi)
    rel_df <- as.data.frame(rel)
    rel_df$factor <- rownames(rel_df)
    rel_df <- rel_df[, c("factor", setdiff(names(rel_df), "factor"))]
    write.csv(rel_df, "pooled_CFA_reliability.csv", row.names = FALSE)


# Modification Indices for MI -------------------------------------------
    mi <- lavaan.mi::modindices.mi(fit_mi, sort. = TRUE)
    write.csv(mi, "pooled_CFA_modindices.csv", row.names = FALSE)
    
    cat("\nDone. All outputs written successfully.\n")

