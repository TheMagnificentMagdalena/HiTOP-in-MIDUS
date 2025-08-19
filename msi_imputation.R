# Packages

packages <- c("dplyr", "tidyverse", "psych", "Matrix", "mice")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

library(dplyr)
library(tidyverse)
library(psych)
library(Matrix)
library(mice)


#Load in data
MIDUS2 <- load("/users/2/march341/multiple_imputation/Midlife in the United States (MIDUS 2), 2004-2006 (ICPSR 4652)/DS0001/04652-0001-Data.rda")
MIDUSrefresher <- load("/users/2/march341/multiple_imputation/Midlife in the United States (MIDUS Refresher 1), 2011-2014 (ICPSR 36532)/DS0001/36532-0001-Data.rda")
MIDUSrestricted36722 <- load("/users/2/march341/multiple_imputation/MIDUS Restricted Data/ICPSR_36722/DS0001/36722-0001-Data-REST.rda")
MIDUSrestricted22840 <- load("/users/2/march341/multiple_imputation/MIDUS Restricted Data/ICPSR_22840/DS0001/22840-0001-Data-REST.rda")

# Making sure the ID's match across the waves
# MIDUS2
da04652.0001 <- da04652.0001 |>
  mutate(ID = M2ID)

# Refresher sample 
da36532.0001 <- da36532.0001 |>
  mutate(ID = MRID)

# MIDUS restriced 22840

da22840.0001 <- da22840.0001 |>
  mutate(ID = M2ID)

# MIDUS refresher restricted 36722
da36722.0001 <- da36722.0001 |>
  mutate(ID = MRID)

# Create the final data set
subset <- bind_rows(da04652.0001,
                    da36532.0001, 
                    da22840.0001, 
                    da36722.0001) |> 
  mutate(
    ID = ID,
    Age = coalesce(B1PAGE_M2, BACRAGE, RA1PRAGE, RAACRAGE),
    Sex = coalesce(B1PRSEX, BACRSEX, RA1PRSEX, RAACRSEX),
    Education = coalesce(B1PB1, BACB1, RA1PB1, RAACB1),
    Race = coalesce(B1PF7A,BACF7A,RA1PF7A,RAACF7A),
    Sexual_Orientation = coalesce(B1SM6, BASG6, RA1SN6, RAASH6),
    Income = coalesce(B1STINC1, BACTINC1, RA1STINC, RAACTINC),
    Moody = coalesce(B1SE6C, BASC6C, RA1SF6C, RAASC7C),
    Organized = coalesce(B1SE6D, BASC6D, RA1SF6D, RAASC7D),
    Worrying = coalesce(B1SE6H, BASC6H, RA1SF6H, RAASC7H),
    Responsible = coalesce(B1SE6I, BASC6I, RA1SF6I, RAASC7I),
    Nervous = coalesce(B1SE6M, BASC6M, RA1SF6M, RAASC7M),
    Careless = coalesce(B1SE6X, BASC6X, RA1SF6X, RAASC7X),
    Aware_Body = coalesce(B1SA9A, BACAS9A, RA1SA9A, RAACAS9A),
    Loud_Noise = coalesce(B1SA9B, BACAS9B, RA1SA9B, RAACAS9B),
    Hate_Temp = coalesce(B1SA9C, BACAS9C, RA1SA9C, RAACAS9C),
    Quick_Sense_Hunger = coalesce(B1SA9D, BACAS9D, RA1SA9D, RAACAS9D),
    Low_Pain_Tolerance = coalesce(B1SA9E, BACAS9E, RA1SA9E, RAACAS9E),
    Headache_Frequency = coalesce(B1SA10A, BACAS10A, RA1SA10A, RAACAS10A),
    Backache_Frequency = coalesce(B1SA10B, BACAS10B, RA1SA10B, RAACAS10B),
    Swea_Frequency = coalesce(B1SA10C, BACAS10C, RA1SA10C, RAACAS10C),
    Irritabiltiy_Frequency = coalesce(B1SA10D, BACAS10D, RA1SA10D, RAACAS10D),
    Hot_Flashes = coalesce(B1SA10E, BACAS10E, RA1SA10E, RAACAS10E),
    Aches_Joints_Stiffness_Frequency = coalesce(B1SA10F, BACAS10F, RA1SA10F, RAACAS10F),
    Falling_Asleep = coalesce(B1SA10G, BACAS10G, RA1SA10G, RAACAS10G),
    Leaking_Urine = coalesce(B1SA10H, BACAS10H, RA1SA10H, RAACAS10H),
    Intercourse_Pain = coalesce(B1SA10I, BACAS10I, RA1SA10I, RAACAS10I),
    Extremities_Aches = coalesce(B1SA10J, BACAS10J, RA1SA10J, RAACAS10J),
    Generalized_Anxiety = coalesce(B1PANXIE, BACANXIE, RA1PANXIE, RAACANXIE),
    Worry_Compare = coalesce(B1PA83, BACA83, RA1PA83, RAACA83),
    Worry_Frequency = coalesce(B1PA84, BACA84, RA1PA84, RAACA84),
    Nervous_Freq_30_Day = coalesce(B1SA24B, BACAS24B, RA1SA20B, RAACAS24B),
    Restless_Freq_30_Day = coalesce(B1SA24C, BACAS24C, RA1SA20C, RAACAS24C),
    Afraid_Freq_30_Day = coalesce(B1SA24H, BACAS24H, RA1SA20H, RAACAS24H),
    Jittery_Freq_30_Day = coalesce(B1SA24I, BACAS24I, RA1SA20I, RAACAS24I),
    Irritable_Freq_30_Day = coalesce(B1SA24J, BACAS24J, RA1SA20J, RAACAS24J),
    People_Mean = coalesce(B1SE7S, BASC7S, RA1SF7S, RAASC8S), 
    Upset_Think_Day = coalesce(B1SE7W, BASC7W, RA1SF7W, RAASC8W), 
    Minor_Setback_Irritate = coalesce(B1SE7X, BASC7X, RA1SF7X, RAASC8X), 
    Frustrated_Freq_30_Day = coalesce(B1SA24N, BACAS24N, RA1SA20N, RAACAS24N),
    Upset_Freq_30_Day = coalesce(B1SA24L, BACAS24L, RA1SA20L, RAACAS24L),
    Others_More_Life = coalesce(B1SE1R, BASC1R, RA1SF1R, RAASC1R),
    Disapoint_Acheivment = coalesce(B1SE1DD, BACG5, RA1SF1DD, RAASC1DD),
    Gave_Up = coalesce(B1SE1GG, BASC1GG, RA1SF1GG, RAASC1GG),
    Mood_Up_Down = coalesce(B1SE7K, BASC7K, RA1SF7K, RAASC8K), 
    Exciting_Look_Forward = coalesce(B1SE7GG, BASC7GG, RA1SF7GG, RAASC8GG), 
    Expect_The_Best = coalesce(B1SE10A, BASC10A, RA1SF10A, RAASC11A),
    If_Go_Wrong_It_Will = coalesce(B1SE10B, BASC10B, RA1SF10B, RAASC11B),
    Optimistic_About_Future = coalesce(B1SE10C, BASC10C, RA1SF10C, RAASC11C),
    Hardly_Expect_Go_My_Way = coalesce(B1SE10D, BASC10D, RA1SF10D, RAASC11D),
    Rarely_Good_Happen = coalesce(B1SE10E, BASC10E, RA1SF10E, RAASC11E),
    Expect_Good_Not_Bad = coalesce(B1SE10F, BASC10F, RA1SF10F, RAASC11F),
    Nothing_Cheer = coalesce(B1SA24A, BACAS24A, RA1SA20A, RAACAS24A),
    Hopeless_Freq = coalesce(B1SA24D, BACAS24D, RA1SA20D, RAACAS24D),
    Everything_Effort = coalesce(B1SA24E, BACAS24E, RA1SA20E, RAACAS24E),
    Worthless_Freq = coalesce(B1SA24F, BACAS24F, RA1SA20F, RAACAS24F),
    Irritable_Freq = coalesce(B1SA24J, BACAS24J, RA1SA20J, RAACAS24J),
    Ashamed_Freq = coalesce(B1SA24K, BACAS24K, RA1SA20K, RAACAS24K),
    Emotional_Problems_Drinking = coalesce(B1SA66A, BACAS66A, RA1SA64A, RAACAS66B),
    Desire_Drinking = coalesce(B1SA66B, BACAS66B, RA1SA64B, RAACAS66C),
    Month_or_more_Drink = coalesce(B1SA66C, BACAS66C, RA1SA64C, RAACAS66D),
    Need_Drink_Same = coalesce(B1SA66D, BACAS66D, RA1SA64D, RAACAS66E),
    Alcohol_More_Than_Intended = coalesce(B1SA67, BACAS67, RA1SA65, RAACAS67),
    Alcohol_Affected_Work = coalesce(B1SA68, BACAS68, RA1SA66, RAACAS68), 
    Life_Day_By_Day = coalesce(B1SE1E, BASC1E, RA1SF1E, RAASC1E),
    Sense_Of_Direction = coalesce(B1SE1K, BASC1K, RA1SF1K, RAASC1K),
    No_Good_Sense = coalesce(B1SE1Q, BASC1Q, RA1SF1Q, RAASC1Q),
    Good_Manage_Responsibility = coalesce(B1SE1T, BASC1T, RA1SF1T, RAASC1T),
    Overwhelmed_Responsibility = coalesce(B1SE1Z, BASC1Z, RA1SF1Y, RAASC1Y),
    Actively_Carry_Plans = coalesce(B1SE1II, BASC1II, RA1SF1II, RAASC1II), 
    Some_Wander_Aimlessly = coalesce(B1SE1OO, BASC1OO, RA1SF1OO, RAASC1OO),
    Weigh_Possibilities = coalesce(B1SE7B,	BASC7B,	RA1SF7B,	RAASC8B),			
    Fun_Earthquake = coalesce(B1SE7D, BASC7D, RA1SF7D, RAASC8D), 
    Like_Think_Over = coalesce(B1SE7F, BASC7F, RA1SF7F, RAASC8F), 
    Keep_Work_Problems = coalesce(B1SE7L, BASC7L, RA1SF7L, RAASC8L), 
    Like_Difficult_Things = coalesce(B1SE7O, BASC7O, RA1SF7O, RAASC8O), 
    Like_Hard_Work = coalesce(B1SE7R, BASC7R, RA1SF7R, RAASC8R), 
    Cautious_Person = coalesce(B1SE7Y, BASC7Y, RA1SF7Y, RAASC8Y), 
    Fun_Tight_Rope = coalesce(B1SE7V, BASC7V, RA1SF7V, RAASC8V), 
    High_Standards = coalesce(B1SE7FF, BASC7FF, RA1SF7FF, RAASC8FF), 
    Plan_Future = coalesce(B1SE12O, BASC11O, RA1SF14O, RAASC12O), 
    Live_Day_By_Day = coalesce(B1SE1E, BASC1E, RA1SF1E, RAASC1E),
    Enjoy_Plans_Future = coalesce(B1SE1CC, BASC1CC, RA1SF1CC, RAASC1CC),
    Difficult_Arrange = coalesce(B1SE1FF,	BASC1FF,	RA1SF1FF,	RAASC1FF),
    Maintaining_Relationships_Hard = coalesce(B1SE1J, BASC1J, RA1SF1J, RAASC1J),
    Dont_Fit_Community = coalesce(B1SE1N, BASC1N, RA1SF1N, RAASC1N),
    Few_Close_Friends = coalesce(B1SE1P, BASC1P, RA1SF1P, RAASC1P),
    Enjoy_Convo = coalesce(B1SE1V, BASC1V, RA1SF1V, RAASC1V),
    No_Warm_Relationship = coalesce(B1SE1HH, BASC1HH, RA1SF1HH, RAASC1HH),
    Like_Time_Alone = coalesce(B1SE7A, BASC7A, RA1SF7A, RAASC8A), 
    Seek_Friends = coalesce(B1SE7C, BASC7C, RA1SF7C, RAASC8C), 
    I_Am_Warm_Person = coalesce(B1SE7H, BASC7H, RA1SF7H, RAASC8H), 
    Prefer_No_Others = coalesce(B1SE7CC, BASC7CC, RA1SF7CC, RAASC8CC), 
    Felt_Close_Others_Freq = coalesce(B1SA26G, BACAS26G, RA1SA22G, RAACAS26G),
    Felt_Belong_Freq = coalesce(B1SA26H, BACAS26H, RA1SA22H, RAACAS26H),
    Trust_Friends = coalesce(B1SE1NN, BASC1NN, RA1SF1NN, RAASC1NN),
    Most_See_Love = coalesce(B1SE1D, BASC1D, RA1SF1D, RAASC1D),
    Enthusiastic_Freq = coalesce(B1SA26I, BACAS26I, RA1SA22I, RAACAS26I),
    Full_Of_Life_Freq = coalesce(B1SA26F, BACAS26F, RA1SA22F, RAACAS26F),
    Extreme_Happy_Freq = coalesce(B1SA26C, BACAS26C, RA1SA22C, RAACAS26C),
    Good_Spirits_Freq = coalesce(B1SA26B, BACAS26B, RA1SA22B, RAACAS26B),
    Cheerful_Freq = coalesce(B1SA26A, BACAS26A, RA1SA22A, RAACAS26A),
    More_Or_Less_Positive = coalesce(B1SA27, BACAS27, RA1SA23, RAACAS27), 
    Lonely_Frequency = coalesce(B1SA24G, BACAS24G, RA1SA20G, RAACAS24G), 
    Angry_Freq = coalesce(B1SA24M, BACAS24M, RA1SA20M, RAACAS24M),
    Frustrated_Freq = coalesce(B1SA24N, BACAS24N, RA1SA20N, RAACAS24N),
    Not_Afraid_Opinion = coalesce(B1SE1A, BASC1A, RA1SF1A, RAASC1A), 
    Confidence_If_Contrary = coalesce(B1SE1S, BASC1S, RA1SF1S, RAASC1S),
    Difficult_Voice_Opinion = coalesce(B1SE1Y, BASC1Y, RA1SF1Y, RAASC1Y),
    Others_Describe_Giving = coalesce(B1SE1BB, BASC1BB, RA1SF1BB, RAASC1BB),
    Angry_Ready_Hit = coalesce(B1SE7I, BASC7I, RA1SF7I, RAASC8I), 
    Effective_Talking_People = coalesce(B1SE7J, BASC7J, RA1SF7J, RAASC8J), 
    Good_Influence_People = coalesce(B1SE7N, BASC7N, RA1SF7N, RAASC8N),
    Enjoy_Hurting_Say_Mean = coalesce(B1SE7T, BASC7T, RA1SF7T, RAASC8T),
    Should_Obey_Law = coalesce(B1SE7U, BASC7U, RA1SF7U, RAASC8U),
    When_Insult_Get_Even = coalesce(B1SE7BB, BASC7BB, RA1SF7BB, RAASC8BB),
    Sometimes_Just_Hit_Someone = coalesce(B1SE7EE, BASC7EE, RA1SF7EE, RAASC8EE),
    Important_Help = coalesce(B1SH16P, BASF1P, RA1SI16P, RAASG1P),
    Not_Happy_If_Friend_Trouble = coalesce(B1SH16Q, BASF1Q, RA1SI16Q, RAASG1Q),
    Moved_By_Other_Hardship = coalesce(B1SH16R, BASF1R, RA1SI16R, RAASG1R),
    Important_Sympathetic = coalesce(B1SH16S, BASF1S, RA1SI16S, RAASG1S),
    Sympathy_Limits = coalesce(B1SH16T, BASF1T, RA1SI16T, RAASG1T), 
    Positive_About_Self = coalesce(B1SE1L, BASC1L, RA1SF1L, RAASC1L),
    Mental_Emotional = coalesce(B1PA2,	BACA2,	RA1PA2,	RAACA2),
    Anhedonia = coalesce(B1PANHED, BACANHED, RA1PANHED, RAACANHED),
    Depressed_Affect = coalesce(B1PDEPAF, BACDEPAF, RA1PDEPAF, RAACDEPAF),
    Activities_Trivial = coalesce(B1SE1W, BASC1W,	RA1SF1W, RAASC1W),
    Self_Attitude = coalesce(B1SE1JJ, BASC1JJ, RA1SF1JJ, RAASC1JJ),
    Approach_Health = coalesce(B1SE5E,	BASC5E,	RA1SF5E,	RAASC6E),
    Choose_Goals = coalesce(B1SE5A,	BASC5A,	RA1SF5A,	RAASC6A),
    Will = coalesce(B1SE12A,	BASC11A,	RA1SF14A,	RAASC12A),
    Change_For_Better = coalesce(B1SE12B,	BASC11B,	RA1SF14B,	RAASC12B),
    Lower_Expectations = coalesce(B1SE12C,	BASC11C, RA1SF14C,	RAASC12C),
    Dissapointment_Low_Goal = coalesce(B1SE12D,	BASC11D,	RA1SF14D,	RAASC12D),
    Releived_No_Responsibility = coalesce(B1SE12F,	BASC11F,	RA1SF14F,	RAASC12F),
    Too_Much_Get_Done = coalesce(B1SE12G,	BASC11G,	RA1SF14G,	RAASC12G),
    Dont_Give_Up = coalesce(B1SE12J,	BASC11J,	RA1SF14J,	RAASC12J),
    Rarely_Give_Up = coalesce(B1SE12K,	BASC11K,	RA1SF14K,	RAASC12K),
    Know_Want_Life = coalesce(B1SE12P,	BASC11P,	RA1SF14P,	RAASC12P),
    Live_One_Day = coalesce(B1SE12Q,	BASC11Q,	RA1SF14Q,	RAASC12Q),
    Helpful_Set_Future_Goal = coalesce(B1SE12R,	BASC11R,	RA1SF14R,	RAASC12R),
    Too_Many_Things = coalesce(B1SE12S,	BASC11S,	RA1SF14S,	RAASC12S),
    Bad_Happens_Prevented = coalesce(B1SE12V,	BASC11V,	RA1SF14V,	RAASC12V),
    No_Sense_Planning = coalesce(B1SE12X,	BASC11X,	RA1SF14X,	RAASC12X),
    Goal_Keep_Benefits = coalesce(B1SE12AA,	BASC11AA,	RA1SF14AA,	RAASC12AA),
    Avoid_Distraction = coalesce(B1SE12EE,	BASC11EE,	RA1SF14EE,	RAASC12EE),
    If_No_Goal_Think_Others = coalesce(B1SE12KK,	BASC11KK,	RA1SF14KK,	RAASC12KK),
    To_Reach_Goals = coalesce(B1SE5B,	BASC5B,	RA1SF5B,	RAASC6B),
    Without_skill = coalesce(B1SE5C,	BASC5C,	RA1SF5C,	RAASC6C),
    When_Difficult = coalesce(B1SE5D,	BASC5D,	RA1SF5D,	RAASC6D),
    More_Sucessful = coalesce(B1SE7P,	BASC7P,	RA1SF7P,	RAASC8P),
    Influenced_By_People = coalesce(B1SE1M,	BASC1M,	RA1SF1M,	RAASC1M),
    Social_Others_Lead = coalesce(B1SE7E,	BASC7E,	RA1SF7E,	RAASC8E),
    Dont_Like_Help = coalesce(B1SE12Y,	BASC11Y,	RA1SF14Y,	RAASC12Y),
    Asking_Help_Naturally = coalesce(B1SE12Z,	BASC11Z,	RA1SF14Z,	RAASC12Z),
    Dont_Feel_Community = coalesce(B1SH16B,	BASF1B,	RA1SI16B,	RAASG1B),
    Feel_Close_Others = coalesce(B1SH16F,	BASF1F,	RA1SI16F,	RAASG1F),
    Community_Source_Comfort = coalesce(B1SH16K,	BASF1K,	RA1SI16K,	RAASG1K),
    People_Kind = coalesce(B1SH16N,	BASF1N,	RA1SI16N,	RAASG1N),
    People_Take_Advantage = coalesce(B1SE7G,	BASC7G,	RA1SF7G,	RAASC8G),
    Others_Turn_To_Me = coalesce(B1SE7DD,	BASC7DD,	RA1SF7DD,	RAASC8DD),
    Status_SAQ = coalesce(B1STATUS,	BASTATUS,	RA1STATUS,	RAASTATUS),
    .keep = "none")

# Apply to all factor columns in the data frame. This function will 
# remove the characters and symbols in each factor, except for the number and factor 
# structure. This is necessary because most factors are coded separately
# across waves (but have the same factor structure), and are spelled incorrectly!!
# be cautious :)

data <- subset %>%
  dplyr::mutate(across(where(is.factor), ~ 
                         factor(gsub("^\\((\\d+)\\).*", "\\1", as.character(.)))
  ))

#### Date Manipulation. Proceed with caution! Perilous mistakes are easy to make ####

#Include only those that have completed the SAQ. This only includes individuals who scored as 2: Interview and SAQ data, or 4: Interview and SAQ data and biometrics. 

data <- data %>%
  filter(Status_SAQ %in% c("2", "4")) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

# Reverse coding variables for interpretation

rc_apply <- function(data, vars, min_val, max_val) {
  v <- intersect(names(data), vars)
  if (length(v)) data[v] <- lapply(data[v], function(x) ifelse(is.na(x), NA, (min_val + max_val) - x))
  data
}

# --------- Groups by scale  ---------

# 1–3
rc_1_3 <- c(
  "Worry_Amount" # 1=a lot more, 3=a little (reverse)
)

# 1–4
rc_1_4 <- c(
  "Worry_Compare","People_Mean","Upset_Think_Day","Minor_Setback_Irritate",
  "Mood_Up_Down","Exciting_Look_Forward",
  "Weigh_Possibilities","Fun_Earthquake","Like_Think_Over","Keep_Work_Problems",
  "Like_Difficult_Things","Like_Hard_Work","Cautious_Person","Fun_Tight_Rope",
  "High_Standards","Plan_Future",
  "Will","Change_For_Better","Lower_Expectations","Dissapointment_Low_Goal","Releived_No_Responsibility",
  "Too_Much_Get_Done","Dont_Give_Up","Rarely_Give_Up","Know_Want_Life","Live_One_Day",
  "Helpful_Set_Future_Goal","Too_Many_Things","Bad_Happens_Prevented","No_Sense_Planning",
  "Goal_Keep_Benefits","Avoid_Distraction","If_No_Goal_Think_Others",
  "Like_Time_Alone","Seek_Friends","Social_Others_Lead","People_Take_Advantage",
  "I_Am_Warm_Person","Prefer_No_Others","Others_Turn_To_Me","Dont_Like_Help","Asking_Help_Naturally",
  "Angry_Ready_Hit","Effective_Talking_People","Good_Influence_People","Enjoy_Hurting_Say_Mean",
  "Should_Obey_Law","When_Insult_Get_Even","Sometimes_Just_Hit_Someone","Outgoing","Helpful", 
  "Moody","Organized","Self_Confident","Friendly","Warm","Worrying","Responsible","Forceful",
  "Lively","Caring","Nervous","Creative","Assertive","Hardworking","Imaginative","Softhearted",
  "Calm","Outspoken","Intelligent","Curious","Active","Careless","Broadminded","Sympathetic", 
  "Talkative","Sophisticated","Adventurous","Dominant"
)

# 1–5
rc_1_5 <- c(
  "Worry_Frequency",
  "Nervous_Freq_30_Day","Restless_Freq_30_Day","Afraid_Freq_30_Day","Jittery_Freq_30_Day",
  "Irritable_Freq_30_Day","Frustrated_Freq_30_Day","Upset_Freq_30_Day",
  "Nothing_Cheer","Hopeless_Freq","Everything_Effort","Worthless_Freq","Irritable_Freq","Ashamed_Freq",
  "Enthusiastic_Freq","Full_Of_Life_Freq","Extreme_Happy_Freq","Good_Spirits_Freq","Cheerful_Freq",
  "Lonely_Frequency","Felt_Close_Others_Freq","Felt_Belong_Freq",
  "Mental_Emotional" # 1=excellent, 5=poor (reverse)
)

# 1–6
rc_1_6 <- c(
  "Headache_Frequency","Backache_Frequency","Swea_Frequency","Irritabiltiy_Frequency",
  "Hot_Flashes","Aches_Joints_Stiffness_Frequency","Falling_Asleep","Leaking_Urine",
  "Intercourse_Pain","Extremities_Aches"
)

# 1–7
rc_1_7 <- c(
  "Positive_About_Self","Others_More_Life","Activities_Trivial","Disapoint_Acheivment","Gave_Up","Self_Attitude",
  "Expect_The_Best","If_Go_Wrong_It_Will","Optimistic_About_Future","Hardly_Expect_Go_My_Way",
  "Rarely_Good_Happen","Expect_Good_Not_Bad","More_Or_Less_Positive",
  "Life_Day_By_Day","Sense_Of_Direction","No_Good_Sense","Good_Manage_Responsibility",
  "Overwhelmed_Responsibility","Actively_Carry_Plans","Some_Wander_Aimlessly","Enjoy_Plans_Future",
  "Difficult_Arrange","Live_Day_By_Day",
  "Maintaining_Relationships_Hard","Influenced_By_People","Dont_Fit_Community","Few_Close_Friends",
  "Enjoy_Convo","No_Warm_Relationship","Dont_Feel_Community","Feel_Close_Others",
  "Community_Source_Comfort","People_Kind","Trust_Friends","Most_See_Love",
  "Not_Afraid_Opinion","Confidence_If_Contrary","Difficult_Voice_Opinion","Others_Describe_Giving",
  "Important_Help","Not_Happy_If_Friend_Trouble","Moved_By_Other_Hardship","Important_Sympathetic","Sympathy_Limits"
)

# Binary 1–2
rc_bin_1_2 <- c(
  "Emotional_Problems_Drinking","Desire_Drinking","Month_or_more_Drink","Need_Drink_Same"
)

# --------- Apply reversals  ---------
data <- rc_apply(data, rc_1_3, 1, 3)
data <- rc_apply(data, rc_1_4, 1, 4)
data <- rc_apply(data, rc_1_5, 1, 5)
data <- rc_apply(data, rc_1_6, 1, 6)
data <- rc_apply(data, rc_1_7, 1, 7)
data <- rc_apply(data, rc_bin_1_2, 1, 2)


#Multiple Imputation
set.seed(123)
init <- mice(data, maxit = 0)
pred_matrix <- init$predictorMatrix
pred_matrix[, c("ID", "Status_SAQ")] <- 0   # don't use as predictors

meth <- init$method
meth[] <- "pmm"                              # PMM for everything
meth[c("ID", "Status_SAQ")] <- ""            # don't impute these two


imputed_data <- mice(
  data,
  m = 20,
  method = meth,  # predictive mean matching, good for mixed data types
  predictorMatrix = pred_matrix,
  seed = 123
)

completed_data_for_comparison <- complete(imputed_data, action = "long", include = TRUE)

write.csv(completed_data_for_comparison, "completed_data_for_comparison_final_2.csv", row.names = FALSE)

completed_data_for_analysis <- complete(imputed_data, action = "long", include = FALSE)

write.csv(completed_data_for_analysis, "completed_data_for_analysis_final_2.csv", row.names = FALSE)






