# ==============================================================
# Measurement Invariance
# Groups tested: Sex and Race
# ==============================================================

# --- Packages ---
if (!requireNamespace("lavaan", quietly=TRUE)) install.packages("lavaan", repos="https://cloud.r-project.org")
library(lavaan)

# --- Input and output ---
data_path  <- "completed_data_for_analysis_final.csv"
output_dir <- "MI_Results_Simple_MIN"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Model ---
cfa_model <- '
  Internalizing =~ Worry_Frequency + Nervous_Freq_30_Day + Restless_Freq_30_Day + 
                   Afraid_Freq_30_Day + Frustrated_Freq_30_Day + Upset_Freq_30_Day + 
                   Hopeless_Freq + Everything_Effort + Worthless_Freq + Irritable_Freq + 
                   Ashamed_Freq + Lonely_Frequency + Angry_Freq

  Anatagonistic =~ People_Mean + Upset_Think_Day + Minor_Setback_Irritate + 
                   Mood_Up_Down + Angry_Ready_Hit + Enjoy_Hurting_Say_Mean + 
                   When_Insult_Get_Even + Sometimes_Just_Hit_Someone + More_Sucessful

  Detachment =~ Others_More_Life + Disapoint_Acheivment + Gave_Up + No_Good_Sense + 
                Difficult_Arrange + Maintaining_Relationships_Hard + Dont_Fit_Community + 
                Few_Close_Friends + No_Warm_Relationship + Activities_Trivial + Self_Attitude

  Somaticizing =~ Backache_Frequency + Swea_Frequency + Hot_Flashes + 
                  Falling_Asleep + Extremities_Aches

  Disinhibited_Externalizing =~ Actively_Carry_Plans + Like_Difficult_Things + 
                  Like_Hard_Work + High_Standards + Plan_Future + Will + Change_For_Better + 
                  Too_Much_Get_Done + Dont_Give_Up + Know_Want_Life + 
                  Goal_Keep_Benefits + Weigh_Possibilities
'

# --- Load + pick ONE imputation  --- #Couldn't figure out MI on multiple imputations...
dat <- read.csv(data_path, check.names = FALSE)

target_imp <- 1
if (".imp" %in% names(dat)) {
  if (!target_imp %in% unique(dat$.imp)) stop("Requested .imp not found.")
  dat <- dat[dat$.imp == target_imp, , drop = FALSE]
}

# --- Items from model ---
pt   <- lavaan::lavaanify(cfa_model, fixed.x = FALSE)
lvs  <- unique(pt$lhs[pt$op == "=~"])
items <- unique(pt$rhs[pt$op == "=~" & !(pt$rhs %in% lvs)])

# Keep only items that exist
items <- intersect(items, names(dat))
if (!length(items)) stop("No model items found in data.")

# Make groups factors if present; drop rows with missing group per test
if ("Sex"  %in% names(dat))  dat$Sex  <- factor(dat$Sex)
if ("Race" %in% names(dat))  dat$Race <- factor(dat$Race)

# Coerce items to ordered (levels from observed values)
to_ord <- function(x) { xi <- as.integer(suppressWarnings(as.numeric(x))); ordered(xi, levels = sort(unique(xi[!is.na(xi)]))) }
for (nm in items) dat[[nm]] <- to_ord(dat[[nm]])

# Reverse code necessary items
rev_items <- intersect(c("Angry_Freq","More_Sucessful"), items)
if (length(rev_items)) {
  rev_ord <- function(x){ xi <- as.numeric(as.character(x)); lo <- min(xi, na.rm=TRUE); hi <- max(xi, na.rm=TRUE); ordered(lo + hi - xi, levels = sort(unique(lo + hi - xi))) }
  for (nm in rev_items) dat[[nm]] <- rev_ord(dat[[nm]])
}

# Drop items that are constant within ANY group level
vary_in_groups <- function(d, itms, gvar) {
  if (!gvar %in% names(d)) return(itms)
  glv <- levels(d[[gvar]])[table(d[[gvar]]) > 0]
  keep <- vapply(itms, function(v) all(sapply(glv, function(g) {
    x <- d[d[[gvar]] == g, v]
    length(unique(na.omit(x))) >= 2
  })), logical(1))
  itms[keep]
}

# Fit 3-step ordinal invariance and save a tiny fit table
run_simple <- function(gvar, label) {
  if (!gvar %in% names(dat)) return(invisible(NULL))
  d <- dat[!is.na(dat[[gvar]]), , drop = FALSE]
  if (nlevels(d[[gvar]]) < 2) return(invisible(NULL))

  itms <- vary_in_groups(d, items, gvar)
  if (!length(itms)) {
    write.csv(data.frame(Message="No varying items within group levels."), file.path(output_dir, paste0(label,"_Fit.csv")), row.names=FALSE)
    return(invisible(NULL))
  }

  fm <- function(fit) {
    want.s <- c("chisq.scaled","df.scaled","pvalue.scaled","cfi.scaled","tli.scaled","rmsea.scaled","srmr")
    out <- try(fitMeasures(fit, want.s), silent=TRUE)
    if (inherits(out,"try-error")) out <- fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
    as.list(out)
  }

  baseArgs <- list(model=cfa_model, data=d, group=gvar,
                   estimator="WLSMV", parameterization="theta",
                   ordered=itms, std.lv=TRUE)

  fit_cfg <- try(do.call(cfa, c(baseArgs, list(group.equal=NULL))), silent=TRUE)
  fit_thr <- try(do.call(cfa, c(baseArgs, list(group.equal="thresholds"))), silent=TRUE)
  fit_met <- try(do.call(cfa, c(baseArgs, list(group.equal=c("thresholds","loadings")))), silent=TRUE)

  rows <- list()
  if (inherits(fit_cfg, "lavaan")) rows[["configural"]]  <- fm(fit_cfg)
  if (inherits(fit_thr, "lavaan")) rows[["thresholds"]]  <- fm(fit_thr)
  if (inherits(fit_met, "lavaan")) rows[["metric"]]      <- fm(fit_met)

  if (!length(rows)) {
    write.csv(data.frame(Message="No models converged."), file.path(output_dir, paste0(label,"_Fit.csv")), row.names=FALSE)
    return(invisible(NULL))
  }

  out <- do.call(rbind, lapply(names(rows), function(nm) {
    c(Model=nm, rows[[nm]][c("chisq.scaled","df.scaled","pvalue.scaled","cfi.scaled","tli.scaled","rmsea.scaled","srmr")])
  }))
  write.csv(out, file.path(output_dir, paste0(label,"_Fit.csv")), row.names=FALSE)

  # Quick deltas if both stages exist
  if (all(c("configural","thresholds") %in% names(rows))) {
    a <- rows[["configural"]]; b <- rows[["thresholds"]]
    d1 <- data.frame(Comparison="configural vs thresholds",
      Delta_CFI = round(as.numeric(b[["cfi.scaled"]]) - as.numeric(a[["cfi.scaled"]]), 3),
      Delta_RMSEA = round(as.numeric(b[["rmsea.scaled"]]) - as.numeric(a[["rmsea.scaled"]]), 3),
      Delta_SRMR = round(as.numeric(b[["srmr"]]) - as.numeric(a[["srmr"]]), 3))
    write.csv(d1, file.path(output_dir, paste0(label,"_Deltas.csv")), row.names=FALSE)
  }
  if (all(c("thresholds","metric") %in% names(rows))) {
    a <- rows[["thresholds"]]; b <- rows[["metric"]]
    d2 <- data.frame(Comparison="thresholds vs metric",
      Delta_CFI = round(as.numeric(b[["cfi.scaled"]]) - as.numeric(a[["cfi.scaled"]]), 3),
      Delta_RMSEA = round(as.numeric(b[["rmsea.scaled"]]) - as.numeric(a[["rmsea.scaled"]]), 3),
      Delta_SRMR = round(as.numeric(b[["srmr"]]) - as.numeric(a[["srmr"]]), 3))
    fn <- file.path(output_dir, paste0(label,"_Deltas.csv"))
    if (file.exists(fn)) write.table(d2, fn, sep=",", row.names=FALSE, col.names=FALSE, append=TRUE) else write.csv(d2, fn, row.names=FALSE)
  }

  invisible(NULL)
}

# --- Run ---
if ("Sex"  %in% names(dat))  run_simple("Sex",  "Sex")
if ("Race" %in% names(dat))  run_simple("Race", "Race")

cat("\nDone. CSVs in: ", normalizePath(output_dir), "\n", sep="")
