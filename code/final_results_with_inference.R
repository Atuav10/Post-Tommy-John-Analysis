#Final Results Script with Inference Metrics

library(dplyr)
library(tidyr)
library(tibble)
library(synthdid)

#--- 1. Load and Prep Data ----------------------------------------
all_pitchers <- read.csv("pitchers_updated.csv")
fip_constants <- read.csv("fip_constant.csv")
all_pitchers <- left_join(all_pitchers, fip_constants, by = c("season_code" = "season"))
all_pitchers$FIP <- (13*all_pitchers$HR + 3*(all_pitchers$BB + all_pitchers$HBP) - 2*all_pitchers$SO)/all_pitchers$IP + all_pitchers$fip_constant

all_pitchers <- select(all_pitchers, bbref_id, Name, Age, season_code, Team, G, GS, IP, ERA, BAbip, WHIP, SO9, LD, PU, GB.FB, StS, StL, SO_perc, FIP)
all_pitchers <- filter(all_pitchers, G >= 20 | GS > 10)
all_pitchers <- filter(all_pitchers, season_code != 2020)

player_surgery_df <- read.csv("treated_players.csv", stringsAsFactors = FALSE)

#Get list of all TJ surgery pitcher names to exclude from donor pools
all_tj_pitchers <- unique(player_surgery_df$player_name)

#--- 2. Initialize Results Data Frame ----------------------------------------
results_df <- data.frame(
  player_name = character(),
  surgery_year = numeric(),
  number_of_donors = numeric(),

  ate_so9_scm = numeric(),
  so9_year1 = numeric(),
  so9_year2 = numeric(),
  so9_year3 = numeric(),
  ate_so9_did = numeric(),
  ate_so9_synthdid = numeric(),
  ate_percentile_so9 = numeric(),
  ate_variance_so9 = numeric(),

  ate_fip_scm = numeric(),
  fip_year1 = numeric(),
  fip_year2 = numeric(),
  fip_year3 = numeric(),
  ate_fip_did = numeric(),
  ate_fip_synthdid = numeric(),
  ate_percentile_fip = numeric(),
  ate_variance_fip = numeric(),

  stringsAsFactors = FALSE
)

#Initialize Donor Weights Data Frame
donor_weights_df <- data.frame(
  player_name = character(),
  donor_name = character(),
  weight_so9 = numeric(),
  weight_fip = numeric(),
  stringsAsFactors = FALSE
)

#--- 3. Function to Run Synthetic Control with All Estimators ----------------------------------------
run_all_estimators <- function(player_name, injury_year, outcome_var, donor_stats_df) {
  tryCatch({
    #Prepare data
    sdid_data <- donor_stats_df %>%
      select(Name, Year = season_code, outcome = !!sym(outcome_var)) %>%
      mutate(treated = ifelse(Name == player_name & Year > injury_year, 1, 0))

    panel_long <- sdid_data %>%
      select(Name, Year, outcome, treated)

    setup <- panel.matrices(panel_long)

    if (setup$N0 == 0 || setup$T0 == 0) {
      stop(sprintf("Invalid setup: N0=%d, T0=%d", setup$N0, setup$T0))
    }

    #Run all three estimators
    tau.sc <- sc_estimate(setup$Y, setup$N0, setup$T0)
    tau.sdid <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
    tau.did <- did_estimate(setup$Y, setup$N0, setup$T0)

    #Calculate year-by-year effects (SC only)
    omega_weights <- attr(tau.sc, 'weights')$omega
    synthetic_control <- t(omega_weights) %*% setup$Y[1:setup$N0, ]
    synthetic_trajectory <- synthetic_control[(setup$T0 + 1):ncol(setup$Y)]
    treated_trajectory <- setup$Y[setup$N0 + 1, (setup$T0 + 1):ncol(setup$Y)]
    yearly_effects <- treated_trajectory - synthetic_trajectory

    year1 <- if(length(yearly_effects) >= 1) yearly_effects[1] else NA
    year2 <- if(length(yearly_effects) >= 2) yearly_effects[2] else NA
    year3 <- if(length(yearly_effects) >= 3) yearly_effects[3] else NA

    #Get donor names and weights
    donor_names <- rownames(setup$Y)[1:setup$N0]
    weights_vector <- as.numeric(omega_weights)

    return(list(
      ate_scm = as.numeric(tau.sc),
      year1 = year1,
      year2 = year2,
      year3 = year3,
      ate_synthdid = as.numeric(tau.sdid),
      ate_did = as.numeric(tau.did),
      donor_names = donor_names,
      weights = weights_vector,
      success = TRUE
    ))

  }, error = function(e) {
    return(list(
      ate_scm = NA,
      year1 = NA,
      year2 = NA,
      year3 = NA,
      ate_synthdid = NA,
      ate_did = NA,
      donor_names = NULL,
      weights = NULL,
      success = FALSE
    ))
  })
}

#--- 4. Function to Run Placebo Tests ----------------------------------------
run_placebo_tests <- function(player_name, injury_year, years_needed, donor_ids, outcome_var, all_tj_pitchers) {

  donor_names_list <- all_pitchers %>%
    filter(bbref_id %in% donor_ids) %>%
    select(Name) %>%
    distinct() %>%
    pull(Name)

  placebo_ates <- c()

  for (placebo_player in donor_names_list) {
    #Create donor pool excluding this placebo player (and ensure no TJ pitchers)
    placebo_donor_stats <- all_pitchers %>%
      filter(season_code %in% years_needed) %>%
      filter(!Name %in% all_tj_pitchers) %>%
      filter(Name == placebo_player | (bbref_id %in% donor_ids & Name != placebo_player)) %>%
      select(bbref_id, Name, season_code, Age, SO9, FIP)

    #Run SC
    result <- tryCatch({
      sdid_data <- placebo_donor_stats %>%
        select(Name, Year = season_code, outcome = !!sym(outcome_var)) %>%
        mutate(treated = ifelse(Name == placebo_player & Year > injury_year, 1, 0))

      panel_long <- sdid_data %>%
        select(Name, Year, outcome, treated)

      setup <- panel.matrices(panel_long)
      tau.sc <- sc_estimate(setup$Y, setup$N0, setup$T0)
      as.numeric(tau.sc)
    }, error = function(e) {
      NA_real_
    })

    if (!is.na(result)) {
      placebo_ates <- c(placebo_ates, result)
    }
  }

  return(placebo_ates)
}

#--- 5. Main Loop Through Players ----------------------------------------
cat("========================================\n")
cat("RUNNING COMPREHENSIVE ANALYSIS\n")
cat("========================================\n\n")

for (i in 1:nrow(player_surgery_df)) {

  player_name <- player_surgery_df$player_name[i]
  injury_year <- player_surgery_df$surgery_year[i]

  cat(sprintf("\n[%d/%d] Processing %s (surgery year: %d)...\n", i, nrow(player_surgery_df), player_name, injury_year))

  tryCatch({
    #Filter to player's data
    player_data <- filter(all_pitchers, Name == player_name)

    if (nrow(player_data) == 0) {
      cat(sprintf("  ✗ No data found\n"))
      next
    }

    #Define pre/post treatment periods
    player_pre <- player_data %>%
      filter(season_code <= injury_year) %>%
      arrange(desc(season_code)) %>%
      head(3)

    player_post <- player_data %>%
      filter(season_code > injury_year) %>%
      arrange(season_code) %>%
      head(3)

    if (nrow(player_pre) < 2 || nrow(player_post) < 1) {
      cat(sprintf("  ✗ Insufficient data (pre=%d, post=%d)\n", nrow(player_pre), nrow(player_post)))
      next
    }

    #Determine required seasons
    years_needed <- unique(c(player_pre$season_code, player_post$season_code))

    #Identify donor pool (exclude ALL TJ surgery pitchers)
    complete_pitchers <- all_pitchers %>%
      filter(!Name %in% all_tj_pitchers) %>%
      filter(season_code %in% years_needed) %>%
      group_by(Name, bbref_id) %>%
      summarize(n_seasons = n_distinct(season_code), .groups = "drop") %>%
      filter(n_seasons == length(years_needed))

    donor_ids <- complete_pitchers$bbref_id

    if (length(donor_ids) < 5) {
      cat(sprintf("  ✗ Donor pool too small (%d donors)\n", length(donor_ids)))
      next
    }

    #Create combined dataset
    donor_stats_df <- all_pitchers %>%
      filter(bbref_id %in% donor_ids | Name == player_name) %>%
      filter(season_code %in% years_needed) %>%
      select(bbref_id, Name, season_code, Age, SO9, FIP)

    cat(sprintf("  Donor pool: %d pitchers | Years: %s\n",
                length(donor_ids), paste(years_needed, collapse=", ")))

    #Run estimators for SO9
    cat("  Running SO9 estimators...\n")
    so9_results <- run_all_estimators(player_name, injury_year, "SO9", donor_stats_df)

    #Run estimators for FIP
    cat("  Running FIP estimators...\n")
    fip_results <- run_all_estimators(player_name, injury_year, "FIP", donor_stats_df)

    if (!so9_results$success && !fip_results$success) {
      cat(sprintf("  ✗ Both SO9 and FIP failed\n"))
      next
    }

    #Run placebo tests for SO9
    cat("  Running SO9 placebo tests...\n")
    so9_placebo_ates <- run_placebo_tests(player_name, injury_year, years_needed, donor_ids, "SO9", all_tj_pitchers)

    #Calculate SO9 percentile and variance
    if (length(so9_placebo_ates) > 0 && !is.na(so9_results$ate_scm)) {
      so9_percentile <- sum(abs(so9_placebo_ates) <= abs(so9_results$ate_scm)) / length(so9_placebo_ates)
      so9_variance <- so9_results$ate_scm / sd(so9_placebo_ates, na.rm = TRUE)
    } else {
      so9_percentile <- NA
      so9_variance <- NA
    }

    #Run placebo tests for FIP
    cat("  Running FIP placebo tests...\n")
    fip_placebo_ates <- run_placebo_tests(player_name, injury_year, years_needed, donor_ids, "FIP", all_tj_pitchers)

    #Calculate FIP percentile and variance
    if (length(fip_placebo_ates) > 0 && !is.na(fip_results$ate_scm)) {
      fip_percentile <- sum(abs(fip_placebo_ates) <= abs(fip_results$ate_scm)) / length(fip_placebo_ates)
      fip_variance <- fip_results$ate_scm / sd(fip_placebo_ates, na.rm = TRUE)
    } else {
      fip_percentile <- NA
      fip_variance <- NA
    }

    #Add to results
    results_df <- rbind(results_df, data.frame(
      player_name = player_name,
      surgery_year = injury_year,
      number_of_donors = length(donor_ids),

      ate_so9_scm = so9_results$ate_scm,
      so9_year1 = so9_results$year1,
      so9_year2 = so9_results$year2,
      so9_year3 = so9_results$year3,
      ate_so9_did = so9_results$ate_did,
      ate_so9_synthdid = so9_results$ate_synthdid,
      ate_percentile_so9 = so9_percentile,
      ate_variance_so9 = so9_variance,

      ate_fip_scm = fip_results$ate_scm,
      fip_year1 = fip_results$year1,
      fip_year2 = fip_results$year2,
      fip_year3 = fip_results$year3,
      ate_fip_did = fip_results$ate_did,
      ate_fip_synthdid = fip_results$ate_synthdid,
      ate_percentile_fip = fip_percentile,
      ate_variance_fip = fip_variance,

      stringsAsFactors = FALSE
    ))

    #Add donor weights - merge SO9 and FIP weights by donor name
    if (!is.null(so9_results$donor_names) && !is.null(fip_results$donor_names)) {
      #Get unique list of all donors from both SO9 and FIP
      all_donor_names <- unique(c(so9_results$donor_names, fip_results$donor_names))

      for (donor in all_donor_names) {
        #Find weight in SO9
        so9_idx <- which(so9_results$donor_names == donor)
        weight_so9 <- if(length(so9_idx) > 0) so9_results$weights[so9_idx] else 0

        #Find weight in FIP
        fip_idx <- which(fip_results$donor_names == donor)
        weight_fip <- if(length(fip_idx) > 0) fip_results$weights[fip_idx] else 0

        donor_weights_df <- rbind(donor_weights_df, data.frame(
          player_name = player_name,
          donor_name = donor,
          weight_so9 = weight_so9,
          weight_fip = weight_fip,
          stringsAsFactors = FALSE
        ))
      }
    }

    cat(sprintf("  ✓ Success! SO9 ATE=%.2f (p=%.3f) | FIP ATE=%.2f (p=%.3f)\n",
                so9_results$ate_scm, 1-so9_percentile,
                fip_results$ate_scm, 1-fip_percentile))

  }, error = function(e) {
    cat(sprintf("  ✗ Error: %s\n", e$message))
  })
}

#--- 6. Save Results ----------------------------------------
cat("\n========================================\n")
cat("FINAL RESULTS\n")
cat("========================================\n")
cat(sprintf("Successfully processed: %d / %d players\n\n", nrow(results_df), nrow(player_surgery_df)))
write.csv(results_df, "final_results_with_inference.csv", row.names = FALSE)
cat("Results saved to: final_results_with_inference.csv\n")

donor_weights_df <- filter(donor_weights_df, weight_so9 > 0 | weight_fip > 0)
write.csv(donor_weights_df, "donor_weights.csv", row.names = FALSE)
cat("Donor weights saved to: donor_weights.csv\n")
cat(sprintf("Total donor-player combinations: %d\n\n", nrow(donor_weights_df)))

#Display summary statistics
if (nrow(results_df) > 0) {
  cat("Summary Statistics:\n")
  cat("\nSO9 Results:\n")
  cat(sprintf("  ATE SCM - Mean: %.2f, Median: %.2f, SD: %.2f\n",
              mean(results_df$ate_so9_scm, na.rm = TRUE),
              median(results_df$ate_so9_scm, na.rm = TRUE),
              sd(results_df$ate_so9_scm, na.rm = TRUE)))
  cat(sprintf("  Significant results (p < 0.05): %d / %d (%.1f%%)\n",
              sum(results_df$ate_percentile_so9 > 0.95, na.rm = TRUE),
              sum(!is.na(results_df$ate_percentile_so9)),
              100 * sum(results_df$ate_percentile_so9 > 0.95, na.rm = TRUE) / sum(!is.na(results_df$ate_percentile_so9))))

  cat("\nFIP Results:\n")
  cat(sprintf("  ATE SCM - Mean: %.2f, Median: %.2f, SD: %.2f\n",
              mean(results_df$ate_fip_scm, na.rm = TRUE),
              median(results_df$ate_fip_scm, na.rm = TRUE),
              sd(results_df$ate_fip_scm, na.rm = TRUE)))
  cat(sprintf("  Significant results (p < 0.05): %d / %d (%.1f%%)\n",
              sum(results_df$ate_percentile_fip > 0.95, na.rm = TRUE),
              sum(!is.na(results_df$ate_percentile_fip)),
              100 * sum(results_df$ate_percentile_fip > 0.95, na.rm = TRUE) / sum(!is.na(results_df$ate_percentile_fip))))
}

cat("\nDone!\n")

mean(results_df$so9_year1)
mean(results_df$so9_year2, na.rm = T)
mean(results_df$so9_year3, na.rm = T)

mean(results_df$fip_year1)
mean(results_df$fip_year2, na.rm = T)
mean(results_df$fip_year3, na.rm = T)

#--- 7. Statistical Significance Tests for Year-by-Year Effects ----------------------------------------
cat("\n========================================\n")
cat("YEAR-BY-YEAR COMPARISON TESTS\n")
cat("========================================\n\n")

cat("SO9 Year-by-Year Comparisons:\n")
cat("----------------------------------------\n")

#SO9: Year 1 vs Year 2
so9_year1_year2 <- results_df %>%
  filter(!is.na(so9_year1) & !is.na(so9_year2))

if (nrow(so9_year1_year2) > 1) {
  so9_test_1_2 <- t.test(so9_year1_year2$so9_year1, so9_year1_year2$so9_year2, paired = TRUE)
  cat(sprintf("SO9 Year 1 vs Year 2 (n=%d):\n", nrow(so9_year1_year2)))
  cat(sprintf("  Year 1 mean: %.3f\n", mean(so9_year1_year2$so9_year1)))
  cat(sprintf("  Year 2 mean: %.3f\n", mean(so9_year1_year2$so9_year2)))
  cat(sprintf("  Mean difference: %.3f\n", mean(so9_year1_year2$so9_year1 - so9_year1_year2$so9_year2)))
  cat(sprintf("  t-statistic: %.3f\n", so9_test_1_2$statistic))
  cat(sprintf("  p-value: %.4f %s\n", so9_test_1_2$p.value,
              ifelse(so9_test_1_2$p.value < 0.05, "***", "")))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n\n", so9_test_1_2$conf.int[1], so9_test_1_2$conf.int[2]))
} else {
  cat("SO9 Year 1 vs Year 2: Insufficient data\n\n")
}

#SO9: Year 2 vs Year 3
so9_year2_year3 <- results_df %>%
  filter(!is.na(so9_year2) & !is.na(so9_year3))

if (nrow(so9_year2_year3) > 1) {
  so9_test_2_3 <- t.test(so9_year2_year3$so9_year2, so9_year2_year3$so9_year3, paired = TRUE)
  cat(sprintf("SO9 Year 2 vs Year 3 (n=%d):\n", nrow(so9_year2_year3)))
  cat(sprintf("  Year 2 mean: %.3f\n", mean(so9_year2_year3$so9_year2)))
  cat(sprintf("  Year 3 mean: %.3f\n", mean(so9_year2_year3$so9_year3)))
  cat(sprintf("  Mean difference: %.3f\n", mean(so9_year2_year3$so9_year2 - so9_year2_year3$so9_year3)))
  cat(sprintf("  t-statistic: %.3f\n", so9_test_2_3$statistic))
  cat(sprintf("  p-value: %.4f %s\n", so9_test_2_3$p.value,
              ifelse(so9_test_2_3$p.value < 0.05, "***", "")))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n\n", so9_test_2_3$conf.int[1], so9_test_2_3$conf.int[2]))
} else {
  cat("SO9 Year 2 vs Year 3: Insufficient data\n\n")
}

cat("FIP Year-by-Year Comparisons:\n")
cat("----------------------------------------\n")

#FIP: Year 1 vs Year 2
fip_year1_year2 <- results_df %>%
  filter(!is.na(fip_year1) & !is.na(fip_year2))

if (nrow(fip_year1_year2) > 1) {
  fip_test_1_2 <- t.test(fip_year1_year2$fip_year1, fip_year1_year2$fip_year2, paired = TRUE)
  cat(sprintf("FIP Year 1 vs Year 2 (n=%d):\n", nrow(fip_year1_year2)))
  cat(sprintf("  Year 1 mean: %.3f\n", mean(fip_year1_year2$fip_year1)))
  cat(sprintf("  Year 2 mean: %.3f\n", mean(fip_year1_year2$fip_year2)))
  cat(sprintf("  Mean difference: %.3f\n", mean(fip_year1_year2$fip_year1 - fip_year1_year2$fip_year2)))
  cat(sprintf("  t-statistic: %.3f\n", fip_test_1_2$statistic))
  cat(sprintf("  p-value: %.4f %s\n", fip_test_1_2$p.value,
              ifelse(fip_test_1_2$p.value < 0.05, "***", "")))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n\n", fip_test_1_2$conf.int[1], fip_test_1_2$conf.int[2]))
} else {
  cat("FIP Year 1 vs Year 2: Insufficient data\n\n")
}

#FIP: Year 2 vs Year 3
fip_year2_year3 <- results_df %>%
  filter(!is.na(fip_year2) & !is.na(fip_year3))

if (nrow(fip_year2_year3) > 1) {
  fip_test_2_3 <- t.test(fip_year2_year3$fip_year2, fip_year2_year3$fip_year3, paired = TRUE)
  cat(sprintf("FIP Year 2 vs Year 3 (n=%d):\n", nrow(fip_year2_year3)))
  cat(sprintf("  Year 2 mean: %.3f\n", mean(fip_year2_year3$fip_year2)))
  cat(sprintf("  Year 3 mean: %.3f\n", mean(fip_year2_year3$fip_year3)))
  cat(sprintf("  Mean difference: %.3f\n", mean(fip_year2_year3$fip_year2 - fip_year2_year3$fip_year3)))
  cat(sprintf("  t-statistic: %.3f\n", fip_test_2_3$statistic))
  cat(sprintf("  p-value: %.4f %s\n", fip_test_2_3$p.value,
              ifelse(fip_test_2_3$p.value < 0.05, "***", "")))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n\n", fip_test_2_3$conf.int[1], fip_test_2_3$conf.int[2]))
} else {
  cat("FIP Year 2 vs Year 3: Insufficient data\n\n")
}

cat("Note: *** indicates p < 0.05 (statistically significant)\n")
cat("Positive differences indicate effect magnitude decreases over time.\n")

#--- 8. Create Year-by-Year Effect Graphs ----------------------------------------
cat("\n========================================\n")
cat("CREATING YEAR-BY-YEAR EFFECT GRAPHS\n")
cat("========================================\n\n")

library(ggplot2)

#SO9 Year-by-Year Graph
so9_yearly_data <- data.frame(
  Year = c(1, 2, 3),
  Mean_Effect = c(
    mean(results_df$so9_year1, na.rm = TRUE),
    mean(results_df$so9_year2, na.rm = TRUE),
    mean(results_df$so9_year3, na.rm = TRUE)
  ),
  N = c(
    sum(!is.na(results_df$so9_year1)),
    sum(!is.na(results_df$so9_year2)),
    sum(!is.na(results_df$so9_year3))
  ),
  SE = c(
    sd(results_df$so9_year1, na.rm = TRUE) / sqrt(sum(!is.na(results_df$so9_year1))),
    sd(results_df$so9_year2, na.rm = TRUE) / sqrt(sum(!is.na(results_df$so9_year2))),
    sd(results_df$so9_year3, na.rm = TRUE) / sqrt(sum(!is.na(results_df$so9_year3)))
  )
)

plot_so9_yearly <- ggplot(so9_yearly_data, aes(x = Year, y = Mean_Effect)) +
  geom_line(color = "#2E86AB", size = 1.5) +
  geom_point(color = "#2E86AB", size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
  labs(
    title = "SO9 Treatment Effect Over Time Post-Surgery",
    subtitle = "Mean year-by-year effects with 95% confidence intervals",
    x = "Years Post-Surgery",
    y = "Treatment Effect (SO9)"
  ) +
  scale_x_continuous(breaks = 1:3, labels = c("Year 1", "Year 2", "Year 3")) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12, face = "bold")
  )

ggsave("so9_yearly_effects.png", plot = plot_so9_yearly, width = 10, height = 6, dpi = 300)
cat("SO9 yearly effects graph saved to: so9_yearly_effects.png\n")

#FIP Year-by-Year Graph
fip_yearly_data <- data.frame(
  Year = c(1, 2, 3),
  Mean_Effect = c(
    mean(results_df$fip_year1, na.rm = TRUE),
    mean(results_df$fip_year2, na.rm = TRUE),
    mean(results_df$fip_year3, na.rm = TRUE)
  ),
  N = c(
    sum(!is.na(results_df$fip_year1)),
    sum(!is.na(results_df$fip_year2)),
    sum(!is.na(results_df$fip_year3))
  ),
  SE = c(
    sd(results_df$fip_year1, na.rm = TRUE) / sqrt(sum(!is.na(results_df$fip_year1))),
    sd(results_df$fip_year2, na.rm = TRUE) / sqrt(sum(!is.na(results_df$fip_year2))),
    sd(results_df$fip_year3, na.rm = TRUE) / sqrt(sum(!is.na(results_df$fip_year3)))
  )
)

plot_fip_yearly <- ggplot(fip_yearly_data, aes(x = Year, y = Mean_Effect)) +
  geom_line(color = "#C73E1D", size = 1.5) +
  geom_point(color = "#C73E1D", size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
  labs(
    title = "FIP Treatment Effect Over Time Post-Surgery",
    subtitle = "Mean year-by-year effects with 95% confidence intervals",
    x = "Years Post-Surgery",
    y = "Treatment Effect (FIP)"
  ) +
  scale_x_continuous(breaks = 1:3, labels = c("Year 1", "Year 2", "Year 3")) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12, face = "bold")
  )

ggsave("fip_yearly_effects.png", plot = plot_fip_yearly, width = 10, height = 6, dpi = 300)
cat("FIP yearly effects graph saved to: fip_yearly_effects.png\n")

cat("\nAll graphs created successfully!\n")

