# This script tests the human ratings against the spectral features.
# Human ratings are stored in the "ratings" folder.
# The script requires "stimuli_comp" from 2-get-spectral-features-kathios.r.
# Please run source ("1-visualise-kathios-stimuli.r") first if you have not already run that script.
# Please run source ("2-get-spectral-features-kathios.r") second if you have not already run that script.

library(parallel)
library(boot)
library(broom)
library(MASS)
library(repmod)
library(performance)
library(car)

# Read in the ratings, extract and label relevant outcomes
ratings <- read_csv("ratings/human-ratings.csv")

ratings <- ratings[,29:ncol(ratings)]

pleasantness <- ratings[, seq(1, ncol(ratings), 3)]
colnames(pleasantness) <- c("C", "E", "F", "I", "J", "K", "A", "B", "D", "G", "H", "L")

liking <- ratings[, seq(2, ncol(ratings), 3)]
colnames(liking) <- c("C", "E", "F", "I", "J", "K", "A", "B", "D", "G", "H", "L")

interest <- ratings[, seq(3, ncol(ratings), 3)]
colnames(interest) <- c("C", "E", "F", "I", "J", "K", "A", "B", "D", "G", "H", "L")

# Normalise within participants, then take the mean across participants
pleasantness <- pleasantness %>%
    rowwise() %>%
    mutate(across(everything(), ~(. - mean(c_across(everything()), na.rm=TRUE)) / sd(c_across(everything()), na.rm=TRUE))) %>% 
    ungroup()

pleasantness_avgs <- pleasantness %>%
    dplyr::summarise(across(everything(), \(x) mean(x, na.rm=TRUE))) %>% 
    pivot_longer(everything(), names_to = "stimulus", values_to = "pleasantness")

liking <- liking %>%
    rowwise() %>%
    mutate(across(everything(), ~(. - mean(c_across(everything()), na.rm=TRUE)) / sd(c_across(everything()), na.rm=TRUE))) %>% 
    ungroup()

interest <- interest %>%
    rowwise() %>%
    mutate(across(everything(), ~(. - mean(c_across(everything()), na.rm=TRUE)) / sd(c_across(everything()), na.rm=TRUE))) %>% 
    ungroup()

liking_avgs <- liking %>%
    dplyr::summarise(across(everything(), \(x) mean(x, na.rm=TRUE))) %>% 
    pivot_longer(everything(), names_to = "stimulus", values_to = "liking")

interest_avgs <- interest %>%
    dplyr::summarise(across(everything(), \(x) mean(x, na.rm=TRUE))) %>% 
    pivot_longer(everything(), names_to = "stimulus", values_to = "interest")

# Join the ratings with the spectral features
dflm <- left_join(stimuli_comp, pleasantness_avgs, by = c("stimulus" = "stimulus"))
dflm <- left_join(dflm, liking_avgs, by = c("stimulus" = "stimulus"))
dflm <- left_join(dflm, interest_avgs, by = c("stimulus" = "stimulus"))

dflm <- dflm %>% dplyr::select(-F01, -F02, -spectrum, -labelled_spectrum, -file)


# Check correlations between variables

dflm <- dflm %>% dplyr::mutate(across(c(spectral_centroid, rough, harmonicity, mean_mps_20_200, mean_mps_30_150, pleasantness, liking, interest), \(x) scale(x)))

cors <- Hmisc::rcorr(dflm %>% dplyr::select(-stimulus, -condition) %>% as.matrix(), type = "spearman")

cor_matrix <- cors$r  # Extract the correlation matrix
p_matrix <- cors$P    # Extract the p-values

# Convert to long format
cor_data <- reshape2::melt(cor_matrix)
p_data <- reshape2::melt(p_matrix)

# Merge correlation and p-value data for plotting
plot_data <- merge(cor_data, p_data, by = c("Var1", "Var2"))
colnames(plot_data) <- c("Var1", "Var2", "Correlation", "P_value")

var_names <- c(
    "spectral_centroid" = "Spectral Centroid",
    "rough" = "Roughness (H&K)",
    "harmonicity" = "Harmonicity (H&P)",
    "mean_mps_20_200" = "20-200Hz AM",
    "mean_mps_30_150" = "30-150Hz AM",
    "pleasantness" = "Pleasantness",
    "liking" = "Liking",
    "interest" = "Interest"
)

# Create the heatmap with ggplot
corplot <- ggplot(data = plot_data, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  scale_x_discrete(labels = var_names) +
  scale_y_discrete(labels = var_names) +
  coord_fixed()

ggsave("outputs/human-ratings-corplot.png", corplot, width = 12, height = 12, dpi = 600)

# Regression models

# BCA function ----
bca_model <- function(formula, data, mode, ..., .nboot = 1000) {

  set.seed(1)

  statistic <- function(data, indices) {
    sampled_data <- data[indices, ]
    warning_flag <- FALSE
    log_file <- NULL
    mean_test_retest_var <- NULL
    if (mode == "ols"){
      mod <- lm(formula = formula, data = sampled_data, ...)

      fit_stats <-
        broom::glance(mod) |>
        dplyr::select(
          "r_squared" = r.squared,
          "adj_r_squared" = adj.r.squared
        ) |>
        unlist()
    } else if (mode == "robust") {
      mod <- MASS::rlm(
        formula = formula,
        data = sampled_data,
        seed = 1,
        method = "MM",
        init = "ls",
        acc = 1e-15,
        maxit = 10000,
        na.action = "na.fail"
      )
      fit_stats <- NULL
    } else if (mode == "sem") {
      mod <- withCallingHandlers(
        expr = lavaan(
          model = formula, data = sampled_data, model.type = "sem", ...),
        warning = function(w) {
          warning_flag <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
      fit_stats <- NULL
    }

    coef_stats <- coef(mod)

    names(coef_stats) <- gsub(
      "(Intercept)", "intercept", names(coef_stats), fixed = TRUE)
    names(coef_stats) <- paste("coef_", names(coef_stats), sep = "")

    if (mode == "sem") {
      test_coefs <- as.numeric(coef_stats[grep("coef_test", names(coef_stats))])
      retest_coefs <- as.numeric(coef_stats[grep("coef_retest", names(coef_stats))])
      if (length(test_coefs) > 0 && length(retest_coefs) > 0) {
        mean_test_retest_var <- mean(c(test_coefs, retest_coefs))
      } else {
        mean_test_retest_var <- 0
      }
      mean_test_retest_var <- rlang::set_names(
        mean_test_retest_var, "mean_test_retest_var")
    }


    coef_stats <- c(fit_stats, coef_stats, mean_test_retest_var)

    if (warning_flag) {
      dir.create("debug_logs", showWarnings = FALSE)
      log_file <- paste0("debug_logs/debug_log_", Sys.time(), ".rds")
      saveRDS(
        list(sampled_data = sampled_data, formula = formula, model = mod),
        log_file
      )
    }

    coef_stats

  }

  if (mode == "ols") {
    non_bootstrapped_model <- lm(formula = formula, data = data, ...)
  } else if (mode == "robust") {
    non_bootstrapped_model <- withCallingHandlers(
      expr = {
        MASS::rlm(
          formula = formula, 
          data = data, 
          seed = 1,
          method = "MM",
          init = "ls",
          acc = 1e-15,
          maxit = 10000,
          na.action = "na.fail"
        )
      },
      warning = function(w) {
        warning_flag <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
  } else if (mode == "sem") {
    non_bootstrapped_model <- lavaan(
      model = formula, data = data, model.type = "sem")
  }

  non_bootstrapped_stats <- statistic(data, indices = seq_len(nrow(data)))

  res_boot <- boot::boot(
    data = data,
    statistic = statistic,
    R = .nboot,
    parallel = "multicore",
    ncpus = detectCores() - 1
  )

  df <- tibble::tibble(
    parameter = names(non_bootstrapped_stats),
    estimate = non_bootstrapped_stats,
    ci_95_lower = NA_real_,
    ci_95_upper = NA_real_
  )

  for (i in seq_len(nrow(df))) {
    res_boot_bca <- boot::boot.ci(res_boot, type = "bca", index = i)
    ci <- res_boot_bca$bca[4:5]
    df$ci_95_lower[i] <- ci[1]
    df$ci_95_upper[i] <- ci[2]
  }

  warning_count <- sum(apply(res_boot$t, 1, function(row) any(is.na(row))))

  list(
    par = df,
    mod = non_bootstrapped_model,
    nboot = .nboot,
    warnings = warning_count
  )
}

# Pleasantness models ----
# Roughness metric: Hutchinson & Knopoff
fit_pleasantness_hk <- lm(
    pleasantness ~ rough + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_pleasantness_hk)
performance::check_model(fit_pleasantness_hk)
car::vif(fit_pleasantness_hk)

fit_pleasantness_hk_ci <- bca_model(
    pleasantness ~ rough + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_pleasantness_hk_ci$par, n = 25)

# Roughness metric: MPS (20-200Hz) ----

fit_pleasantness_mps_20_200 <- lm(
    pleasantness ~ mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_pleasantness_mps_20_200)
performance::check_model(fit_pleasantness_mps_20_200)
performance::check_outliers(fit_pleasantness_mps_20_200)
car::vif(fit_pleasantness_mps_20_200)

fit_pleasantness_mps_20_200_ci <- bca_model(
    pleasantness ~ mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_pleasantness_mps_20_200_ci$par, n = 25)

# Roughness metric: MPS (30-150Hz) ----
fit_pleasantness_mps_30_150 <- lm(
    pleasantness ~ mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_pleasantness_mps_30_150)
performance::check_model(fit_pleasantness_mps_30_150)
performance::check_outliers(fit_pleasantness_mps_30_150)
car::vif(fit_pleasantness_mps_30_150)

fit_pleasantness_mps_30_150_ci <- bca_model(
    pleasantness ~ mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_pleasantness_mps_30_150_ci$par, n = 25)

# Roughness metric: H&K and MPS (20-200Hz) ----
fit_pleasantness_hk_mps_20_200 <- lm(
    pleasantness ~ rough + mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_pleasantness_hk_mps_20_200)
performance::check_model(fit_pleasantness_hk_mps_20_200)
performance::check_outliers(fit_pleasantness_hk_mps_20_200)
car::vif(fit_pleasantness_hk_mps_20_200)

fit_pleasantness_hk_mps_20_200_ci <- bca_model(
    pleasantness ~ rough + mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_pleasantness_hk_mps_20_200_ci$par, n = 25)

# Roughness metric: H&K and MPS (30-150Hz) ----
fit_pleasantness_hk_mps_30_150 <- lm(
    pleasantness ~ rough + mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_pleasantness_hk_mps_30_150)
performance::check_model(fit_pleasantness_hk_mps_30_150)
performance::check_outliers(fit_pleasantness_hk_mps_30_150)
car::vif(fit_pleasantness_hk_mps_30_150)

fit_pleasantness_hk_mps_30_150_ci <- bca_model(
    pleasantness ~ rough + mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_pleasantness_hk_mps_30_150_ci$par, n = 25)

# Liking models ----
# Roughness metric: Hutchinson & Knopoff
fit_liking_hk <- lm(
    liking ~ rough + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_liking_hk)
performance::check_model(fit_liking_hk)
car::vif(fit_liking_hk)

fit_liking_hk_ci <- bca_model(
    liking ~ rough + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_liking_hk_ci$par, n = 25)

# Roughness metric: MPS (20-200Hz) ----
fit_liking_mps_20_200 <- lm(
    liking ~ mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_liking_mps_20_200)
performance::check_model(fit_liking_mps_20_200)
performance::check_outliers(fit_liking_mps_20_200)
car::vif(fit_liking_mps_20_200)

fit_liking_mps_20_200_ci <- bca_model(
    liking ~ mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_liking_mps_20_200_ci$par, n = 25)

# Roughness metric: MPS (30-150Hz) ----
fit_liking_mps_30_150 <- lm(
    liking ~ mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_liking_mps_30_150)
performance::check_model(fit_liking_mps_30_150)
performance::check_outliers(fit_liking_mps_30_150)
car::vif(fit_liking_mps_30_150)

fit_liking_mps_30_150_ci <- bca_model(
    liking ~ mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_liking_mps_30_150_ci$par, n = 25)

# Roughness metric: H&K and MPS (20-200Hz) ----
fit_liking_hk_mps_20_200 <- lm(
    liking ~ rough + mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_liking_hk_mps_20_200)
performance::check_model(fit_liking_hk_mps_20_200)
performance::check_outliers(fit_liking_hk_mps_20_200)
car::vif(fit_liking_hk_mps_20_200)

fit_liking_hk_mps_20_200_ci <- bca_model(
    liking ~ rough + mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_liking_hk_mps_20_200_ci$par, n = 25)

# Roughness metric: H&K and MPS (30-150Hz) ----
fit_liking_hk_mps_30_150 <- lm(
    liking ~ rough + mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_liking_hk_mps_30_150)
performance::check_model(fit_liking_hk_mps_30_150)
performance::check_outliers(fit_liking_hk_mps_30_150)
car::vif(fit_liking_hk_mps_30_150)

fit_liking_hk_mps_30_150_ci <- bca_model(
    liking ~ rough + mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_liking_hk_mps_30_150_ci$par, n = 25)

# Interest models ----
# Roughness metric: Hutchinson & Knopoff
fit_interest_hk <- lm(
    interest ~ rough + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_interest_hk)
performance::check_model(fit_interest_hk)
car::vif(fit_interest_hk)

fit_interest_hk_ci <- bca_model(
    interest ~ rough + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_interest_hk_ci$par, n = 25)

# Roughness metric: MPS (20-200Hz) ----
fit_interest_mps_20_200 <- lm(
    interest ~ mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_interest_mps_20_200)
performance::check_model(fit_interest_mps_20_200)
performance::check_outliers(fit_interest_mps_20_200)
car::vif(fit_interest_mps_20_200)

fit_interest_mps_20_200_ci <- bca_model(
    interest ~ mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_interest_mps_20_200_ci$par, n = 25)

# Roughness metric: MPS (30-150Hz) ----
fit_interest_mps_30_150 <- lm(
    interest ~ mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_interest_mps_30_150)
performance::check_model(fit_interest_mps_30_150)
performance::check_outliers(fit_interest_mps_30_150)
car::vif(fit_interest_mps_30_150)

fit_interest_mps_30_150_ci <- bca_model(
    interest ~ mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_interest_mps_30_150_ci$par, n = 25)

# Roughness metric: H&K and MPS (20-200Hz) ----
fit_interest_hk_mps_20_200 <- lm(
    interest ~ rough + mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_interest_hk_mps_20_200)
performance::check_model(fit_interest_hk_mps_20_200)
performance::check_outliers(fit_interest_hk_mps_20_200)
car::vif(fit_interest_hk_mps_20_200)

fit_interest_hk_mps_20_200_ci <- bca_model(
    interest ~ rough + mean_mps_20_200 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_interest_hk_mps_20_200_ci$par, n = 25)

# Roughness metric: H&K and MPS (30-150Hz) ----
fit_interest_hk_mps_30_150 <- lm(
    interest ~ rough + mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm
)
summary(fit_interest_hk_mps_30_150)
performance::check_model(fit_interest_hk_mps_30_150)
performance::check_outliers(fit_interest_hk_mps_30_150)
car::vif(fit_interest_hk_mps_30_150)

fit_interest_hk_mps_30_150_ci <- bca_model(
    interest ~ rough + mean_mps_30_150 + harmonicity + spectral_centroid,
    data = dflm,
    mode = "ols",
    na.action = "na.fail",
    .nboot = 10000
)
print(fit_interest_hk_mps_30_150_ci$par, n = 25)

