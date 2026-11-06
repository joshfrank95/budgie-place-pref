# This script extracts spectral features from the stimuli used in the Kathios study.
# The frequency-amplitude spectra extracted in 1-visualise-kathios-stimuli.r are required.
# Please run source("1-visualise-kathios-stimuli.r") first if you have not already run that script.

if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("pmcharrison/hrep")
devtools::install_github("pmcharrison/incon")

library(hrep)
library(incon)
library(vegan)

# We now calculate some spectral features,
# first normalising the spectrum to the maximum amplitude.

stimuli_labelled <- stimuli_labelled %>%
    dplyr::mutate(labelled_spectrum = map(
        labelled_spectrum, normalise_spectrum))

# Spectral centroid: from trimmed FFT,
# for comparison with Kathios et al.

get_spectral_centroid <- function(spectrum) {
    spectrum %>%
        dplyr::filter(freq >= 25) %>%
        dplyr::mutate(weighted = amp * freq) %>%
        dplyr::summarize(spectral_centroid = sum(weighted) / sum(amp)) %>%
        dplyr::pull(spectral_centroid)
}

spectral_centroids <- stimuli_labelled %>%
    dplyr::mutate(spectral_centroid = map_dbl(labelled_spectrum, get_spectral_centroid)) %>% 
    dplyr::select(stimulus, spectral_centroid)

View(spectral_centroids)

# Roughness: Hutchinson & Knopoff (1978) method.
# This is the best-performing method for correspondence
# with human aesthetic judgments in prior literature
# (see Harrison & Pearce, 2020; Eerola & Lahdelma, 2021).

get_rough <- function(spectrum){
    spectrum <- spectrum %>% dplyr::filter(freq >= 25)
    spect <- hrep::sparse_fr_spectrum(list(frequency = spectrum$freq, amplitude = spectrum$amp))
    rough <- dycon::roughness_hutch(spect)
    return(rough)
}

rough_vals <- stimuli_labelled %>%
    dplyr::mutate(rough = map_dbl(labelled_spectrum, get_rough)) %>% 
    dplyr::select(stimulus, rough)

View(rough_vals)

# Harmonicity: Harrison & Pearce (2018) method.
# This is the best-performing method that can account
# for the actual spectrum of the stimulus (Harrison & Pearce, 2020);
# an alternative (Stolzenburg, 2015; Eerola & Lahdelma, 2021) imposes
# a harmonic spectrum above the fundamental/s.

#Calculate a chord's harmonicity using either the model of Harrison & Pearce (2018).
# The input x must be a sparse_pi_spectrum.
.get_harmonicity <- function(x, mode) {
  pc_spectrum <- hrep::smooth_pc_spectrum(x)
  pitch_profile <- get_pitch_profile(pc_spectrum)
  eval_pitch_profile(pitch_profile, mode)
}

get_harmonicity <- memoise::memoise(.get_harmonicity)

get_pitch_profile <- function(x) {
  stopifnot(is(x, "smooth_pc_spectrum"))
  x <- as.numeric(x)
  array_dim <- length(x)
  template <- hrep::smooth_pc_spectrum(
    hrep::pi_chord(60),
    array_dim = array_dim
  )
  res <- har18::sweep_template(x, template)
  hrep::.smooth_pc_spectrum(res)
}

eval_pitch_profile <- function(profile, mode) {
  if (mode == "harrison") {
    kl_div_from_uniform(as.numeric(profile))
  } else {
    stop("unrecognised mode: ", mode)
  }
}

kl_div_from_uniform <- function(x) {
  probs <- x / sum(x)
  n <- length(probs)
  uniform_probs <- 1 / n
  non_zero_probs <- probs[probs > 0]
  sum(
    non_zero_probs * log(non_zero_probs / uniform_probs, base = 2)
  )
}

hfun <- function(spectrum){
    spectrum <- spectrum %>% dplyr::filter(freq >= 25)
    spect <- hrep::sparse_pi_spectrum(list(pitch = freq_to_midi(spectrum$freq), amplitude = spectrum$amp))
    h_val <- get_harmonicity(spect, "harrison")
    return(h_val)
}

freq_to_midi <- function(freq) {
    69 + 12 * log2(freq / 440)
}

harmonicity_vals <- stimuli_labelled %>%
    dplyr::mutate(harmonicity = map_dbl(labelled_spectrum, hfun)) %>% 
    dplyr::select(stimulus, harmonicity)

View(harmonicity_vals)

# MPS: extracted using mps-extractor in this repository
# Settings: all default, except:
# - return_percentage: this expresses the power in the specified AM range as a percentage of total spectrogram power
# - rms_weighted: this weights the AM per time-slice by the corresponding RMS value
# - tmr_low 20 OR 30 (depending on which MPS you are extracting)
# - tmr_high 200 OR 150 (depending on which MPS you are extracting)


# MPS 20-200Hz:
mps_20_200 <- read_csv('outputs/mps-outputs/mps_20_200.csv')

mps_20_200 <- mps_20_200 %>%
    dplyr::group_by(file) %>%
    dplyr::summarise(mean_mps_20_200 = mean(mps_roughness)) %>% 
    dplyr::mutate(stimulus = basename(file)) %>% 
    dplyr::mutate(stimulus = str_remove(stimulus, "_trimmed.wav")) %>% 
    dplyr::mutate(stimulus = str_remove(stimulus, "stim")) %>% 
    dplyr::select(stimulus, mean_mps_20_200)

View(mps_20_200)

# MPS 30-150Hz:
mps_30_150 <- read_csv('outputs/mps-outputs/mps_30_150.csv')

mps_30_150 <- mps_30_150 %>%
    dplyr::group_by(file) %>%
    dplyr::summarise(mean_mps_30_150 = mean(mps_roughness)) %>% 
    dplyr::mutate(stimulus = basename(file)) %>% 
    dplyr::mutate(stimulus = str_remove(stimulus, "_trimmed.wav")) %>% 
    dplyr::mutate(stimulus = str_remove(stimulus, "stim")) %>% 
    dplyr::select(stimulus, mean_mps_30_150)

View(mps_30_150)

# Join feature values with prior data.

stimuli_labelled <- stimuli_labelled %>% 
    dplyr::mutate(
        condition = case_when(
            stimulus == "C" | stimulus == "E" | stimulus == "I" | stimulus == "F" | stimulus == "J" | stimulus == "K" ~ "pleasant",
            stimulus == "A" | stimulus == "G" | stimulus == "L" | stimulus == "H" | stimulus == "B" | stimulus == "D" ~ "unpleasant"
        )
    )

stimuli_comp <- stimuli_labelled %>%
    dplyr::left_join(spectral_centroids, by = "stimulus") %>%
    dplyr::left_join(rough_vals, by = "stimulus") %>%
    dplyr::left_join(harmonicity_vals, by = "stimulus") %>% 
    dplyr::left_join(mps_20_200, by = "stimulus") %>%
    dplyr::left_join(mps_30_150, by = "stimulus")

# Optional: save as RDS to avoid later re-computation.
# saveRDS(stimuli_comp, "stimuli_comp.rds")

View(stimuli_comp)

means_per_condition <- stimuli_comp %>%
    dplyr::group_by(condition) %>%
    dplyr::summarize(
        mean_spectral_centroid = mean(spectral_centroid),
        mean_rough = mean(rough),
        mean_harmonicity = mean(harmonicity),
        grand_mean_mps_20_200 = mean(mean_mps_20_200),
        grand_mean_mps_30_150 = mean(mean_mps_30_150),
        sd_spectral_centroid = sd(spectral_centroid),
        sd_rough = sd(rough),
        sd_harmonicity = sd(harmonicity),
        sd_mps_20_200 = sd(mean_mps_20_200),
        sd_mps_30_150 = sd(mean_mps_30_150)
    )

View(means_per_condition)

# Do conditions differ?

numerics <- stimuli_comp %>%
    dplyr::select(-condition, -stimulus, -spectrum, -labelled_spectrum, -F01, -F02, -file) %>% 
    dplyr::mutate(across(everything(), scale))

nums_ranked <- apply(numerics, 2, rank)

# Same result with either ranked or unranked data.

set.seed(123)
p_res <- adonis2(
    numerics ~ condition,
    data = stimuli_comp,
    method = "euclidean",
    permutations = 999
)

set.seed(123)
disp <- betadisper(
    dist(numerics), stimuli_comp$condition,
)

set.seed(123)
permutest(disp, permutations = 999)

set.seed(123)
p_res_ranked <- adonis2(
    nums_ranked ~ condition,
    data = stimuli_comp,
    method = "euclidean",
    permutations = 999
)

set.seed(123)
disp_ranked <- betadisper(
    dist(nums_ranked), stimuli_comp$condition,
)

set.seed(123)
permutest(disp_ranked, permutations = 999)

u_test_res <- tibble(
    variable = character(),
    test_statistic = numeric(),
    p_value_raw = numeric(),
    p_value_adjusted = numeric()
)
for (v in names(numerics)) {
    result <- wilcox.test(stimuli_comp[[v]] ~ stimuli_comp$condition)
    u_test_res <- u_test_res %>%
        add_row(
            variable = v,
            test_statistic = result$statistic,
            p_value_raw = result$p.value,
            p_value_adjusted = p.adjust(result$p.value, method = "bonferroni", n = 5)
        )
}
