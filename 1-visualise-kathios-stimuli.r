# This script is used to visualise the stimuli used in the Kathios study.
# The stimuli are stored in the "stimuli-kathios" folder.

library(tidyverse)
library(tuneR)
library(seewave)
library(cowplot)

#### LOAD DATA ####

# The stimuli have been trimmed to remove silence

stimulus_dir <- "stimuli-kathios"

files <- list.files(stimulus_dir, pattern = "\\.wav$", full.names = TRUE)

# For each file, we read the .wav,
# perform an FFT, and extract peaks above 1% amplitude.
# We set window length to the length of the file for maximum resolution.
# The peak-picking process is quite slow, so expect a minute or two of waiting.

stimuli <- tibble(
    file = files,
    stimulus = str_remove(
        file,
        "stimuli-kathios/stim")) %>% 
        dplyr::mutate(
            stimulus = str_remove(
                stimulus,
                "_trimmed.wav")
        ) %>% 
    dplyr::mutate(
        spectrum = purrr::map(file, function(file) {
            w <- tuneR::readWave(file)
            res <- seewave::spec(w, f = 44100, wl = length(w), plot = FALSE)
            peaks <- seewave::fpeaks(res, threshold = 0.01)
            return(tibble::as_tibble(peaks))
        })
    )

# Now we add the fundamental frequencies for each stimulus,
# taken from Kathios et al.

stimuli <- stimuli %>%
    dplyr::mutate(
        stimulus = str_remove(stimulus, ".wav"),
        F01 = case_when(
            stimulus == "C" ~ 100,
            stimulus == "E" ~ 125,
            stimulus == "I" ~ 110, 
            stimulus == "F" ~ 110,
            stimulus == "J" ~ 125,
            stimulus == "K" ~ 440,
            stimulus == "A" ~ 110,
            stimulus == "G" ~ 1000,
            stimulus == "L" ~ 1745,
            stimulus == "H" ~ 1745,
            stimulus == "B" ~ 125,
            stimulus == "D" ~ 950
    ),
    F02 = case_when(
        stimulus == "C" ~ 124.5,
        stimulus == "E" ~ 156,
        stimulus == "I" ~ 330,
        stimulus == "F" ~ 880,
        stimulus == "J" ~ 126,
        stimulus == "K" ~ 880,
        stimulus == "A" ~ 1050,
        stimulus == "G" ~ 1050,
        stimulus == "L" ~ 1799,
        stimulus == "H" ~ 1799,
        stimulus == "B" ~ 130,
        stimulus == "D" ~ 1000
    )
)

# We now label the spectrum with its harmonics.
# This is done by checking if the frequency is
# a (rough) multiple of one of the fundamental frequencies.
# We also convert the frequencies from kHz to Hz,
# and round the frequencies to the nearest integer
# to avoid minor errors.

spectA <- stimuli$spectrum[[1]]
spectA <- spectA %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 110 == 0 ~ "110",
        round(freq, 0) %% 1050 == 0 ~ "1050",
        TRUE ~ NA_character_
    )
)

spectB <- stimuli$spectrum[[2]]
spectB <- spectB %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 125 == 0 ~ 125,
        round(freq, 0) %% 130 == 0 ~ 130,
        TRUE ~ NA_real_
    )
)

#Some generosity needed here...
spectC <- stimuli$spectrum[[3]]
spectC <- spectC %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 100 == 0 ~ 100,
        round(freq, 0) %% 124.5 > 2 | round(freq, 0) %% 124.5 < 122.5 ~ 124.5,
        TRUE ~ NA_real_
    )
)

spectD <- stimuli$spectrum[[4]]
spectD <- spectD %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 950 == 0 ~ 950,
        round(freq, 0) %% 1000 == 0 ~ 1000,
        TRUE ~ NA_real_
    )
)

spectE <- stimuli$spectrum[[5]]
spectE <- spectE %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 125 == 0 ~ 125,
        round(freq, 0) %% 156 == 0 ~ 156,
        TRUE ~ NA_real_
    )
)

spectF <- stimuli$spectrum[[6]]
spectF <- spectF %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 880 < 1 ~ 880,
        round(freq, 0) %% 110 < 1 ~ 110,
        TRUE ~ NA_real_
    )
)

spectG <- stimuli$spectrum[[7]]
spectG <- spectG %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 1000 == 0 ~ 1000,
        round(freq, 0) %% 1050 == 0 ~ 1050,
        TRUE ~ NA_real_
    )
)

spectH <- stimuli$spectrum[[8]]
spectH <- spectH %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 1745 == 0 ~ 1745,
        round(freq, 0) %% 1799 == 0 ~ 1799,
        TRUE ~ NA_real_
    )
)

spectI <- stimuli$spectrum[[9]]
spectI <- spectI %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 330 == 0 ~ 330,
        round(freq, 0) %% 110 == 0 ~ 110,
        TRUE ~ NA_real_
    )
)

spectJ <- stimuli$spectrum[[10]]
spectJ <- spectJ %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 125 == 0 ~ 125,
        round(freq, 0) %% 126 == 0 ~ 126,
        TRUE ~ NA_real_
    )
)

spectK <- stimuli$spectrum[[11]]
spectK <- spectK %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 880 == 0 ~ 880,
        round(freq, 0) %% 440 == 0 ~ 440,
        TRUE ~ NA_real_
    )
)

spectL <- stimuli$spectrum[[12]]
spectL <- spectL %>% dplyr::mutate(
    freq = freq * 1000,
    fundamental = case_when(
        round(freq, 0) %% 1745 == 0 ~ 1745,
        round(freq, 0) %% 1799 == 0 ~ 1799,
        TRUE ~ NA_real_
    )
)

# Combine the labelled spectra with the rest of the data.

stimuli_labelled <- stimuli %>%
    dplyr::mutate(
        labelled_spectrum = purrr::map(stimulus, ~ {
            if (.x == "A") spectA
            else if (.x == "B") spectB
            else if (.x == "C") spectC
            else if (.x == "D") spectD
            else if (.x == "E") spectE
            else if (.x == "F") spectF
            else if (.x == "G") spectG
            else if (.x == "H") spectH
            else if (.x == "I") spectI
            else if (.x == "J") spectJ
            else if (.x == "K") spectK
            else if (.x == "L") spectL
            else NULL
        })
    )

# We have sidebands at harmonics of the modulation frequency.

identify_sidebands <- function(spectrum, freq_tol = 1, amp_rel_tol = 0.2, amp_abs_tol = 0) {
    '
    This function identifies sidebands in a spectrum.
    It searches for frequencies that are (roughly) equidistant
    from a harmonic carrier,
    and within a certain amplitude tolerance of one another.
    It returns a spectrum with the sidebands labelled.
    Sidebands are limited to ± 20-200Hz as specified in Kathios et al.
    '
    carriers <- spectrum %>% dplyr::filter(!is.na(fundamental)) %>% dplyr::pull(freq)

    non_harm <- spectrum %>% dplyr::filter(is.na(fundamental)) %>%
        dplyr::select(freq, amp)


    sideband_info <- purrr::map_dfr(carriers, function(fc) {
        if (nrow(non_harm) < 2) return(tibble::tibble())

        tmp <- non_harm %>%
            dplyr::mutate(
                delta = round(freq, 0) - round(fc, 0),
                side  = dplyr::if_else(delta < 0, "L", "R"),
                key   = round(abs(delta) / freq_tol)
            ) %>%
            dplyr::filter(key > 0)

        L <- tmp %>% dplyr::filter(side == "L") %>%
            dplyr::mutate(lower_freq = freq, lower_amp = amp, deltaL = -delta) %>% 
            dplyr::select(key, lower_freq, lower_amp, deltaL)
        R <- tmp %>% dplyr::filter(side == "R") %>%
            dplyr::mutate(upper_freq = freq, upper_amp = amp, deltaR = delta) %>% 
            dplyr::select(key, upper_freq, upper_amp, deltaR)

        pairs <- dplyr::inner_join(L, R, by = "key") %>%
            dplyr::mutate(
                delta_hz = (deltaL + deltaR) / 2,
                amp_diff = abs(lower_amp - upper_amp),
                amp_allow = pmax(amp_abs_tol, amp_rel_tol * ((lower_amp + upper_amp) / 2))
            ) %>%
            dplyr::filter(amp_diff <= amp_allow) %>%
            dplyr::filter(19 < delta_hz & delta_hz < 201) %>%
            dplyr::group_by(key) %>%
            dplyr::slice_min(amp_diff, n = 1, with_ties = FALSE) %>%
            dplyr::ungroup()

        if (nrow(pairs) == 0) return(tibble::tibble())

        pairs %>%
            dplyr::transmute(
                carrier = fc,
                lower_freq, upper_freq
            )
    })

    # Check if any sidebands were found
    if(nrow(sideband_info) == 0) {
        return(spectrum)
    }
    
    spectrum %>%
        dplyr::mutate(
            fundamental = dplyr::case_when(
                freq %in% sideband_info$lower_freq | freq %in% sideband_info$upper_freq ~ "sideband",
                TRUE ~ as.character(fundamental)
            )
        )
}

sideA <- identify_sidebands(stimuli_labelled$labelled_spectrum[[1]] %>% filter(freq >= 23))
sideB <- identify_sidebands(stimuli_labelled$labelled_spectrum[[2]] %>% filter(freq >= 23))
sideC <- identify_sidebands(stimuli_labelled$labelled_spectrum[[3]] %>% filter(freq >= 23))
sideD <- identify_sidebands(stimuli_labelled$labelled_spectrum[[4]] %>% filter(freq >= 23))
sideE <- identify_sidebands(stimuli_labelled$labelled_spectrum[[5]] %>% filter(freq >= 23))
sideF <- identify_sidebands(stimuli_labelled$labelled_spectrum[[6]] %>% filter(freq >= 23))
sideG <- identify_sidebands(stimuli_labelled$labelled_spectrum[[7]] %>% filter(freq >= 23))
sideH <- identify_sidebands(stimuli_labelled$labelled_spectrum[[8]] %>% filter(freq >= 23))
sideI <- identify_sidebands(stimuli_labelled$labelled_spectrum[[9]] %>% filter(freq >= 23))
sideJ <- identify_sidebands(stimuli_labelled$labelled_spectrum[[10]] %>% filter(freq >= 23))
sideK <- identify_sidebands(stimuli_labelled$labelled_spectrum[[11]] %>% filter(freq >= 23))
sideL <- identify_sidebands(stimuli_labelled$labelled_spectrum[[12]] %>% filter(freq >= 23))

# We now plot the spectra for visual inspection.

plot_func <- function(spectrum) {

    numeric_funds <- spectrum %>% 
            dplyr::filter(!is.na(as.numeric(fundamental))) %>% 
            dplyr::pull(fundamental) %>% 
            unique() %>% 
            as.numeric() %>% 
            sort()

    color_mapping <- c("sideband" = "green", "NA" = "grey")
    if(length(numeric_funds) == 2) {
        color_mapping[as.character(numeric_funds[1])] <- "blue"
        color_mapping[as.character(numeric_funds[2])] <- "red"
    }

    plot <- ggplot(spectrum %>% filter(freq >= 23),
     aes(x = freq, y = amp, fill = factor(fundamental))) +
        geom_col(width = max(spectrum$freq) / 200) +
        theme_minimal() +
        theme() +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        coord_cartesian(expand = FALSE) +
        scale_fill_manual(values = color_mapping, na.value = "grey",
        name = "Frequency type",
        labels = c("Harmonic (F1)", "Harmonic (F2)", "Sideband", "None"))

    return(plot)
}

# Normalise the spectrum to the maximum amplitude,
# and remove frequencies below the audible range (ish).

normalise_spectrum <- function(spectrum) {
    spectrum %>%
        dplyr::filter(freq >= 25) %>%
        dplyr::mutate(amp = amp / max(amp))
}

p1 <- plot_func(normalise_spectrum(sideA)) + ggtitle("A")
p2 <- plot_func(normalise_spectrum(sideB)) + ggtitle("B")
p3 <- plot_func(normalise_spectrum(sideC)) + ggtitle("C")
p4 <- plot_func(normalise_spectrum(sideD)) + ggtitle("D")   
p5 <- plot_func(normalise_spectrum(sideE)) + ggtitle("E")
p6 <- plot_func(normalise_spectrum(sideF)) + ggtitle("F")
p7 <- plot_func(normalise_spectrum(sideG)) + ggtitle("G")
p8 <- plot_func(normalise_spectrum(sideH)) + ggtitle("H")
p9 <- plot_func(normalise_spectrum(sideI)) + ggtitle("I")
p10 <- plot_func(normalise_spectrum(sideJ)) + ggtitle("J")
p11 <- plot_func(normalise_spectrum(sideK)) + ggtitle("K")
p12 <- plot_func(normalise_spectrum(sideL)) + ggtitle("L")

legend <- cowplot::get_legend(p1)

full <- cowplot::plot_grid(p1 + theme(legend.position = "none"),
 p2 + theme(legend.position = "none"), 
 p3 + theme(legend.position = "none"), 
 p4 + theme(legend.position = "none"), 
 p5 + theme(legend.position = "none"), 
 p6 + theme(legend.position = "none"), 
 p7 + theme(legend.position = "none"), 
 p8 + theme(legend.position = "none"), 
 p9 + theme(legend.position = "none"), 
 p10 + theme(legend.position = "none"),
 p11 + theme(legend.position = "none"),
 p12 + theme(legend.position = "none"), 
 legend,
 ncol = 4)

full_labelled <- cowplot::ggdraw(full) +
    theme(plot.margin = margin(2, 2, 2, 2, "cm")) +
    draw_label("Frequency (Hz)", x = 0.5, y = 0.2, size = 12) +
    draw_label("Amplitude", x = -0.05, y = 0.5, size = 12, angle = 90)

ggsave("outputs/frequency_analyses.png", full_labelled, width = 12, height = 12, dpi = 600)
