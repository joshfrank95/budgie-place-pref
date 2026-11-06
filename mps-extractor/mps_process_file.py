import librosa
import numpy as np
import csv
from mps_core import (
    spectrogram,
    mps,
    strengths_by_temporal_bin,
    mps_roughness,
    # window_rms,
)
def process_file(
    wav_path,
    tmr_low,
    tmr_high,
    rms_weighted,
    return_percentage=False,
    mode="roughness",
    normalize=True,
    Norm=True,
    sr_target=44100,
    window_duration=0.50,
    window_hop_size=0.25,
    spec_sample_rate=1000,
    freq_spacing=50,
    min_freq=0,
    max_freq=20000,
    nstd=6,
    cmplx=True,
    mps_window=0.50,
    # wav_path,
    # win_sec,
    # hop_sec,
):

    rows = []
    weighted_candidates = []  # for roughness
    tm_vectors = []  # list of (col_sums, rms, start, end)
    

    try:
        samples, sample_rate = librosa.load(wav_path, sr = None, mono = True)

        if sample_rate != sr_target:
            # print(f"  Resampling {wav_path} from {sample_rate} Hz to {sr_target} Hz...")
            samples = librosa.resample(y=samples, orig_sr=sample_rate, target_sr=sr_target)
            sample_rate = sr_target
            # print(f"  Resampling complete.")

        original_duration = len(samples) / float(sample_rate)
        # print(f"  Original duration (after resampling): {original_duration:.3f}s")

    except Exception as e:
        print(f"ERROR loading/resampling {wav_path}: {e}")
        # Optionally create an empty/error CSV here
        return

    if normalize:
        maxAmp = np.abs(samples).max()
        if maxAmp == 0:
            print(f"  Warning: Max amplitude is zero for {wav_path}. Normalization will have no effect.")
            maxAmp = 1.0 # Avoid division by zero later, signal is already zeros
    else:
        maxAmp = 1.0 # Set to 1.0 if not normalizing

    total_samples = len(samples) # Use original length here
    sample_rate = float(sample_rate)
    window_size = int(window_duration * sample_rate)

    if window_size > total_samples:
        window_size = total_samples
        print(f"  Warning: Original signal length ({total_samples} samples) is shorter than window duration ({window_duration}s). Using full length as window.")
        hop_size_samples = 0
        num_windows = 1
    elif window_size == 0:
        print(f"  Error: Calculated window size is 0 for {wav_path}. Skipping multi-window processing.")
        rows.append((str(wav_path), 'N/A', 'N/A', 'Window Size Error'))
        num_windows = 0 # Prevent loop from running
    else:
        hop_size_samples = int(window_hop_size * sample_rate)
        hop_size_samples = max(1, hop_size_samples)
        if total_samples < window_size:
            num_windows = 1
            hop_size_samples = 0
        else:
            num_windows = 1 + int(np.floor((total_samples - window_size) / hop_size_samples))

    # print(f"  Calculated num_windows = {num_windows} (Window Size: {window_size}, Hop Size: {hop_size_samples})")

    for window_num in range(num_windows):
                    start_idx = window_num * hop_size_samples
                    end_idx = start_idx + window_size
                    end_idx = min(end_idx, total_samples)
                    audio_window = samples[start_idx:end_idx] # Extract from ORIGINAL signal
                    start_time = start_idx / sample_rate
                    end_time = end_idx / sample_rate

                    # Extract window data from potentially padded 'samples'
                    window_data = samples[start_idx:end_idx]

                    # Check for unexpectedly short windows (can happen at the very end if total_samples is not multiple of hop_size)
                    # We process these potentially shorter final windows.
                    # if end_idx - start_idx < window_size and window_num != num_windows -1 :
                    #      print(f"  Skipping unexpectedly short window {window_num + 1} (indices {start_idx}-{end_idx}, size {end_idx-start_idx}, expected {window_size})")
                    #      continue

                    # --- Apply Normalization Per Window ---
                    if normalize and maxAmp > 0: # Use maxAmp calculated from ORIGINAL signal
                        signal_data = window_data.astype(float) / maxAmp
                    else:
                        signal_data = window_data.astype(float)
                    # --- End Normalization ---

                    # print(f"\n--- Processing Window {window_num + 1}/{num_windows} for {wav_path} ---")
                    # print(f"  Sample indices: {start_idx} to {end_idx}")

                    # Calculate RMS for the window
                    rms_value = np.sqrt(np.mean(signal_data**2))
                    silence_threshold = 1e-6

                    if rms_value < silence_threshold:
                        # print(f"  Skipping window {window_num+1} - Near-silent audio detected (RMS: {rms_value:.6f})")
                        rows.append((str(wav_path), start_time, end_time, 'Silent Window'))
                        continue # Skip to the next window

                    # --- Process this audio window ---
                    try:
                        # print(f"--- Calculating Spectrogram for window {window_num + 1} ---")
                        # Use signal_data (normalized window)
                        t, f, spec, rms_spec = spectrogram(signal_data, sample_rate, spec_sample_rate, freq_spacing, min_freq, max_freq, nstd, cmplx)

                        # --- Log Conversion & Cleaning ---
                        # print(f"\n--- DEBUG Log Spectrogram [Window {window_num + 1}] ---")
                        spec_abs = np.abs(spec)
                        epsilon = 1e-12
                        with np.errstate(divide='ignore'):
                            spec_log = 20 * np.log10(spec_abs + epsilon)

                        spec_log_neginfs = np.isneginf(spec_log).sum()
                        if spec_log_neginfs > 0:
                            finite_vals = spec_log[np.isfinite(spec_log)]
                            min_finite_val = np.min(finite_vals) if finite_vals.size > 0 else -100
                            floor_db = min(min_finite_val - 50, -300)
                            # print(f"  Replacing {spec_log_neginfs} -Inf values with {floor_db:.2f} dB.")
                            spec_log[np.isneginf(spec_log)] = floor_db
                        # --- End Cleaning ---


                        # --- MPS Calculation & Debug ---
                        # print(f"\n--- DEBUG MPS Input [Window {window_num + 1}] ---")
                        if np.isnan(spec_log).any() or np.isinf(spec_log).any():
                            # print(f"  ERROR: NaNs or Infs in spec_log! Skipping MPS for window {window_num + 1}.")
                            rows.append((str(wav_path), start_time, end_time, 'Pre-MPS NaN/Inf Error'))
                            continue # Skip to next window

                        spec_time_step = (t[1] - t[0]) if len(t) > 1 else (1.0 / spec_sample_rate)
                        mps_window_points = int(round(mps_window / spec_time_step))
                        if mps_window_points % 2 == 0: mps_window_points += 1

                        if spec_log.shape[1] < mps_window_points:
                            # print(f"  ERROR: Spectrogram time points ({spec_log.shape[1]}) < required MPS window points ({mps_window_points}). Skipping MPS for window {window_num + 1}.")
                            rows.append((str(wav_path), start_time, end_time, 'Spectrogram Too Short'))
                            continue # Skip to next window

                        wf, canonical_wt, mps_powAvg = mps(spec_log, f, t, window=mps_window, Norm=Norm, return_percentage=return_percentage)


                        # --- MPS Validity Check ---
                        mps_valid = False
                        if mps_powAvg is not None and isinstance(mps_powAvg, np.ndarray) and mps_powAvg.size > 0:
                            if not np.isnan(mps_powAvg).any() and not np.isinf(mps_powAvg).any():
                                mps_valid = True

                        # --- Get colsums ---
                        col_sums = strengths_by_temporal_bin(mps_powAvg)

                        # # --- Get roughness ---
                        # rough = mps_roughness(col_sums, canonical_wt, tmr_low, tmr_high, weights=None, rms_weighted = False)

                        # # --- Return results ---
                        # rows.append((str(wav_path), start_time, end_time, float(rough)))

                        if mode == "roughness":
                            rough = mps_roughness(
                                 col_sums,
                                 canonical_wt,
                                 tmr_low,
                                 tmr_high,
                                 weights=None
                            )
                            if rms_weighted:
                                 weighted_candidates.append(
                                     (str(wav_path), start_time, end_time, float(rough), rms_value)
                                 )
                            else:
                                 rows.append((str(wav_path), start_time, end_time, float(rough)))
                        elif mode.startswith("tmr_"):
                            # store vector and rms for later weighting logic
                            tm_vectors.append((col_sums, rms_value, start_time, end_time))


                    except Exception as e:
                        print(f"ERROR processing window {window_num + 1} for {wav_path}: {e}")
                        rows.append((str(wav_path), start_time, end_time, 'Error Processing Window'))
                        continue # Skip to next window

    # --------------------------------------------------------------
    # Apply normalised RMS weighting exactly once, after all windows
    # --------------------------------------------------------------
    if rms_weighted and weighted_candidates:
        max_rms = max(r for *_ , r in weighted_candidates)
        if max_rms == 0:
            max_rms = 1.0  # avoid divide-by-zero (all-silent edge case)

        for f, st, et, ro, rv in weighted_candidates:
            rows.append((f, st, et, ro * (rv / max_rms)))


        # ------------------------
    # Handle tmr_summary / tmr_file output after all windows
    if mode.startswith("tmr_") and tm_vectors:
        # determine RMS weighting
        max_rms = max(rms for _, rms, _, _ in tm_vectors)
        vecs = []
        for v, rms_v, st, et in tm_vectors:
            weight = (rms_v / max_rms) if rms_weighted else 1.0
            vecs.append((v * weight, st, et))
            if mode == "tmr_summary":
                rows.append((str(wav_path), st, et, canonical_wt.tolist(), (v * weight).tolist()))
        if mode == "tmr_file":
            mean_vec = np.mean(np.stack([vv for vv, *_ in vecs]), axis=0)
            rows.append((str(wav_path), len(vecs), canonical_wt.tolist(), mean_vec.tolist()))

    return rows
