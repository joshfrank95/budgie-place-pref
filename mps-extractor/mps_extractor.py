#!/usr/bin/env python3
"""
Compute MPS-based roughness per window for *all* audio files in a folder
and save one CSV with:

file_path,start_time,end_time,rms,mps_roughness
"""
from pathlib import Path
from argparse import ArgumentParser
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd
import soundfile as sf

from mps_core import (
    spectrogram,
    mps,
    strengths_by_temporal_bin,
    mps_roughness,
    # window_rms,
)
from mps_process_file import process_file

# ---------------------------------------------------------------------



# ---------------------------------------------------------------------
def main():
    ap = ArgumentParser()
    ap.add_argument("--audio_dir", required=True, help="Folder with wav/mp3 files")
    ap.add_argument("--csv_out", required=True, help="Output CSV path")
    ap.add_argument("--window", type=float, default=0.50, help="window length (s)")
    ap.add_argument("--hop", type=float, default=0.25, help="hop length (s)")
    ap.add_argument("--sr", type=int, default=44100, help="target sample-rate (Hz)")
    ap.add_argument("--tmr_low", type=float, default=11, help="low |wt| cutoff (Hz)")
    ap.add_argument("--tmr_high", type=float, default=65, help="high |wt| cutoff (Hz)")
    ap.add_argument("--rms_weighted", action="store_true", help="weight by window RMS")
    ap.add_argument("--return_percentage", action="store_true", help="Return MPS power as percentage of total spectrogram energy")
    ap.add_argument("--mode", choices=["roughness", "tmr_summary", "tmr_file"], default="roughness", help="Output mode: roughness scalar or temporal-modulation summary")
    ap.add_argument("--mps_window", type=float, default=0.50, help="MPS window length (s)")
    ap.add_argument(
        "--n_jobs", type=int, default=-1, help="Joblib worker count (-1 = all cores)"
    )
    args = ap.parse_args()

    audio_paths = sorted(
        p
        for p in Path(args.audio_dir).rglob("*")
        if p.suffix.lower() in {".wav", ".mp3"}
    )
    if not audio_paths:
        ap.error("No .wav or .mp3 files found")

    proc = delayed(process_file)

    results = Parallel(args.n_jobs)(
        proc(
            wav_path=p,
            window_duration=args.window,
            window_hop_size=args.hop,
            sr_target=args.sr,
            tmr_low=args.tmr_low,
            tmr_high=args.tmr_high,
            rms_weighted=args.rms_weighted,
            mps_window=args.mps_window,
            return_percentage=args.return_percentage,
            mode=args.mode,
        )
        for p in tqdm(audio_paths, desc="Audio")
    )

    flat = [row for sub in results for row in sub]
    if args.mode == "roughness":
        columns = ["file", "start_time", "end_time", "mps_roughness"]
    elif args.mode == "tmr_summary":
        columns = ["file", "start_time", "end_time", "wt_bins", "energies"]
    else:  # tmr_file
        columns = ["file", "n_windows", "wt_bins", "mean_energies"]

    out_df = pd.DataFrame(flat, columns=columns)
    if args.mode == "roughness":
        out_df["mps_roughness"] = pd.to_numeric(
            out_df["mps_roughness"], errors="coerce"
        )
    out_df.to_csv(args.csv_out, index=False)
    print(f"Wrote {len(out_df)} rows → {args.csv_out}")


if __name__ == "__main__":
    import librosa  # deferred heavy import
    main()