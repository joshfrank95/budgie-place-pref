import numpy as np
import librosa
from mps_gaussian import gaussian_stft   # new file holding the ported code
from copy import deepcopy
from numpy.fft import fftshift, ifftshift
from scipy.fftpack import fft, ifft, fftfreq, fft2, ifft2
from mps_gaussian import (gaussian_window, gaussian_stft)
# ---------- helpers --------------------------------------------------

def spectrogram(s, sample_rate, spec_sample_rate = 1000, freq_spacing = 50, min_freq=0, max_freq=20000, nstd=6,  cmplx = True):
    """
        Given a sound pressure waveform, s, compute the complex spectrogram. 
        See documentation on gaussian_stft for additional details.

        Returns: 
          t, freq, timefreq, rms
              t: array of time values to use as x axis
              freq: array of frequencies to use as y axis
              timefreq: the spectrogram (a time-frequency represnetaion)
              rms : the time varying average
              
        Arguments:
            REQUIRED:
                s: sound pressssure waveform
                sample_rate: sampling rate for s in Hz
                spec_sample_rate: sampling rate for the output spectrogram in Hz. This variable sets the overlap for the windows in the STFFT.
                freq_spacing: the time-frequency scale for the spectrogram in Hz. This variable determines the width of the gaussian window. 
            
            OPTIONAL            
                complex = False: returns the absolute value
                use min_freq and max_freq to save space
                nstd = number of standard deviations of the gaussian in one window.
    """
     
    # We need units here!!
    increment = 1.0 / spec_sample_rate
    window_length = nstd / (2.0*np.pi*freq_spacing)
    t,freq,timefreq,rms = gaussian_stft(s, sample_rate, window_length, increment, nstd=nstd, min_freq=min_freq, max_freq=max_freq)

    # rms = spec.std(axis=0, ddof=1)
    if cmplx == False:
        timefreq = np.abs(timefreq)
        
    return t, freq, timefreq, rms


def mtfft(spectrogram, df, dt, Log=False, return_percentage=False):
    """
        Compute the 2d modulation power and phase for a given time frequency slice.
        return temporal_freq,spectral_freq,mps_pow,mps_phase
    """

    #take the 2D FFT and center it
    smps = fft2(spectrogram)
    smps = fftshift(smps)

    #compute the log amplitude
    mps_pow = np.abs(smps)**2
    # if Log:
    #     mps_pow = 10*np.log10(mps_pow)
    if return_percentage:
        total_power = np.sum(mps_pow)
        if total_power > 0:
            mps_pow = 100.0 * mps_pow / total_power

    #compute the phase
    mps_phase = np.angle(smps)

    #compute the axes
    nf = mps_pow.shape[0]
    nt = mps_pow.shape[1]

    spectral_freq = fftshift(fftfreq(nf, d=df[1]-df[0]))
    temporal_freq = fftshift(fftfreq(nt, d=dt[1]-dt[0]))

    """
    nb = sdata.shape[1]
    dwf = np.zeros(nb)
    for ib in range(int(np.ceil((nb+1)/2.0))+1):
        posindx = ib
        negindx = nb-ib+2
        print 'ib=%d, posindx=%d, negindx=%d'% (ib, posindx, negindx )
        dwf[ib]= (ib-1)*(1.0/(df*nb))
        if ib > 1:
            dwf[negindx] =- dwf[ib]

    nt = sdata.shape[0]
    dwt = np.zeros(nt)
    for it in range(0, int(np.ceil((nt+1)/2.0))+1):
        posindx = it
        negindx = nt-it+2
        print 'it=%d, posindx=%d, negindx=%d' % (it, posindx, negindx)
        dwt[it] = (it-1)*(1.0/(nt*dt))
        if it > 1 :
            dwt[negindx] = -dwt[it]

    spectral_freq = dwf
    temporal_freq = dwt
    """

    return spectral_freq, temporal_freq, mps_pow, mps_phase

def mps(spectrogram, df, dt, window=None, Norm=True, tShift_overlap=10, return_percentage=False):
    """
    Calculates the modulation power spectrum using overlapp and add method with a gaussian window of length window in s
    Assumes that spectrogram is in dB.  df and dt are the axis of spectrogram.
    """
    
    # Debugging paramenter
    debugPlot = False

    # Resolution of spectrogram in DB
    dbRES = 50 
    
    # Check the size of the spectrogram vs dt
    nt = dt.size
    nf = df.size
    if spectrogram.shape[1] != nt and spectrogram.shape[0] != nf:   
        # print('Error in mps. Expected  %d bands in frequency and %d points in time' % (nf, nt))
        # print('Spectrogram had shape %d, %d' % spectrogram.shape)
        return 0, 0, 0
        
    # Z-score the flattened spectrogram is Norm is True
    sdata = deepcopy(spectrogram)

    maxdata = sdata.max()
    mindata = maxdata - dbRES
    sdata[sdata < mindata] = mindata
    
    if Norm:
        # maxdata = sdata.max()
        # mindata = maxdata - dbRES
        # sdata[sdata< mindata] = mindata
        # print(f"  Min value of sdata: {sdata.min():.3f}")
        # print(f"  Max value of sdata: {sdata.max():.3f}")
        sdata -= sdata.mean()
        sdata /= sdata.std()
        
    if window == None:
        window = dt[-1]/10.0 # Default window size if none provided

    # --- Robust nWindow and dt_step calculation ---
    if nt > 1:
        # Calculate actual step from input dt, assuming it's mostly consistent
        d_step = dt[1] - dt[0]
    elif nt == 1:
         # Cannot calculate step, MPS not meaningful
         print("Warning in mps: Spectrogram has only one time point.")
         return np.asarray([]), np.asarray([]), np.asarray([])
    else: # nt == 0
         print("Warning in mps: Spectrogram has zero time points.")
         return np.asarray([]), np.asarray([]), np.asarray([])

    # Check if requested window is longer than total duration
    total_duration = dt[-1] - dt[0] + d_step # Add one step for total duration
    if window > total_duration:
        print('Warning in mps: Requested MPS window size (%.3f s) is greater than spectrogram temporal extent (%.3f s).' % (window, total_duration))
        print('mps will be calculated with windows of %d points or %.3f s' % (nt, total_duration) )
        nWindow = nt # Use all available points
    else:
        # Calculate nWindow based on duration and step, using round() for consistency
        nWindow_float = window / d_step
        nWindow = int(round(nWindow_float))

    # Ensure nWindow is at least 1 and make it odd size for a symmetric window
    nWindow = max(1, nWindow)
    if nWindow % 2 == 0:
        nWindow += 1

    pad_size = (nWindow - 1) // 2 # Still need pad_size for indexing
    # --- End of robust calculation ---

    # --- Calculate canonical wt once ---
    # Assuming mtfft uses fftfreq like this:
    # Note: If mtfft uses a different convention (e.g., different d), adjust this.
    canonical_wt = fftshift(fftfreq(nWindow, d=d_step))
    canonical_wf = None # wf depends on df, calculate inside loop or from first chunk?
                       # Let's get it from the first valid mtfft call.
    # ---

    # Generate the Gaussian window using the calculated nWindow
    alpha = 6 # Standard deviations for Gaussian window
    gt, wg = gaussian_window(nWindow, alpha)
    
    # Determine chunk overlap/shift (e.g., 1/3 of the window)
    # Ensure shift is at least 1 point
    tShift = max(1, int(nWindow / tShift_overlap))
    
    # Initialize average MPS power spectrum
    # Need the shape, get it from a dummy mtfft call or assume based on df and nWindow
    # Let's initialize after the first successful chunk instead.
    mps_powAvg = None
    nchunks = 0

    # Pad the spectrogram data temporally to handle windowing at edges
    # pad_size = (nWindow - 1) // 2
    # padded_sdata = np.pad(sdata, ((0, 0), (pad_size, pad_size)), mode='edge') # Pad time axis

    # Loop through the *original* time points as centers for the windows
    for tmid_idx in range(pad_size, nt - pad_size, tShift):
        
        # # Calculate indices in the *padded* data array
        # tstart_padded = tmid_idx # Center in original maps to start of window in padded
        # tend_padded = tstart_padded + nWindow

        # Calculate indices directly in the *original* sdata array
        tstart = tmid_idx - pad_size
        tend = tmid_idx + pad_size + 1 # Slice goes up to, but not including, tend
        
        # # Ensure indices are within the bounds of the padded array
        # if tend_padded > padded_sdata.shape[1]:
        #      # This shouldn't happen if padding is correct, but check just in case
        #      print(f"Warning: Skipping chunk at tmid_idx {tmid_idx} due to boundary issue.")
        #      continue

        # Extract the windowed spectrogram segment from padded data
        # wSpect = deepcopy(padded_sdata[:, tstart_padded:tend_padded])

        # Extract the windowed spectrogram segment from original data
        # No need for deepcopy here if we apply window carefully
        wSpect_segment = sdata[:, tstart:tend]

        # Check if the extracted window has the correct size
        if wSpect_segment.shape[1] != nWindow:
            print(f"Error: Window size mismatch. Expected {nWindow}, got {wSpect_segment.shape[1]} at tmid_idx {tmid_idx}. Skipping chunk.")
            continue # Skip this chunk

        # Apply the Gaussian window (create a temporary copy to avoid modifying sdata)
        wSpect = wSpect_segment * wg[np.newaxis, :] # Apply window using broadcasting

        nchunks += 1 # Increment chunk count *after* validation

        # # Multiply the spectrogram segment by the Gaussian window
        # # Ensure wg has the same length as wSpect's time dimension
        # if len(wg) != wSpect.shape[1]:
        #      # This should not happen if nWindow calculation and slicing are correct
        #      print(f"Error: Gaussian window size {len(wg)} doesn't match spectrogram segment size {wSpect.shape[1]}. Skipping chunk.")
        #      nchunks -= 1 # Decrement count as this chunk is invalid
        #      continue

        # for fInd in range(nf):
        #     wSpect[fInd,:] = wSpect[fInd,:] * wg # Apply window row by row

        # Get the 2d FFT using mtfft
        # Pass the consistent d_step instead of a dt slice if mtfft supports it.
        # Assuming mtfft needs the dt slice for now:
        # We need a time vector corresponding to wSpect, but mtfft only uses its step.
        # Let's pass the original dt slice corresponding to the *center* of the window.
        # This is complex. A better mtfft would take nt and d_step.
        # Sticking with original approach for now: pass dt slice of length nWindow
        # This slice is only used by mtfft to determine the step 'd'.

        dt_slice = dt[tstart:tend]
        # Save wSpect to CSV
        # csv_filename = os.path.join('h5files', '_wSpect_chunk.csv')
        # np.savetxt(csv_filename, wSpect, delimiter=',')
        wf_chunk, wt_chunk, mps_pow, mps_phase = mtfft(wSpect, df, dt_slice, return_percentage=return_percentage) # Pass dummy dt slice
        # Save mps_pow to CSV
        # csv_filename = os.path.join('h5files', '_mps_pow_chunk{nchunks}.csv')
        # np.savetxt(csv_filename, mps_pow, delimiter=',')

        # Initialize or add to the average power spectrum
        if mps_powAvg is None:
            mps_powAvg = np.zeros_like(mps_pow)
            # Get the spectral frequency axis from the first valid chunk
            canonical_wf = wf_chunk
            # Basic check: Ensure wt from this chunk matches canonical wt
            if not np.allclose(wt_chunk, canonical_wt):
                 print(f"Warning: wt mismatch detected in chunk {nchunks}. Using canonical wt.")
                 # This indicates an issue either here or in mtfft

        # Check if mtfft returned valid dimensions matching mps_powAvg
        if mps_pow.shape != mps_powAvg.shape:
             print(f"Error: Shape mismatch between mps_pow {mps_pow.shape} and mps_powAvg {mps_powAvg.shape}. Skipping chunk.")
             # If this was the first chunk, reset mps_powAvg
             if nchunks == 1:
                 mps_powAvg = None
                 canonical_wf = None
             nchunks -= 1 # Decrement count as this chunk is invalid
             continue

        mps_powAvg += mps_pow # Add subsequent valid chunks

    # Average the MPS power, handle division by zero if no chunks were processed
    if nchunks > 0 and mps_powAvg is not None:
        mps_powAvg /= nchunks
    elif mps_powAvg is None: # No valid chunks processed
        print("Warning: No valid chunks processed for MPS calculation.")
        # Return empty arrays or zeros matching expected output structure
        return np.asarray([]), np.asarray([]), np.asarray([]) # Return empty arrays

    # Return the final averaged MPS and the *canonical* axes
    if canonical_wf is None: # Should not happen if nchunks > 0
         print("Error: canonical_wf not set despite processed chunks.")
         return np.asarray([]), np.asarray([]), np.asarray([])

    return canonical_wf, canonical_wt, mps_powAvg

def get_mps(t, freq, spec):
    "Computes the MPS of a spectrogram (idealy a log-spectrogram) or other REAL time-freq representation"
    mps = fftshift(fft2(spec))
    amps = np.real(mps * np.conj(mps))
    nf = mps.shape[0]
    nt = mps.shape[1]
    wfreq = fftshift(fftfreq(nf, d=freq[1] - freq[0]))
    wt = fftshift(fftfreq(nt, d=t[1] - t[0]))
    return wt, wfreq, mps, amps


def strengths_by_temporal_bin(mps_pow):
    """Sum over spectral axis → 1-D vector (temporal bins)."""
    return mps_pow.sum(axis=0)


def mps_roughness(col_sums, wt_bins, tmr_low=11, tmr_high=65, weights=None):
    """
    Weighted (or un-weighted) energy between |wt| in [tmr_low, tmr_high] Hz.
    """
    mask = (np.abs(wt_bins) >= tmr_low) & (np.abs(wt_bins) <= tmr_high)
    vec  = col_sums[mask]
    if weights is not None:
        vec = vec * weights
    return float(vec.sum())


# def window_rms(x):
#     return float(np.sqrt(np.mean(x.astype(float) ** 2)))