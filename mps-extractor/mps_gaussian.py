import numpy as np
from scipy.fftpack import fft, fftfreq, next_fast_len, ifft
from numpy.lib.stride_tricks import sliding_window_view
from abc import ABCMeta,abstractmethod
def gaussian_window(N, nstd):
    """
        Generate a Gaussian window of length N and standard deviation nstd.
    """
    hnwinlen = (N + (1-N%2)) // 2
    gauss_t = np.arange(-hnwinlen, hnwinlen+1, 1.0)
    gauss_std = float(N) / float(nstd)
    gauss_window = np.exp(-gauss_t**2 / (2.0*gauss_std**2)) / (gauss_std*np.sqrt(2*np.pi))
    return gauss_t,gauss_window

class ComplexSpectrumEstimator(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def estimate(self, signal, sample_rate, start_time, end_time):
        return NotImplementedError('Use a subclass of PowerSpectrumEstimator!')

    @abstractmethod
    def get_frequencies(self, signal_length, sample_rate):
        return NotImplementedError('Use a subclass of PowerSpectrumEstimator!')

class GaussianSpectrumEstimator(ComplexSpectrumEstimator):

    def __init__(self, nstd=6):
        ComplexSpectrumEstimator.__init__(self)
        self.nstd = nstd
        self._gauss_window_cache = {}

    def get_frequencies(self, signal_length, sample_rate):
        signal_length = next_fast_len(signal_length)
        freq = fftfreq(signal_length, d=1.0/sample_rate)
        nz = freq >= 0.0
        return freq[nz]

    def _get_gauss_window(self, nwinlen):
        if nwinlen in self._gauss_window_cache:
            return self._gauss_window_cache[nwinlen]
        else:
            if nwinlen % 2 == 0:
                nwinlen += 1
            hnwinlen = nwinlen // 2
            gauss_t = np.arange(-hnwinlen, hnwinlen+1, 1.0)
            gauss_std = float(nwinlen) / float(self.nstd)
            gauss_window = np.exp(-gauss_t**2 / (2.0*gauss_std**2)) / (gauss_std*np.sqrt(2*np.pi))
            self._gauss_window_cache[nwinlen] = gauss_window
            return gauss_window

    def estimate(self, signal, sample_rate, start_time, end_time):
        nwinlen = len(signal)
        gauss_window = self._get_gauss_window(nwinlen)

        fft_len = next_fast_len(len(signal))
        #window the signal and take the FFT
        windowed_slice = signal[:fft_len]*gauss_window[:fft_len]
        s_fft = fft(windowed_slice, n=fft_len, overwrite_x=1)
        freq = fftfreq(fft_len, d=1.0/sample_rate)
        nz = freq >= 0.0

        return freq[nz],s_fft[nz]

def timefreq(s, sample_rate, window_length, increment, spectrum_estimator, min_freq=0, max_freq=None, zero_pad=True):
    """
        Compute a time-frequency representation of the signal s.

        s: the raw waveform.
        sample_rate: the sample rate of the waveform
        increment: the spacing in seconds between points where the spectrum is computed, i.e. inverse of the spectrogram sample rate
        spectrum_estimator: an instance of PowerSpectrumEstimator
        min_freq: the minimum frequency to analyze (Hz)
        max_freq: the maximum frequency to analyze (Hz)
        window_length: The length in seconds of the window to use for each segment. If None, then the spectrum_estimator
            object is queried.

        Returns t,freq,spec,rms:

        t: the time axis of the spectrogram
        freq: the frequency axis of the spectrogram
        tf: the time-frequency representation
    """

    if max_freq is None:
        max_freq = sample_rate / 2.0

    if window_length is None:
        window_length = spectrum_estimator.get_window_length()

    #compute lengths in # of samples
    nwinlen = int(sample_rate*window_length)
    if nwinlen % 2 == 0:
        nwinlen += 1
    hnwinlen = nwinlen // 2
    assert len(s) > nwinlen, "len(s)=%d, nwinlen=%d" % (len(s), nwinlen)

    # get the values for the frequency axis by estimating the spectrum of a dummy slice
    full_freq = spectrum_estimator.get_frequencies(nwinlen, sample_rate)
    freq_index = (full_freq >= min_freq) & (full_freq <= max_freq)
    freq = full_freq[freq_index]
    nfreq = freq_index.sum()

    nincrement = int(np.round(sample_rate*increment))
    if zero_pad:
        # print 'len(s)=%d, nwinlen=%d, hwinlen=%d, nincrement=%d, nwindows=%d' % (len(s), nwinlen, hnwinlen, nincrement, nwindows)
        #pad the signal with zeros
        zs = np.zeros([len(s) + 2*hnwinlen])
        zs[hnwinlen:-hnwinlen] = s
        windows = sliding_window_view(zs, nwinlen, axis=0)[::nincrement]
        window_centers = windows[:, hnwinlen]
        nwindows = len(windows)
    else:
        windows = sliding_window_view(s, nwinlen, axis=0)[::nincrement]
        window_centers = windows[:, hnwinlen]
        nwindows = len(window_centers)
        assert nwindows > 0, "nwindows=0, len(s)=%d, nwinlen=%d, nincrement=%d, window_centers=%s" % (len(s), nwinlen, nincrement, str(window_centers))
        zs = s
        assert window_centers.min() >= hnwinlen, "window_centers.minmax=(%d,%d), hnwinlen=%d, len(s)=%d" % (window_centers.min(), window_centers.max(), hnwinlen, len(s))
        assert window_centers.max() < len(s)-hnwinlen, "window_centers.minmax=(%d,%d), hnwinlen=%d, len(s)=%d" % (window_centers.min(), window_centers.max(), hnwinlen, len(s))

    #take the FFT of each segment, padding with zeros when necessary to keep window length the same
    #tf = np.zeros([nfreq, nwindows], dtype='complex')
    tf = np.zeros([nfreq, nwindows], dtype='complex')
    for k, window in enumerate(windows):
        si = window_centers[k] - hnwinlen
        ei = window_centers[k] + hnwinlen + 1
        spec_freq,est = spectrum_estimator.estimate(window, sample_rate, si/sample_rate, ei/sample_rate)
        findex = (spec_freq <= max_freq) & (spec_freq >= min_freq)
        #print 'k=%d' % k
        #print 'si=%d, ei=%d' % (si, ei)
        #print 'spec_freq.shape=',spec_freq.shape
        #print 'tf.shape=',tf.shape
        #print 'est.shape=',est.shape
        tf[:, k] = est[findex]

    # Note that the desired spectrogram rate could be slightly modified
    t = np.arange(0, nwindows, 1.0) * float(nincrement)/sample_rate

    return t, freq, tf

def gaussian_stft(s, sample_rate, window_length, increment, min_freq=0,
                  max_freq=None, nstd=6, zero_pad=True):
    """
        Compute a gaussian-windowed short-time fourier transform representation of the signal s.

        s: the raw waveform.
        sample_rate: the sample rate of the waveform
        window_length: The length in seconds of the window to use for each segment.
        increment: the spacing in seconds between points where the spectrum is computed, i.e. inverse of the spectrogram sample rate
        min_freq: the minimum frequency to analyze (Hz)
        max_freq: the maximum frequency to analyze (Hz) (if None then the nyquist frequency will be used)

        Returns t,freq,spec,rms:

        t: the time axis of the spectrogram
        freq: the frequency axis of the spectrogram
        tf: the time-frequency representation
    """

    spectrum_estimator = GaussianSpectrumEstimator(nstd=nstd)
    t, freq, tf = timefreq(s, sample_rate, window_length, increment,
                           spectrum_estimator=spectrum_estimator,
                           min_freq=min_freq, max_freq=max_freq,
                           zero_pad=zero_pad)
    ps = np.abs(tf)
    rms = ps.sum(axis=0)
    return t, freq, tf, rms