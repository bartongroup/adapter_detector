import os

import numpy as np
from ont_fast5_api.fast5_file import Fast5File


def mad_scaling(signal):
    '''
    Scale a signal using median absolute deviation method
    '''
    shift = np.median(signal)
    scale = np.median(np.abs(signal - shift))
    return (signal - shift) / scale


def get_fast5_fiveprime(fast5_fn, signal_size):
    '''
    Open fast5 file and return final signal_size measurements
    (corresponding to 5' end of an RNA signal). Signals
    are MAD scaled before the end is cropped.
    '''
    with Fast5File(f5_fn) as f5:
        signal = f5.get_raw_data(scale=True)
        signal = mad_scaling(signal)
        sig_len = len(signal)
        if sig_len >= size:
            fiveprime = signal[sig_len - size:]
        else:
            fiveprime = np.zeros(size)
            fiveprime[size - len(signal):] = signal
    return fiveprime