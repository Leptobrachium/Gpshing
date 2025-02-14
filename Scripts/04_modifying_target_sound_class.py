###########################################################################
# Cane Toad Training Data Modifying
#
# This script used to apply a low-pass filter to the
# cane toad training audio data. Filtering out frequencies above 1300 Hz helps 
# reduce interference from other amphibians and insect noise, enhancing the 
# classifier's robustness.
#
# References:
# - MacCallum et al. (2011)
#
# Dependencies: os, numpy, scipy.signal, soundfile
#
# Note: Update file paths and parameter values as appropriate.
###########################################################################


import os
import numpy as np
from scipy.signal import butter, filtfilt
import soundfile as sf

# Define the cutoff frequency (Hz)
cutoff_frequency = 1300

# Define the Butterworth filter order
filter_order = 3

# Define a reduction factor for frequencies above the cutoff (if desired)
reduction_factor = 0.7  # Adjust as needed

# Path to the directory containing the audio files (assumed to be in .wav format)
audio_dir = r"path/to/wav"

# Create a new directory for saving modified (filtered) audio files
save_dir = os.path.join(audio_dir, "modified_Cane_toad")
os.makedirs(save_dir, exist_ok=True)

# Get a list of WAV audio files in the directory
wav_files = [file for file in os.listdir(audio_dir) if file.endswith('.wav')]

# Loop through each audio file
for file_name in wav_files:
    print(f"Processing file: {file_name}")
    file_path = os.path.join(audio_dir, file_name)

    # Load the audio file (audio data and sampling rate)
    audio, sr = sf.read(file_path)

    # Normalize the audio signal to the range [-1, 1]
    audio = audio / np.max(np.abs(audio))

    # Compute the Nyquist frequency (half the sampling rate)
    nyquist = 0.5 * sr

    # Compute the normalized cutoff frequency for the filter
    cutoff_norm = cutoff_frequency / nyquist

    # Design the Butterworth low-pass filter coefficients
    b, a = butter(filter_order, cutoff_norm, btype='low')

    # Apply the low-pass filter to the audio signal
    filtered_audio = filtfilt(b, a, audio)

    # Optionally, reduce the amplitude of signal components above the cutoff
    filtered_audio[audio > cutoff_frequency] *= reduction_factor

    # Create the modified file name and path
    modified_file_name = 'modified_' + file_name
    modified_file_path = os.path.join(save_dir, modified_file_name)

    # Save the filtered audio file in WAV format
    sf.write(modified_file_path, filtered_audio, sr, format='WAV')
    print(f"Saved modified file as: {modified_file_path}")