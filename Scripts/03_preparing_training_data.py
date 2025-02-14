###########################################################################
# Cane Toad Training Dataset Preparation 
#
# This script used to sgement labeled audio files 
# into three-second audio segments for use with BirdNET. The training 
# dataset comprises both true positives (cane toad calls) and repeated 
# false-positive sounds (see Table 1 in the manuscript).
#
# References:
# - Hu and Wang (2007) for the AudioSegment module.
#
# Dependencies: os, pandas, pydub
#
# Note: Update file paths and parameter values with those appropriate for your data.
###########################################################################

import os
import pandas as pd
from pydub import AudioSegment

def split_audio(input_audio, trim_list, output_folder):
    """
    Splits an audio file into three-second segments based on provided trim information.

    Parameters:
    - input_audio: str, path to the input audio file.
    - trim_list: list of dictionaries, each containing:
          - "start_time": float, start time in seconds.
          - "file_name": str, base name for the output file.
          - "name": str, additional identifier for the output file.
    - output_folder: str, directory where the output audio segments will be saved.
    """
    # Load the audio file using pydub
    audio = AudioSegment.from_file(input_audio)
    
    for trim_info in trim_list:
        start_time_sec = trim_info["start_time"]
        # Define a three-second segment (start_time to start_time + 3)
        end_time_sec = start_time_sec + 3  
        file_name = trim_info["file_name"]
        name = trim_info["name"]
        # Construct output file name by combining file_name and name
        output_audio_name = f"{file_name}_{name}.wav"
        output_audio = os.path.join(output_folder, output_audio_name)
        # Extract the segment (convert seconds to milliseconds)
        trim_audio = audio[start_time_sec * 1000:end_time_sec * 1000]
        # Export the trimmed segment as a WAV file
        trim_audio.export(output_audio, format="wav")

# Update these paths as needed
csv_file_path = "path/to/csv"
audio_folder_path = "path/to/Canetoad"
output_folder_path = "path/to/Canetoad_Ouput_wav"

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder_path):
    os.makedirs(output_folder_path)

# Read the CSV file containing segmentation metadata
df = pd.read_csv(csv_file_path)


# Iterate through each row in the CSV file
for index, row in df.iterrows():
    file_name = row["file_name"]
    name = row["name"]
    # Construct the input audio file path (assumes .wav format)
    input_audio_path = os.path.join(audio_folder_path, f"{file_name}.wav")
    
    # Prepare trim information as a dictionary (wrapped in a list)
    trim_info = [row.to_dict()]
    
    # Split the audio file and save the three-second segment
    split_audio(input_audio_path, trim_info, output_folder_path)

print("Audio splitting completed.")
