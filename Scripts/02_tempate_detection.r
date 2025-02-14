###########################################################################
# Cane Toad Template Detection
#
# This script demonstrates template detection using the binMatch function 
# from the monitoR package, more details in Hafner and Katz, 2018.
#
# Workflow Overview:
# 1. Set parameters and load required libraries.
# 2. Load the detection template and set the detection cutoff.
# 3. Specify the audio file to analyze.
# 4. Calculate bin match scores, detect peaks, and extract detections.
# 5. Save outputs (detections) for further validation.
#
# References:
# - Hafner and Katz (2018)
#
# Dependencies: monitoR
#
# Note: Update file paths, template names, and parameter values as needed.
###########################################################################

# Load required library
library(monitoR)

#############################
# 1. Set Parameters & Working Directory
#############################

# Set your working directory (update this path as needed)
setwd("your/working/directory")

# Define the template file name and detection parameters
templateName <- "cane_toad_template.bt"  # Binary template file (stored locally)
detectionCutoff <- 8                     # Example cutoff score

# Directory where the binary template is stored (update as needed)
template_dir <- "./calltemplates"

# Directory to save detection outputs
output_dir <- "./TemplateDetections"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#############################
# 2. Load Detection Template
#############################

# Read the binary template from the specified directory
binTemplate <- readBinTemplates(files = templateName, dir = template_dir, ext = "bt")

# Extract the template ID (filename without extension)
template_id <- gsub(".bt", "", templateName)

# Set the detection cutoff score for the template
binTemplate@templates[[template_id]]@score.cutoff <- detectionCutoff

#############################
# 3. Specify Audio File to Analyze
#############################

# Replace with the path to the audio file you wish to analyze
audio_to_analyse <- "path/to/audio_file.wav"

# Extract recording name (without extension) to use in output file names.
# (Assumes the filename starts with a date-time stamp in "YYYYMMDDHHMMSS" format.)
Recording <- tools::file_path_sans_ext(basename(audio_to_analyse))

#############################
# 4. Calculate Scores and Detect Peaks
#############################

# Calculate bin match scores using the detection template
binscores <- binMatch(survey = audio_to_analyse,
                      templates = binTemplate,
                      quiet = TRUE,          # Suppress status messages
                      time.source = "filename", # Use filename for time extraction
                      write.wav = FALSE)     # Do not write out wave files

# Identify peaks in the score data
peaks <- findPeaks(binscores)

# Extract detection information from the peaks
detections <- getDetections(peaks)

#############################
# 5. Save Outputs for Further Validation
#############################

if (nrow(detections) == 0) {
  # Create an empty table if no detections were found
  emptyTable <- data.frame(Selection = character(),
                           `Begin Time (s)` = character(),
                           `End Time (s)` = character(),
                           `Low Freq (Hz)` = character(),
                           `High Freq (Hz)` = character(),
                           score = character(),
                           template = character(),
                           recording = character(),
                           dateTime = character(),
                           check.names = FALSE)
  
  write.table(emptyTable,
              file = file.path(output_dir, paste0(Recording, "_detections.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
} else {
  # Process detection results if available
  detections$Selection <- seq_len(nrow(detections))
  
  # Use the template's duration to calculate begin and end times
  template_duration <- binTemplate@templates[[template_id]]@duration
  detections$`Begin Time (s)` <- detections$time - (template_duration / 2)
  detections$`End Time (s)` <- detections$time + (template_duration / 2)
  
  # Add fixed frequency range (adjust if needed)
  detections$`Low Freq (Hz)` <- 0
  detections$`High Freq (Hz)` <- 11025
  
  # Record the audio file name for reference
  detections$recording <- Recording
  
  # Calculate dateTime assuming the recording name is in "YYYYMMDDHHMMSS" format
  detections$dateTime <- as.POSIXct(Recording, format = "%Y%m%d%H%M%S") + detections$time
  
  # Save the detection table to a text file
  write.table(detections,
              file = file.path(output_dir, paste0(Recording, "_detections.txt")),
              quote = FALSE, sep = "\t", row.names = FALSE)
}
