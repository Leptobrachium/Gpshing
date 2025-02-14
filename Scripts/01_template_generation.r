###########################################################################
# Cane Toad Acoustic Template Generation 
#
# This script demonstrates the process of generating training data for 
# cane toad call detection using the monitoR package, 
# more details in Hafner and Katz, 2018.
#
# Workflow Overview:
# 1. View a spectrogram of a sample call.
# 2. Generate a binary-point acoustic template (3-second snippet) 
#    based on a typical cane toad call (see Fig. S2 in the manuscript).
# 3. Adjust amplitude cutoff thresholds.
# 4. Save outputs (detections, templates, plots) for further use.
#
# References:
# - Hafner and Katz (2018)
#
# Dependencies: monitoR, tuneR
#
# Note: Replace file paths and parameter values with those appropriate for your data.
###########################################################################

# Load required libraries
library(monitoR)
library(tuneR)

# Set your working directory (update this path as needed)
setwd("your/working/directory")

#############################
# 1. View Sample Spectrogram
#############################

# Replace with the path to your sample call .wav file
viewSpec("path/to/sample_call.wav")

#####################################
# 2. Generate a Binary Acoustic Template
#####################################

# Define a temporary file path for the template
template_fp <- file.path(tempdir(), "template.wav")

# Load the audio file containing a typical cane toad call
# Replace with the correct path to your template audio file
template <- readWave("path/to/template_audio.wav")

# Save the audio file to the temporary directory (required for template creation)
writeWave(template, template_fp)

# Create a binary template using a 3-second snippet with target frequency range (e.g., from 2 to 5 seconds, frequency = 0-2000Hz)
wbt <- makeBinTemplate(template_fp, amp.cutoff = "i", overlap = 50, 
                       frq.lim = c(0.0, 2.0), t.lim = c(2, 5), 
                       select = 'rectangle', name = "cane_toad_template")

####################################
# 3. Adjust Amplitude Cutoff Thresholds
####################################

# Adjust the amplitude cutoff thresholds for the template
templateCutoff(wbt) <- c(10, 5)

#########################################
# 4. Save Template for Further Use
#########################################

# Save the binary template for future use
writeBinTemplates(wbt, dir = ".", ext = "bt", parallel = FALSE)

# Save a plot of the binary template as a PNG file
png(filename = "template_plot.png")
plot(wbt)
dev.off()