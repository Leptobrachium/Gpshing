## Cane Toad Acoustic Classifier: A Free Tool for Monitoring Invasive Cane Toads in Australia

**Overview**

This repository provides a freely available and user-friendly acoustic classifier for detecting the invasive cane toad (*Rhinella marina*) across Australia. The classifier was developed using BirdNET, a deep-learning-based species recognition algorithm, and trained on a large-scale passive acoustic monitoring (PAM) dataset from the Australian Acoustic Observatory (A2O).

**Workflow**

The classifier development follows a structured workflow, as illustrated below:

1. Data Generation

- Cane toad calls were extracted from field recordings.

- Training and testing datasets were generated using template detection methods.

Classifier Training

The training dataset was used to develop a custom classifier in BirdNET.
Additional non-target vocalizations and background noise were included to enhance performance.

Preliminary Assessment

The trained classifier was tested on a subset of audio data from selected A2O sites.
Performance metrics such as precision, recall, and accuracy were evaluated.

Audio Analysis

The classifier was applied to the full dataset across 40 A2O sites (~89 years of audio).
Probabilistic thresholds were computed to assess detection confidence.

Classifier Performance Assessment

Validation was conducted using manual verification of randomly sampled detections.
Performance varied across sites, with accuracy exceeding 90% at many locations.

Detections Optimization (Optional)

Aggregated time-series features (ATF) were used to refine detection outputs.
A conditional inference tree (CIT) model was applied to improve precision and recall.
Optimized thresholds significantly reduced false positives compared to universal confidence thresholds.

## Data Availability
The trained classifier, scripts, and dataset are available in this repository.

The full Australian Acoustic Observatory (A2O) dataset can be accessed at [A2O website](https://data.acousticobservatory.org/).

**CaneToad_BirdNET_Classifier_V072024.tflite**: 
The TensorFlow Lite model file for the cane toad classifier.

**Species_List for_the_Classifier_V072024.txt**:
A species list recognized by the classifier, including the cane toad.

To use this classifier, you need to have BirdNET installed. Follow the installation and setup instructions provided in the [BirdNET-Analyzer GitHub repository](https://github.com/kahst/BirdNET-Analyzer).

The full training dataset used to develop this classifier is available on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13826911.svg)](https://doi.org/10.5281/zenodo.13826911)

## Contact
For questions or suggestions, feel free to contact me at franco.leung@my.jcu.edu.au



