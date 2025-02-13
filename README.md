# Cane Toad Acoustic Classifier: A Free Tool for Monitoring Invasive Cane Toads in Australia

---

## Overview

This repository provides a freely available and user-friendly acoustic classifier for detecting the invasive cane toad (*Rhinella marina*) across Australia. The classifier was developed using BirdNET, a deep-learning-based species recognition algorithm, and trained on a large-scale passive acoustic monitoring (PAM) dataset from the Australian Acoustic Observatory (A2O).

## Workflow

The classifier development follows a structured workflow:

<p align="center">
  <img src="https://raw.githubusercontent.com/Leptobrachium/Gpshing/main/Classifier%20Development%20Workflow.jpg" alt="Workflow Diagram" width="500" height="650">
</p>

### **1. Data Generation**
- Cane toad calls were extracted from A2O data.
- Training and testing datasets were generated using template detection methods.

### **2. Classifier Training**
- The training dataset was used to develop a custom classifier in BirdNET.
- Additional non-target vocalizations and background noise were included to enhance performance.

### **3. Preliminary Assessment**
- The trained classifier was tested on a subset of audio data from selected A2O sites.
- Performance metrics such as precision, recall, and accuracy were evaluated across multiple temporal scales.

### **4. Audio Analysis**
- The classifier was applied to the full dataset.
- Probabilistic thresholds were computed to assess detection confidence.

### **5. Classifier Performance Assessment**
- Validation was conducted using manual verification of randomly sampled detections.

### **6. Detections Optimization (Optional)**
- Aggregated time-series features (ATF) were used to refine detection outputs.
- A conditional inference tree (CIT) model was applied to improve precision and recall.

---

## **Data Availability**
The trained classifier, scripts, and dataset are available in this repository.

The full Australian Acoustic Observatory (A2O) dataset can be accessed at:  
🔗 [A2O website](https://data.acousticobservatory.org/)

- **CaneToad_BirdNET_Classifier_V072024.tflite**: The TensorFlow Lite model file for the cane toad classifier.
- **Species_List_for_the_Classifier_V072024.txt**: A species list recognized by the classifier, including the cane toad.

To use this classifier, you need to have BirdNET installed. Follow the installation and setup instructions provided in the [BirdNET-Analyzer GitHub repository](https://github.com/kahst/BirdNET-Analyzer).

The full training dataset used to develop this classifier is available on Zenodo:

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.13826911">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13826911.svg" alt="DOI">
  </a>
</p>

---

## **Contact**
For questions or suggestions, feel free to contact me at **franco.leung@my.jcu.edu.au**
