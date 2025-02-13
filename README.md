<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cane Toad Acoustic Classifier</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 20px;
            padding: 20px;
            background-color: #f4f4f4;
        }
        h1, h2 {
            color: #333;
        }
        ul {
            padding-left: 20px;
        }
        p {
            margin-bottom: 10px;
        }
        .container {
            max-width: 800px;
            margin: auto;
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0px 0px 10px rgba(0,0,0,0.1);
        }
        .workflow-image {
            text-align: center;
        }
        .workflow-image img {
            width: 70%;
            max-width: 500px;
        }
        .doi-badge {
            display: flex;
            justify-content: center;
            margin-top: 20px;
        }
    </style>
</head>
<body>

<div class="container">
    <h1>Cane Toad Acoustic Classifier: A Free Tool for Monitoring Invasive Cane Toads in Australia</h1>

    <h2>Overview</h2>
    <p>
        This repository provides a freely available and user-friendly acoustic classifier for detecting the invasive cane toad 
        (<i>Rhinella marina</i>) across Australia. The classifier was developed using BirdNET, a deep-learning-based species 
        recognition algorithm, and trained on a large-scale passive acoustic monitoring (PAM) dataset from the Australian Acoustic 
        Observatory (A2O).
    </p>

    <h2>Workflow</h2>

    <div class="workflow-image">
        <img src="images/workflow_diagram.jpg" alt="Workflow Diagram">
    </div>

    <p>The classifier development follows a structured workflow, as illustrated below:</p>

    <h3>1. Data Generation</h3>
    <ul>
        <li>Cane toad calls were extracted from A2O data.</li>
        <li>Training and testing datasets were generated using template detection methods.</li>
    </ul>

    <h3>2. Classifier Training</h3>
    <ul>
        <li>The training dataset was used to develop a custom classifier in BirdNET.</li>
        <li>Additional non-target vocalizations and background noise were included to enhance performance.</li>
    </ul>

    <h3>3. Preliminary Assessment</h3>
    <ul>
        <li>The trained classifier was tested on a subset of audio data from selected A2O sites.</li>
        <li>Performance metrics such as precision, recall, and accuracy were evaluated across multiple temporal scales.</li>
    </ul>

    <h3>4. Audio Analysis</h3>
    <ul>
        <li>The classifier was applied to the full dataset.</li>
        <li>Probabilistic thresholds were computed to assess detection confidence.</li>
    </ul>

    <h3>5. Classifier Performance Assessment</h3>
    <ul>
        <li>Validation was conducted using manual validation of randomly sampled detections.</li>
    </ul>

    <h3>6. Detections Optimization (Optional)</h3>
    <ul>
        <li>Aggregated time-series features (ATF) were used to refine detection outputs.</li>
        <li>A conditional inference tree (CIT) model was applied to improve precision and recall.</li>
    </ul>

    <h2>Data Availability</h2>
    <p>The trained classifier, scripts, and dataset are available in this repository.</p>
    <p>The full Australian Acoustic Observatory (A2O) dataset can be accessed at <a href="https://data.acousticobservatory.org/">A2O website</a>.</p>

    <ul>
        <li><b>CaneToad_BirdNET_Classifier_V072024.tflite</b>: The TensorFlow Lite model file for the cane toad classifier.</li>
        <li><b>Species_List for_the_Classifier_V072024.txt</b>: A species list recognized by the classifier, including the cane toad.</li>
    </ul>

    <p>To use this classifier, you need to have BirdNET installed. Follow the installation and setup instructions provided in the 
        <a href="https://github.com/kahst/BirdNET-Analyzer">BirdNET-Analyzer GitHub repository</a>.</p>

    <p>The full training dataset used to develop this classifier is available on Zenodo:</p>

    <div class="doi-badge">
        <a href="https://doi.org/10.5281/zenodo.13826911">
            <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13826911.svg" alt="DOI">
        </a>
    </div>

    <h2>Contact</h2>
    <p>For questions or suggestions, feel free to contact me at <b>franco.leung@my.jcu.edu.au</b></p>
</div>

</body>
</html>
