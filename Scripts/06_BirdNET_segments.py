#!/usr/bin/env python3
###########################################################################
# Detection Segmentation & Extraction
#
# Adopted and modified from BirdNET-Analyzer.
# For more details, please visit:
#   https://kahst.github.io/BirdNET-Analyzer/
#
#
# This script extracts audio segments corresponding to BirdNET detections.
# Detections were treated as putative observations and sampled across
# confidence score intervals. For example, we randomly sampled 270 cane toad 
# detections per site (validating 30 per 0.1 interval) and then manually validated
# a subset using Kaleidoscope Lite.
#
#
# References:
# - Kahl, S., Wood, C. M., Eibl, M., & Klinck, H. (2021). BirdNET: A deep
#   learning solution for avian diversity monitoring. Ecological Informatics, 61, 101236.
#
# Dependencies: argparse, os, multiprocessing, numpy, audio, config, utils
#
# Note: Customize the paths via command-line arguments.
#
# Example usage:
#   python segments.py --audio "path/to/your/audio/" --results "path/to/your/results/" \
#       --o "path/to/your/output/"
###########################################################################

import argparse
import os
from multiprocessing import Pool
import numpy as np

import audio
import config as cfg
import utils

# Set numpy random seed from configuration
np.random.seed(cfg.RANDOM_SEED)


def detectRType(line: str):
    """Detects the type of result file based on its first line.
    
    Args:
        line: The first line of the result file.
    
    Returns:
        A string indicating the result file type: "table", "r", "kaleidoscope",
        "csv", or "audacity".
    """
    if line.lower().startswith("selection"):
        return "table"
    elif line.lower().startswith("filepath"):
        return "r"
    elif line.lower().startswith("indir"):
        return "kaleidoscope"
    elif line.lower().startswith("start (s)"):
        return "csv"
    else:
        return "audacity"


def parseFolders(apath: str, rpath: str, allowed_result_filetypes: list[str] = ["txt", "csv"]) -> list[dict]:
    """Recursively reads audio files and BirdNET result files from provided folders.
    
    Args:
        apath: Path to search for audio files.
        rpath: Path to search for result files.
        allowed_result_filetypes: Allowed file extensions for result files.
    
    Returns:
        A list of dictionaries in the form {"audio": path_to_audio, "result": path_to_result}.
    """
    data = {}
    # Normalize path separators
    apath = apath.replace("/", os.sep).replace("\\", os.sep)
    rpath = rpath.replace("/", os.sep).replace("\\", os.sep)

    # Retrieve all audio files using allowed file types (as defined in cfg)
    for root, _, files in os.walk(apath):
        for f in files:
            if f.rsplit(".", 1)[-1].lower() in cfg.ALLOWED_FILETYPES:
                key = f.rsplit(".", 1)[0]
                data[key] = {"audio": os.path.join(root, f), "result": ""}

    # Retrieve all result files (only those containing ".BirdNET." in the name)
    for root, _, files in os.walk(rpath):
        for f in files:
            if f.rsplit(".", 1)[-1] in allowed_result_filetypes and ".BirdNET." in f:
                key = f.split(".BirdNET.", 1)[0]
                if key in data:
                    data[key]["result"] = os.path.join(root, f)

    # Convert to list and filter to include only entries with a result file.
    flist = [entry for entry in data.values() if entry["result"]]
    print(f"Found {len(flist)} audio files with valid result files.")
    return flist


def sample_segments_by_confidence(segments, samples_per_bin=30):
    """
    Samples up to samples_per_bin segments for each 0.1 confidence interval between 0.1 and 1.0.
    
    Args:
        segments: List of detection segment dictionaries.
        samples_per_bin: Number of segments to sample per confidence bin.
    
    Returns:
        A list of sampled segments.
    """
    sampled = []
    # Loop over confidence bins: [0.1, 0.2), [0.2, 0.3), ..., [0.9, 1.0)
    for i in range(1, 10):
        lower = i * 0.1
        upper = lower + 0.1
        bin_segs = [seg for seg in segments if lower <= seg["confidence"] < upper]
        if len(bin_segs) > samples_per_bin:
            np.random.shuffle(bin_segs)
            sampled.extend(bin_segs[:samples_per_bin])
        else:
            sampled.extend(bin_segs)
    return sampled


def parseFiles(flist: list[dict], samples_per_bin=30):
    """
    Extracts and groups segments for all files, then samples detections per site.
    
    This function first extracts all segments from the audio/result file pairs.
    It then groups the segments by site (derived as the parent folder name of the audio file),
    and for each site and for each confidence interval (0.1 to 1.0) it samples up to samples_per_bin
    segments. Finally, the sampled segments are regrouped by audio file.
    
    Args:
        flist: List of dictionaries {"audio": path_to_audio, "result": path_to_result}.
        samples_per_bin: Number of segments to sample per confidence bin per site.
    
    Returns:
        A list of tuples: ((audio file, [segments]), ...).
    """
    # Group segments by site (assume site = parent folder of the audio file)
    site_segments: dict[str, list] = {}
    for f in flist:
        afile = f["audio"]
        rfile = f["result"]
        segments = findSegments(afile, rfile)
        # Determine site from parent directory name of the audio file
        site = os.path.basename(os.path.dirname(afile))
        for seg in segments:
            seg["site"] = site
            site_segments.setdefault(site, []).append(seg)
    
    # For each site, sample segments for each 0.1 confidence bin
    sampled_by_site = {}
    for site, segs in site_segments.items():
        sampled_by_site[site] = sample_segments_by_confidence(segs, samples_per_bin)

    # Regroup sampled segments by audio file
    segments_by_audio: dict[str, list] = {}
    seg_cnt = 0
    for site, segs in sampled_by_site.items():
        for seg in segs:
            segments_by_audio.setdefault(seg["audio"], []).append(seg)
            seg_cnt += 1

    print(f"After confidence bin sampling, found {seg_cnt} segments in {len(segments_by_audio)} audio files.")
    # Convert to list of tuples: (audio file, segments)
    flist = [tuple(item) for item in segments_by_audio.items()]
    return flist


def findSegments(afile: str, rfile: str):
    """
    Extracts detection segments from an audio file using its result file.
    
    Args:
        afile: Path to the audio file.
        rfile: Path to the BirdNET result file.
    
    Returns:
        A list of dictionaries, each containing:
        {"audio": afile, "start": start, "end": end, "species": species, "confidence": confidence}.
    """
    segments: list[dict] = []
    lines = utils.readLines(rfile)
    rtype = detectRType(lines[0])

    for i, line in enumerate(lines):
        if rtype == "table" and i > 0:
            d = line.split("\t")
            start = float(d[4])
            end = float(d[5])
            species = d[-2]
            confidence = float(d[-1])
        elif rtype == "audacity":
            d = line.split("\t")
            start = float(d[0])
            end = float(d[1])
            species = d[2].split(", ")[1]
            confidence = float(d[-1])
        elif rtype == "r" and i > 0:
            d = line.split(",")
            start = float(d[1])
            end = float(d[2])
            species = d[4]
            confidence = float(d[5])
        elif rtype == "kaleidoscope" and i > 0:
            d = line.split(",")
            start = float(d[3])
            end = float(d[4]) + start
            species = d[5]
            confidence = float(d[7])
        elif rtype == "csv" and i > 0:
            d = line.split(",")
            start = float(d[0])
            end = float(d[1])
            species = d[3]
            confidence = float(d[4])
        else:
            continue

        if confidence >= cfg.MIN_CONFIDENCE:
            segments.append({
                "audio": afile,
                "start": start,
                "end": end,
                "species": species,
                "confidence": confidence
            })
    return segments


def extractSegments(item: tuple[tuple[str, list], float, dict]):
    """
    Extracts and saves each audio segment as a separate file.
    
    Args:
        item: A tuple containing ((audio file path, segments), segment length, config).
    """
    afile = item[0][0]
    segments = item[0][1]
    seg_length = item[1]
    cfg.setConfig(item[2])

    print(f"Extracting segments from {afile}")

    try:
        sig, _ = audio.openAudioFile(afile, cfg.SAMPLE_RATE)
    except Exception as ex:
        print(f"Error: Cannot open audio file {afile}", flush=True)
        utils.writeErrorLog(ex)
        return

    for seg_cnt, seg in enumerate(segments, 1):
        try:
            start = int(seg["start"] * cfg.SAMPLE_RATE)
            end = int(seg["end"] * cfg.SAMPLE_RATE)
            offset = ((seg_length * cfg.SAMPLE_RATE) - (end - start)) // 2
            start = max(0, start - offset)
            end = min(len(sig), end + offset)

            if end > start:
                seg_sig = sig[start:end]
                outpath = os.path.join(cfg.OUTPUT_PATH, seg["species"])
                os.makedirs(outpath, exist_ok=True)

                seg_name = "{:.3f}_{}_{}.wav".format(
                    seg["confidence"],
                    seg_cnt,
                    os.path.splitext(os.path.basename(seg["audio"]))[0]
                )
                seg_path = os.path.join(outpath, seg_name)
                audio.saveSignal(seg_sig, seg_path)
        except Exception as ex:
            print(f"Error: Cannot extract segments from {afile}.", flush=True)
            utils.writeErrorLog(ex)
            return False
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract segments from audio files based on BirdNET detections."
    )
    parser.add_argument("--audio", default="path/to/your/audio/",
                        help="Path to folder containing audio files.")
    parser.add_argument("--results", default="path/to/your/results/",
                        help="Path to folder containing BirdNET result files.")
    parser.add_argument("--o", default="path/to/your/output/",
                        help="Output folder path for extracted segments.")
    parser.add_argument("--min_conf", type=float, default=0.1,
                        help="Minimum confidence threshold (between 0.01 and 0.99).")
    parser.add_argument("--samples_per_bin", type=int, default=30,
                        help="Number of segments to sample per 0.1 confidence bin per site.")
    parser.add_argument("--seg_length", type=float, default=3.0,
                        help="Length (in seconds) of the extracted segments.")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of CPU threads to use for parallel processing.")

    args = parser.parse_args()

    # Parse audio and result folders
    cfg.FILE_LIST = parseFolders(args.audio, args.results)
    cfg.OUTPUT_PATH = args.o
    cfg.CPU_THREADS = int(args.threads)
    cfg.MIN_CONFIDENCE = max(0.01, min(0.99, float(args.min_conf)))

    # Parse file list and prepare list of segments by grouping by site and sampling by confidence bins.
    cfg.FILE_LIST = parseFiles(cfg.FILE_LIST, samples_per_bin=int(args.samples_per_bin))
    flist = [
        (entry, max(cfg.SIG_LENGTH, float(args.seg_length)), cfg.getConfig())
        for entry in cfg.FILE_LIST
    ]

    # Extract segments using parallel processing if specified
    if cfg.CPU_THREADS < 2:
        for entry in flist:
            extractSegments(entry)
    else:
        with Pool(cfg.CPU_THREADS) as p:
            p.map(extractSegments, flist)