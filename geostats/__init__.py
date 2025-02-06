# -*- coding: utf-8 -*-
"""The geostats package is a python module that provides an interface to compute spatial statistics based on a shapefile for various datasets."""
# Import metadata
from .metadata import *
from .main import *
import os
import requests
import hashlib
import concurrent.futures
from tqdm import tqdm
import py7zr
import shutil
import time
import subprocess
import re

# ✅ Define MD5 hashes and expected sizes for each file (matches Zenodo)
FILE_DETAILS = {
    "geostats-data.7z.001": {
        "url": "https://zenodo.org/record/14714910/files/geostats-data.7z.001",
        "size": 5368709120,  # 5.4 GB in bytes
        "md5": "4bb0f51d22b2bfd576a731972e21ea7e"
    },
    "geostats-data.7z.002": {
        "url": "https://zenodo.org/record/14714910/files/geostats-data.7z.002",
        "size": 5368709120,
        "md5": "2013a2e6191cf57d7da9a0872c3e3fbc"
    },
    "geostats-data.7z.003": {
        "url": "https://zenodo.org/record/14714910/files/geostats-data.7z.003",
        "size": 5368709120,
        "md5": "27494df491cf20b7f27df60bb91119e3"
    },
    "geostats-data.7z.004": {
        "url": "https://zenodo.org/record/14714910/files/geostats-data.7z.004",
        "size": 5108605577,  # 5.1 GB in bytes
        "md5": "642fd3bb0aa43dc63914527137f813a8"
    },
}

def calculate_md5(filename):
    """Calculate the MD5 checksum of a file (matches Zenodo's hash)."""
    md5_hash = hashlib.md5()
    with open(filename, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()

def download_file_with_verification(url, filename, expected_size, expected_hash, max_retries=3):
    """Download a file, verify its MD5, and retry if incorrect."""
    for attempt in range(max_retries):
        print(f"📥 Downloading {filename} (Attempt {attempt+1}/{max_retries})")
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))

        # ✅ If file exists, verify size & hash
        if os.path.exists(filename):
            actual_size = os.path.getsize(filename)
            if actual_size == expected_size:
                calculated_hash = calculate_md5(filename)
                if calculated_hash == expected_hash:
                    print(f"✅ {filename} is already correctly downloaded. Skipping.")
                    return
                else:
                    print(f"❌ MD5 mismatch for {filename}. Redownloading.")
                    os.remove(filename)

        # ✅ Download with progress bar
        block_size = 1024
        with tqdm(total=total_size, unit='iB', unit_scale=True, desc=f"Downloading {filename}", leave=False) as t:
            with open(filename, 'wb') as f:
                for data in response.iter_content(block_size):
                    t.update(len(data))
                    f.write(data)

        actual_size = os.path.getsize(filename)
        if actual_size != expected_size:
            print(f"⚠️ Size mismatch for {filename}. Retrying...")
            os.remove(filename)
            time.sleep(2 ** attempt)  # Exponential backoff
            continue

        calculated_hash = calculate_md5(filename)
        if calculated_hash != expected_hash:
            print(f"❌ MD5 mismatch for {filename}. Retrying...")
            os.remove(filename)
            time.sleep(2 ** attempt)  # Exponential backoff
            continue

        print(f"✅ Successfully downloaded and verified {filename}.")
        return

    print(f"❌ Failed to download {filename} after {max_retries} attempts.")

def download_7z_parts():
    """Download and verify all 7z parts."""
    temp_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '7z_data')
    os.makedirs(temp_dir, exist_ok=True)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(
                download_file_with_verification,
                file_info["url"],
                os.path.join(temp_dir, filename),
                file_info["size"],
                file_info["md5"]
            ): filename
            for filename, file_info in FILE_DETAILS.items()
        }
        for future in concurrent.futures.as_completed(futures):
            filename = futures[future]
            try:
                future.result()
            except Exception as e:
                print(f"❌ Error downloading {filename}: {e}")

    return temp_dir

def extract_7z_parts(temp_dir):
    """Extract 7z parts using system-installed 7z with a real-time progress bar."""
    home_dir = os.path.expanduser("~")
    first_part = os.path.join(temp_dir, "geostats-data.7z.001")

    if not os.path.exists(first_part):
        print(f"❌ {first_part} not found. Extraction failed.")
        return

    print(f"🛠 Extracting {first_part} to {home_dir} using 7z CLI with multi-threading and progress tracking...")

    try:
        # Run 7z with multi-threading (-mmt=on) and real-time progress tracking
        command = ["7z", "x", first_part, f"-o{home_dir}", "-mmt=on", "-bsp1", "-bso0", "-y"]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

        # Initialize progress bar
        progress_bar = tqdm(total=100, unit="%", dynamic_ncols=True, leave=True)

        for line in iter(process.stdout.readline, ''):
            sys.stdout.write(line)  # Force output to be displayed in real-time
            sys.stdout.flush()  # Flush buffer immediately

            match = re.search(r'(\d+)%', line)
            if match:
                progress = int(match.group(1))
                progress_bar.n = progress
                progress_bar.refresh()

        process.wait()
        progress_bar.close()

        if process.returncode == 0:
            print(f"✅ Extraction completed successfully to {home_dir}")
            cleanup_temp_files(temp_dir)
        else:
            print(f"❌ Extraction failed with error code {process.returncode}")

    except Exception as e:
        print(f"❌ Error during extraction: {e}")

def cleanup_temp_files(temp_dir):
    """Delete downloaded 7z files and temporary directory."""
    try:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            print(f"🗑 Deleted temporary directory: {temp_dir}")
    except Exception as e:
        print(f"❌ Error while deleting temp files: {e}")

def has_subdirectories(directory):
    """Check if a directory contains subdirectories."""
    return any(os.path.isdir(os.path.join(directory, item)) for item in os.listdir(directory))

# ✅ Main execution
path = os.getenv('HOME') + '/geostats-data/'
if os.path.exists(path) and has_subdirectories(path):
    print('✅ Data already present with subdirectories, skipping download.')
else:
    print('📥 Downloading required datasets (>25GB)...')
    print('This may take a few minutes.')
    print('This will be done only once.')
    print('Please be patient.')
    temp_dir = download_7z_parts()

    print('📥 Extracting required datasets (>25GB)...')
    print('This may take a few minutes.')
    print('This will be done only once.')
    print('Please be patient.')
    extract_7z_parts(temp_dir)
