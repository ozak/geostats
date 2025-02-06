# -*- coding: utf-8 -*-
"""The geostats package is a Python module that provides an interface to compute spatial statistics based on a shapefile for various datasets."""

import os
import requests
import concurrent.futures
from tqdm import tqdm
import py7zr
import shutil  # Used for deleting folders

# Securely load Zenodo token from environment variable
ZENODO_TOKEN = os.getenv("ZENODO_TOKEN")
if not ZENODO_TOKEN:
    raise ValueError("⚠️ ERROR: ZENODO_TOKEN environment variable is not set! Please set it before running.")

def download_file_with_progress(url, filename):
    """Download a file with a progress bar using Zenodo authentication."""
    headers = {"Authorization": f"Bearer {ZENODO_TOKEN}"}
    response = requests.get(url, stream=True, headers=headers)

    if response.status_code == 403:
        print(f"❌ Access denied: Make sure your Zenodo token has read permissions.")
        return
    elif response.status_code != 200:
        print(f"❌ Error {response.status_code} downloading {filename}: {response.text}")
        return

    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024  # 1 KiB

    # Use tqdm for the progress bar
    with tqdm(total=total_size, unit='iB', unit_scale=True, desc=f"Downloading {filename}", leave=False) as t:
        with open(filename, 'wb') as f:
            for data in response.iter_content(block_size):
                t.update(len(data))
                f.write(data)

def download_7z_parts():
    """Download 7z split parts in parallel from Zenodo (private access)."""
    temp_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), '7z_data')
    os.makedirs(temp_dir, exist_ok=True)

    # Define split 7z files with new naming convention
    base_url = "https://zenodo.org/api/records/14714910/files/"
    parts = [f"geostats-data.7z.{str(i).zfill(3)}" for i in range(1, 5)]  # Adjust range as needed

    # Construct full URLs
    files = {part: base_url + part for part in parts}

    # Parallel downloading
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(download_file_with_progress, url, os.path.join(temp_dir, filename)): filename
            for filename, url in files.items()
        }
        for future in concurrent.futures.as_completed(futures):
            filename = futures[future]
            try:
                future.result()
                print(f"✅ Downloaded {filename} successfully.")
            except Exception as e:
                print(f"❌ Error downloading {filename}: {e}")

    return temp_dir  # Return temp directory path for cleanup

def extract_7z_parts(temp_dir):
    """Extract downloaded 7z parts into the home directory using py7zr and clean up after extraction."""
    home_dir = os.path.expanduser("~")

    # Locate the first part (py7zr only needs .001 to extract the full set)
    first_part = os.path.join(temp_dir, "geostats-data.7z.001")
    if not os.path.exists(first_part):
        print(f"❌ Error: {first_part} not found. Ensure all parts are downloaded.")
        return

    print(f"🛠 Extracting {first_part} to {home_dir}...")

    # Extract using py7zr
    try:
        with py7zr.SevenZipFile(first_part, mode='r') as archive:
            archive.extractall(path=home_dir)
        print(f"✅ Extraction completed successfully to {home_dir}")
        
        # Clean up downloaded files and temp directory
        cleanup_temp_files(temp_dir)

    except Exception as e:
        print(f"❌ Error during extraction: {e}")

def cleanup_temp_files(temp_dir):
    """Delete the downloaded 7z files and temporary directory."""
    try:
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)  # Deletes entire directory and contents
            print(f"🗑 Deleted temporary directory: {temp_dir}")
    except Exception as e:
        print(f"❌ Error while deleting temp files: {e}")

# Check whether data is present or not
path = os.getenv('HOME') + '/geostats-data/'
if not os.path.exists(path):
    print('📥 Downloading required datasets (>25GB).')
    print('This may take a few minutes.')
    print('This will be done only once.')
    print('Please be patient.')
    
    temp_dir = download_7z_parts()  # Step 1: Download in parallel
    extract_7z_parts(temp_dir)      # Step 2: Extract using py7zr and clean up
else:
    print('Data already present, no need to download.')