import requests
import re
import time
from urllib.parse import urljoin
from pathlib import Path
import gzip
import shutil
import os

# Base URLs for downloading genomes/ (change this url for other organism groups)
BASE_URLS = {
    "protozoa": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/"
}

# Function to log messages
def log_message(message, log_file):
    with open(log_file, 'a') as file:
        file.write(message + "\n")

def find_assembly_directory(organism_name):
    max_retries = 3  # Number of retries
    retry_wait = 5  # Initial wait time in seconds, will increase exponentially

    for domain, base_url in BASE_URLS.items():
        for directory in ["reference", "representative", "latest_assembly_versions"]:
            dir_url = urljoin(base_url, f"{organism_name}/{directory}/")
            for attempt in range(max_retries):
                try:
                    response = requests.get(dir_url)
                    if response.ok:
                        match = re.findall(r'href="(GCF_[^"/]+)/"', response.text)
                        if match:
                            return urljoin(dir_url, match[-1] + '/'), domain
                except requests.exceptions.ConnectionError as e:
                    if attempt < max_retries - 1:
                        time.sleep(retry_wait)
                        retry_wait *= 2  # Double the wait time for the next retry
                    else:
                        log_message(f"Failed to connect to {dir_url} after {max_retries} attempts.", log_file)
                except Exception as e:
                    log_message(f"Error accessing {dir_url}: {e}", log_file)
    return None, None

def download_and_extract_file(url, organism_name, domain, log_file):
    try:
        file_name = f"{organism_name}_{domain}_{Path(url).name}"
        download_path = Path(file_name)
        extracted_file_path = download_path.with_suffix('')
        if extracted_file_path.exists():
            print(f"File {extracted_file_path} already exists. Skipping download for {organism_name} from {domain}.")
            return True
        response = requests.get(url, stream=True)
        if response.ok:
            with open(download_path, 'wb') as f:
                shutil.copyfileobj(response.raw, f)
            with gzip.open(download_path, 'rb') as f_in:
                with open(extracted_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(download_path)
            print(f"Downloaded and extracted genome to {extracted_file_path} for {organism_name} from {domain}.")
            return True
    except Exception as e:
        log_message(f"Error downloading or extracting for {organism_name}: {e}", log_file)
        return False
    return False

def download_genomes(organisms_file, log_file):
    with open(organisms_file, 'r') as file:
        organisms = file.read().splitlines()

    for organism_name in organisms:
        assembly_dir_url, domain = find_assembly_directory(organism_name)
        if assembly_dir_url:
            file_url = urljoin(assembly_dir_url, f"{Path(assembly_dir_url).name}_genomic.gbff.gz")
            if not download_and_extract_file(file_url, organism_name, domain, log_file):
                log_message(f"Failed to download genome for {organism_name}.", log_file)

if __name__ == "__main__":
    log_file_path = 'download_log.txt'
    organisms_file_path = 'list_organisms.txt'
    download_genomes(organisms_file_path, log_file_path)

