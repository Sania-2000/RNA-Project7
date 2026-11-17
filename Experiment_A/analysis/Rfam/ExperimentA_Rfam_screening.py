# Project 7 — Experiment A: Rfam screening of 45 RNA sequences
# Reconstructed from the Colab code used for the completed screening.
# Requires Google Colab and an internet connection to install Infernal / download Rfam.

import hashlib
import shutil
import subprocess
import zipfile
from pathlib import Path
from google.colab import files

# Upload the existing 45-RNA FASTA directly, or a ZIP containing it.
print("Select Task5_45RNA_Rfam_batch.fasta (or a ZIP containing that FASTA)")
uploaded = files.upload()
filename = next(iter(uploaded))
folder = Path("rfam_inputs")
folder.mkdir(exist_ok=True)
if filename.lower().endswith(".zip"):
    with zipfile.ZipFile(filename) as archive:
        archive.extractall(folder)
else:
    shutil.copy(filename, folder / Path(filename).name)

fasta_files = sorted(
    p for p in folder.rglob("*")
    if p.is_file() and p.suffix.lower() in {".fasta", ".fa", ".fna"}
)
if not fasta_files:
    raise ValueError("No FASTA input was found")

combined = Path("all_45_rnas.fasta")
combined.write_text(
    "".join(p.read_text().strip() + "\n" for p in fasta_files),
    encoding="utf-8",
)
headers = [line for line in combined.read_text().splitlines() if line.startswith(">")]
if len(headers) != 45 or len(set(headers)) != 45:
    raise ValueError(f"Expected 45 unique FASTA records; found {len(headers)}")
print("Number of sequences:", len(headers))
print("FASTA SHA256:", hashlib.sha256(combined.read_bytes()).hexdigest())

# Install Infernal and download the current public Rfam covariance-model collection.
subprocess.run(["apt-get", "-qq", "update"], check=True, stdout=subprocess.DEVNULL)
subprocess.run(["apt-get", "-qq", "install", "-y", "infernal"], check=True,
               stdout=subprocess.DEVNULL)
rfam_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
subprocess.run(["wget", "-q", "--show-progress", "-O", "Rfam.cm.gz", rfam_url], check=True)
subprocess.run(["gzip", "-dkf", "Rfam.cm.gz"], check=True)
subprocess.run(["cmpress", "-F", "Rfam.cm"], check=True)
print("Rfam.cm SHA256:", hashlib.sha256(Path("Rfam.cm").read_bytes()).hexdigest())

# Thresholds: each Rfam model's curated gathering cutoff (GA).
# Save both the machine-readable hit table and the full scan log.
subprocess.run([
    "cmscan", "--cut_ga", "--cpu", "2", "--noali",
    "--tblout", "Rfam_45_results.tbl",
    "-o", "Rfam_45_full_output.txt",
    "Rfam.cm", str(combined)
], check=True)

print("Rfam screening completed.")
print("Hit lines:", sum(
    bool(line.strip()) and not line.startswith("#")
    for line in Path("Rfam_45_results.tbl").read_text().splitlines()
))
# Download both records to keep the query-level log and reportable hits.
files.download("Rfam_45_results.tbl")
files.download("Rfam_45_full_output.txt")
