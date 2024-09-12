#!/usr/bin/env python3
import os
import csv
import argparse
import logging

logging.basicConfig(level=logging.INFO)
def parse_args():
    parser = argparse.ArgumentParser(description="Scan a folder for paired-end FASTQ files and generate a sample sheet in CSV format (written for Novogene folder structure).")
    parser.add_argument("path", help="The path to the folder to scan.")
    parser.add_argument("output", help="The path to the output CSV file.")
    return parser.parse_args()

def main():
    args = parse_args()

    # Initialize the dictionary of files
    files = {}

    # Recursively scan the path for files
    for root, dirs, filenames in os.walk(args.path):
        for filename in filenames:
            if filename.endswith(".fq.gz"):
                logging.info(f"Found file: {filename}")
                filepath = os.path.join(root, filename)
                # Extract the sample ID from the folder name
                sample_id = os.path.basename(root)
                # Add the file to the dictionary
                if sample_id in files:
                    files[sample_id].append(filepath)
                else:
                    files[sample_id] = [filepath]
    logging.info(f"Found {len(files)} samples.")
    # Verify that all files exist
    for filepaths in files.values():
        for filepath in filepaths:
            if not os.path.exists(filepath):
                raise ValueError(f"File not found: {filepath}")
    # print(files)
    # Write the list of files to a new CSV file
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample_id", "R1", "R2"])
        for sample_id, filepaths in files.items():
            # Sort the file paths by name
            filepaths.sort()
            print(f"Sample ID: {sample_id}, Number of files: {len(filepaths)}")

            # Group files into pairs (R1 and R2)
            for i in range(0, len(filepaths), 2):
                if i + 1 < len(filepaths):
                # Determine which file is R1 and which is R2
                    if filepaths[i].endswith("_1.fq.gz") and filepaths[i+1].endswith("_2.fq.gz"):
                        r1_path, r2_path = filepaths[i], filepaths[i+1]
                    elif filepaths[i].endswith("_2.fq.gz") and filepaths[i+1].endswith("_1.fq.gz"):
                        r1_path, r2_path = filepaths[i+1], filepaths[i]
                    else:
                        print(f"Warning: Unexpected file naming for {sample_id}: {filepaths[i]}, {filepaths[i+1]}")
                        continue
                # Write the row to the CSV file
                    print(f"Writing: {sample_id}, {r1_path}, {r2_path}")
                    writer.writerow([sample_id, r1_path, r2_path])
                else:
                    # Handle the case of an odd number of files
                    print(f"Warning: Unpaired file for {sample_id}: {filepaths[i]}")
                    writer.writerow([sample_id, filepaths[i], ""])
if __name__ == "__main__":
    main()