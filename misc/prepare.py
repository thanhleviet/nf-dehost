#!/usr/bin/env python3
import os
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Scan a folder for paired-end FASTQ files and generate a sample sheet in CSV format.")
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
                filepath = os.path.join(root, filename)
                # Extract the sample ID from the folder name
                sample_id = os.path.basename(root)
                # Add the file to the dictionary
                if sample_id in files:
                    files[sample_id].append(filepath)
                else:
                    files[sample_id] = [filepath]

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
            # print(sample_id, filepaths)
            # Check if there are two files with the same sample ID
            if len(filepaths) == 2:
                # Determine which file is R1 and which is R2
                r1_path = filepaths[0] if filepaths[0].endswith("_1.fq.gz") else filepaths[1]
                r2_path = filepaths[1] if filepaths[1].endswith("_2.fq.gz") else filepaths[0]
                # Write the row to the CSV file
                print(sample_id, r1_path, r2_path)
                writer.writerow([sample_id, r1_path, r2_path])
            elif len(filepaths) == 4:
                for i in range(1,3):
                    # Determine which file is R1 and which is R2
                    _index = i if i == 2 else 0
                    r1_path = filepaths[0+_index] if filepaths[0].endswith("_1.fq.gz") else filepaths[1]
                    r2_path = filepaths[1+_index] if filepaths[1].endswith("_2.fq.gz") else filepaths[0]
                    # print(sample_id,r1_path,r2_path)
                    writer.writerow([sample_id, r1_path, r2_path])
if __name__ == "__main__":
    main()