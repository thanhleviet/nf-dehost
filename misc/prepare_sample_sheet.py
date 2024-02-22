#!/usr/bin/env python
# -*- coding: utf-8 -*-
# require python 3.10 or later

__author__ = "Thanh Le Viet"
__email__ = "thanh@cloudbinfies.com"
__version__ = "0.0.1"

from collections import Counter
from itertools import chain
import re
from pathlib import Path
import argparse

def scan_gz_files(directory):
    """
    Scan all files ending with .gz in the specified directory.

    Parameters:
        directory (str): Path to the directory to scan.

    Returns:
        list: A list of paths to files ending with .gz.
    """
    gz_files = []
    directory_path = Path(directory)
    for file in directory_path.rglob('*.gz'):
        gz_files.append(file)
    return gz_files


def read_input(filename: str) -> list[str]:
    """
    Read the contents of a file and return each line as a list of strings.

    Args:
    filename (str): The name of the file to read.

    Returns:
    list[str]: List of strings, each representing a line from the file.
    """
    with open(filename, "r") as f:
        return f.read().splitlines()


def get_kmer_list(text: str, min_kmer: int) -> list:
    """
    Generate a list of k-mers from the input text.

    Args:
    text (str): The input text
    min_kmer (int): The minimum length of k-mer to include in the list

    Returns:
    list: A list of k-mers
    """
    kmer_list = []
    for i in range(len(text)):
        kmer = text[i:]
        if len(kmer) >= min_kmer:
            kmer_list.append(kmer)
    return kmer_list


def extract_longest_common_substring(file_list: list[str], min_kmer: int) -> list[str]:
    """
    Extracts the k-mer list from each file in the file list and returns a flattened list of all k-mers.

    Args:
    file_list (list[str]): List of file names
    min_kmer (int): Minimum length of k-mer

    Returns:
    list[str]: Flattened list of all k-mers extracted from the files
    """
    # Get the k-mer list from each file in the file list
    substring_list = [get_kmer_list(x, min_kmer) for x in file_list]

    # Flatten the list of k-mer lists
    flatten_substring_list = list(chain.from_iterable(substring_list))

    return flatten_substring_list


def remove_extension(file_name: str, pattern: str) -> str:
    return re.sub(pattern, "", file_name)

if __name__ == "__main__":
    # fastq_files = read_input("fastq_files.csv")
    # # print(get_kmer_list(fastq_files[0], 50))
    # number_of_fastq_files = len(fastq_files)
    # rs = extract_longest_common_substring(fastq_files, 10)
    # common_substrings = Counter(rs).most_common(number_of_fastq_files)
    # sorted_common_substring = sorted(common_substrings, key=lambda x: len(x[0]), reverse=True)
    # patterns = "|".join(k[0] for k in sorted_common_substring[:2])
    # name = [re.sub(patterns, "", _name).split("/")[1] for _name in fastq_files]
    # print(name)
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="Scan directory for .gz files")

    # Add positional argument for directory path
    parser.add_argument("directory", type=str, help="Path to the directory to scan")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to scan .gz files
    gz_files = scan_gz_files(args.directory)
    fastq_files = [str(fq) for fq in gz_files]
    # Print the list of .gz files found
    if gz_files:
        number_of_fastq_files = len(fastq_files)
        rs = extract_longest_common_substring(fastq_files, 10)
        common_substrings = Counter(rs).most_common(number_of_fastq_files)
        sorted_common_substring = sorted(common_substrings, key=lambda x: len(x[0]), reverse=True)
        patterns = "|".join(k[0] for k in sorted_common_substring[:2])
        if len(fastq_files) % 2 == 0:
            print("sample_id,R1,R2")
            for i in range(0, len(fastq_files), 2):
                file_name = [remove_extension(Path(fastq_files[i]).name, patterns), remove_extension(Path(fastq_files[i + 1]).name, patterns)]
                if file_name[0] == file_name[1]:
                    print(f"{file_name[0]},{fastq_files[i]},{fastq_files[i + 1]}")
    else:
        print("No .gz files found in the specified directory.")