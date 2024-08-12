import os
import json
import csv
import argparse
import pandas as pd  # Import pandas for Excel file support

def extract_json_content(log_content):
    json_start = log_content.find("[")
    json_end = log_content.rfind("]") + 1
    json_content = log_content[json_start:json_end]
    return json_content

def process_log_files(log_directory, output_file, output_format):
    data = []

    for filename in os.listdir(log_directory):
        if filename.endswith(".log"):
            file_path = os.path.join(log_directory, filename)
            with open(file_path, "r") as file:
                log_content = file.read()

            json_content = extract_json_content(log_content)
            log_data = json.loads(json_content)

            data.append([
                filename.replace(".log", ""),
                log_data[0]["reads_in"],
                log_data[0]["reads_out"],
                log_data[0]["reads_removed"],
                log_data[0]["reads_removed_proportion"]
            ])

    data.sort(key=lambda x: x[4], reverse=True)

    if output_format.lower() == "csv":
        with open(output_file, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Sample_ID", "Reads_In", "Reads_Out", "Reads_Removed", "Reads_Removed_Proportion"])
            writer.writerows(data)
    elif output_format.lower() == "excel":
        df = pd.DataFrame(data, columns=["Sample_ID", "Reads_In", "Reads_Out", "Reads_Removed", "Reads_Removed_Proportion"])
        df.to_excel(output_file, index=False)

    print(f"Extracted data has been written to {output_file}.")

def main():
    parser = argparse.ArgumentParser(description="Extract information from log files.")
    parser.add_argument("log_directory", help="Directory containing the log files.")
    parser.add_argument("output_file", help="Output file.")
    parser.add_argument("--format", choices=["csv", "excel"], default="csv", help="Output format: csv (default) or excel.")
    args = parser.parse_args()

    process_log_files(args.log_directory, args.output_file, args.format)

if __name__ == "__main__":
    main()