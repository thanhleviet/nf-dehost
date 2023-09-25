import argparse
import csv

def main(input_file):
    with open(input_file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        rows = [row for row in reader]

    cleaned_rows = [[row[0].replace('_1.non_host', '').replace('_EKDL230014950-1A_HFT3VDSX7_L2_1.fq', ''), row[1]] for row in rows]

    merged_rows = []
    i = 0
    while i < len(cleaned_rows):
        removed_human = 0
        if i+1 < len(cleaned_rows):
            cleaned_row1 = cleaned_rows[i]
            cleaned_row2 = cleaned_rows[i+1]
            cleaned_row1[0] = cleaned_row1[0].replace('_1.non_host', '').replace('_EKDL230014950-1A_HFT3VDSX7_L2_1.fq', '')
            cleaned_row2[0] = cleaned_row2[0].replace('_1.non_host', '').replace('_EKDL230014950-1A_HFT3VDSX7_L2_1.fq', '')
            if cleaned_row1[0] == cleaned_row2[0]:
                merged_rows.append([cleaned_row1[0], cleaned_row1[1], cleaned_row2[1], int(cleaned_row2[1]) - int(cleaned_row1[1])])
                i += 2
            else:
                merged_rows.append([cleaned_row1,0,0,removed_human])
                i += 1
        else:
            merged_rows.append([cleaned_rows[i],0,0,removed_human])
            i += 1

    with open('output.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sample_Name', 'Project_ID', 'Non_Host_Reads', 'Original_Reads', 'Removed_human'])
        for row in merged_rows:
            writer.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Group every two consecutive lines in a CSV file.')
    parser.add_argument('input_file', type=str, help='path to the input CSV file')
    args = parser.parse_args()

    main(args.input_file)