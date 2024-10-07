import argparse
import pandas as pd

class ExcelMerger:
    def __init__(self, file_path):
        self.file_path = file_path
        self.sheet_name = "Total reads and bases"

    def merge_sheets(self):
        df = pd.read_excel(self.file_path, sheet_name=self.sheet_name)

        sum_df = df.groupby("sample_id").agg({'Total_reads_number': 'sum', 'Total_yield_base': 'sum'})
        sum_df.reset_index(inplace=True)

        with pd.ExcelWriter(self.file_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            sum_df.to_excel(writer, sheet_name=self.sheet_name, index=False, header=True)

def main():
    parser = argparse.ArgumentParser(description="Merge sheets in an Excel file.")
    parser.add_argument("--file_path", help="Path to the Excel file.")
    args = parser.parse_args()

    excel_merger = ExcelMerger(args.file_path)

    excel_merger.merge_sheets()

if __name__ == "__main__":
    main()
