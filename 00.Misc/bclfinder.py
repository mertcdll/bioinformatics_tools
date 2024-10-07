import os
import re
import argparse

class BCLFinder:
    def __init__(self, directory):
        self.directory = directory

    def find_bcl(self):
        folder_pattern = r'\d{6}_([A-Z\d]+)'
        pattern = re.compile(folder_pattern)
        bcl_folders = []

        for root, dirs, files in os.walk(self.directory):
            for dir_name in dirs:
                folder_path = os.path.join(root, dir_name)
                runinfo_path = os.path.join(folder_path, 'RunInfo.xml')

                if pattern.match(dir_name) and os.path.isfile(runinfo_path):
                    
                    run_name = dir_name.split()[0]    
                    split_run = run_name.split("_")
                    date = split_run[0]
                    device_id = split_run[1]
                    try:
                        if split_run[2] == "RUO":
                            run_code = "_".join([split_run[2], split_run[3]])
                            cflowcell_id = split_run[4]

                            if "-" in cflowcell_id:
                                flowcell_id = "-".join([cflowcell_id.split("-")[0], cflowcell_id.split("-")[1]])
                            else:
                                flowcell_id = cflowcell_id

                        else:
                            run_code = split_run[2]
                            cflowcell_id = split_run[3]
                            
                            if "-" in cflowcell_id:
                                flowcell_id = "-".join([cflowcell_id.split("-")[0], cflowcell_id.split("-")[1]])
                            else:
                                flowcell_id = cflowcell_id
                    except IndexError:
                        print(f"Indexerror in folder: {folder_path}")
                        continue

                    bcl_dict = {
                        "date" : date,
                        "device_id": device_id,
                        "experiment_id": run_code,
                        "flow_cell_id": flowcell_id,
                        "run_name": run_name,
                        "real_path": folder_path
                    }
                    
                    bcl_folders.append(bcl_dict)
        
        return bcl_folders


def main():
    parser = argparse.ArgumentParser(description='Find BCL run folders in a directory and its subdirectories.')
    parser.add_argument("-d", metavar='--directory', help='The directory to search for files.')

    args = parser.parse_args()

    finder = BCLFinder(args.d)
    result = finder.find_bcl()

    for entry in result:
        print(entry)

if __name__ == '__main__':
    main()


