### Scanning all files with argparser using terminal
### It has a function that can check all files in a directory if they are truncated or not, and creates a text file about their statuses


import os
import gzip
import zlib
import argparse

def file_validation_checker(directory, extension):
    output_file = os.path.join(directory, 'file_statuses.txt')  # Output file name
    result = []  

    
    for target_file in os.listdir(directory): #List all directory and iterate inside
        if target_file.endswith(extension):
            f_path = os.path.join(directory, target_file)
            f_status = target_file + ": "

            try: # add the try except block here
                with gzip.open(f_path, 'rb') as file:
                    file.read()  
                f_status += "Valid File"
            except (gzip.BadGzipFile, zlib.error):
                f_status += "Invalid or Truncated File"
            except IOError:
                f_status += "Read Error"

            result.append(f_status)

    
    with open(output_file, 'w') as res_file: # write the results
        for f_status in result:
            res_file.write(f_status + "\n")

    print("File statuses are in:", output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='It checks compressed files in a specified directory')
    parser.add_argument('directory', help='Directory path which includes the files to be checked')
    parser.add_argument('extension', help='File extension - It can be ".gz",".bed.gz", ".fastq.gz" etc.')

    args = parser.parse_args()

    file_validation_checker(args.directory, args.extension)

