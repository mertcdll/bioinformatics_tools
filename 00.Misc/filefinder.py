import os
import argparse
import time

class FileFinder:
    def __init__(self, directory, extensions):
        self.directory = directory
        self.extensions = extensions

    def find_files(self):
        start_time = time.time() 

        found_files = []
        for root, dirs, files in os.walk(self.directory):
            for file in files:
                for extension in self.extensions:
                    if file.endswith(extension):
                        file_path = os.path.join(root, file)
                        found_files.append({'file_name': file, 'real_path': file_path})

        end_time = time.time() 
        elapsed_time = end_time - start_time
        print(f"Execution time: {elapsed_time:.2f} seconds")

        return found_files

def main():
    parser = argparse.ArgumentParser(description='Find files with specific extensions in a directory and its subdirectories.')
    parser.add_argument("-d", metavar='--directory', help='The directory to search for files.')
    parser.add_argument("-e", metavar='--extensions', nargs='+', help='Space-separated list of file extensions to search for.')

    args = parser.parse_args()

    finder = FileFinder(args.d, args.e)
    result = finder.find_files()

    for entry in result:
        print(entry)

if __name__ == '__main__':
    main()



