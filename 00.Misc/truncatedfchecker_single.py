import gzip
import zlib
import argparse

def check_gzip_validity(file_path):
    try:
        with gzip.open(file_path, 'rb') as file:
            file.read()
        print("File read successfully.")
    except gzip.BadGzipFile:
        print("The file is truncated or invalid.")
    except zlib.error:
        print("The file is truncated or invalid.")
    except IOError:
        print('Error reading file')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='It checks a single compressed file in a specified directory')
    parser.add_argument('file_path', help='Location of the file')

    args = parser.parse_args()

    check_gzip_validity(args.file_path)