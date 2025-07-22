''' The following scipt takes in a directory path as input, and 
    unzips all the gzip and zip files in that directory.
    It will not overwrite existing files.
    Usage: python unzip_if_needed.py <directory_path>'''

import os
import sys
import gzip
import zipfile
import pandas as pd

sys.path.insert(1, '../../')
sys.path.insert(1, '../')
sys.path.insert(1, '../../../')

from src.hp import unzip_file
from src.hp import ungzip_file

def main(path):
    if not os.path.isdir(path):
        print(f"Invalid path: {path}")
        return

    for root, _, files in os.walk(path):
        for filename in files:
            full_path = os.path.join(root, filename)

            if filename.endswith('.gzip'):
                dest_path = full_path[:-5]  # remove .gzip
                if not os.path.exists(dest_path):
                    ungzip_file(full_path, dest_path)
                else:
                    print(f"Skipped existing file: {dest_path}")

            elif filename.endswith('.zip'):
                unzip_file(full_path, root)

            if filename.endswith('.gz'):
                dest_path = full_path[:-3]  # remove .gz
                if not os.path.exists(dest_path):
                    ungzip_file(full_path, dest_path)
                else:
                    print(f"Skipped existing file: {dest_path}")

    # finally formatting the data for R:
    W = pd.read_csv(f"{path}salmon_raw_counts_for_way_pipeline_whites.tsv", index_col=0, sep='\t')
    W.to_csv(f"{path}salmon_raw_counts_for_way_pipeline_whites.tsv", sep='\t')
    B = pd.read_csv(f"{path}salmon_raw_counts_for_way_pipeline.tsv", index_col=0, sep='\t')
    B.to_csv(f"{path}salmon_raw_counts_for_way_pipeline.tsv", sep='\t')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python unzip_if_needed.py <directory_path>")
    else:
        main(sys.argv[1])