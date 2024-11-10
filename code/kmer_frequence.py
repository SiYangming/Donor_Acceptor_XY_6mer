import pandas as pd
import os
import argparse

def Write_kmer_frequence(input_file, output_file, gzip=True):
    features = pd.read_csv(input_file, compression='gzip', header=None, index_col = 0)
    freqs = features.apply(pd.value_counts)
    freqs.to_csv(output_file, compression='gzip', na_rep = "NA")
    
def Write_kmer_frequence_BigFile(input_file, output_file, gzip=True):
    features = pd.read_csv(input_file, compression='gzip', header=None, index_col = 0, iterator=True)
    loop = True
    chunkSize = 1000000
    i = 1
    while loop:
        try:
            if i==1:
                chunk = features.get_chunk(chunkSize)
                freqs = chunk.apply(pd.value_counts)
            else:
                chunk = features.get_chunk(chunkSize)
                freqs = freqs + chunk.apply(pd.value_counts)
            i = i + 1
        except StopIteration:
            loop = False
            print("Iteration is stopped.")
    freqs.to_csv(output_file, compression='gzip', na_rep = "NA")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess data.")
    parser.add_argument("--input_file", help="")
    parser.add_argument("--output_prefix", help="")
    parser.add_argument("--file_size", help="", default="Small")
    args = parser.parse_args()
    
    # test = pd.read_csv("merAnalysis/test.csv.gz", compression='gzip', header=None, index_col = 0)
    # testT = test.T
    # testT[0,]
    # freqs = test.apply(pd.value_counts)
    # freqs.to_csv("merAnalysis/test_freqs.csv.gz", compression='gzip',)
    # test = pd.read_csv("merAnalysis/altAcc_6mer_features.csv.gz", compression='gzip', header=None, index_col = 0)
    # freqs = test.apply(pd.value_counts)
    # path_head, path_ext = os.path.split("merAnalysis/test.csv.gz")
    # file_head, file_ext = os.path.splitext("merAnalysis/test.csv.gz")
    # freqs.to_csv("merAnalysis/test_freqs.csv.gz", compression='gzip', na_rep = "NA")
    
    
    input_file = args.input_file
    file_size = args.file_size
    path, filename = os.path.split(input_file)
    output_file = args.output_prefix + filename.replace(".csv.gz","_freqs.csv.gz")
    if file_size == "Big":
        Write_kmer_frequence_BigFile(input_file, output_file)
    else:
        Write_kmer_frequence(input_file, output_file)
    

