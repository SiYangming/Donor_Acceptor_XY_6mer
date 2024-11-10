import os
import os.path
import csv
import argparse
### Reading FASTA sequences from a file, step 1

def read_FASTA_strings(filename):
    with open(filename) as file:
        return file.read().split('>')[1:]
    
### Reading FASTA sequences from a file, step 2

def read_FASTA_entries(filename):
    return [seq.partition('\n') for seq in read_FASTA_strings(filename)]
    
    
### Reading FASTA sequences from a file, step 3

def read_FASTA_sequences(filename):
    return [[seq[0].split(" ")[0], seq[2].replace('\n', '')]           # delete newlines
             for seq in read_FASTA_entries(filename)]

def write_csv_hearders(file_name, headers):
    with open(file_name, 'w', newline='', encoding="utf-8")as f:
        f_csv = csv.writer(f)
        f_csv.writerow(headers)
        #f_csv.writerows(headers)

def write_csv_rows(file_name, rows):
    with open(file_name, 'a+', newline='', encoding='utf-8')as f:
        f_csv = csv.writer(f)
        f_csv.writerow(list(map(str, rows)))
    f.close()

def feature_extract(dnas, seq_type = "dna", region=50,writeToFile = True, file = "features.csv", k = 6):
    featureMatrix = []
    if seq_type != "rna" and seq_type != "dna":
        print("Only support rna or dna~~")
        exit(0)
        
    for i in range(len(dnas)):
        if seq_type == "rna":
            # dna = [dnas[i][0], dnas[i][1][:50] + dnas[i][1][-50:]]
            dna = [dnas[i][0], dnas[i][1][:region] + dnas[i][1][-region:]]
        else:
            dna = dnas[i]
        if writeToFile:
            write_csv_rows(file, kmers_extract(dna, k))
        else:
            featureMatrix.append(kmers_extract(dna, k))
    
    if writeToFile:
        None
    else:
        return featureMatrix
    
def kmers_extract(dna, k = 6):
    features = []
    for i in range(len(dna[1])-int(k)+1):
        features.append(dna[1][i:i+int(k)])
    features.insert(0, dna[0])
    return features
  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess data.")
    parser.add_argument("--input_sequence", help="")
    parser.add_argument("--output_features", help="")
    parser.add_argument("--kmer", help="")
    parser.add_argument("--seq_type", default="dna",help="")
    parser.add_argument("--region",help="", default=50,type=int)
    args = parser.parse_args()
    
    # seqs = read_FASTA_sequences("features_seq/altAcc.fa")
    # # remember delete file bdfore running
    # filename = "features_seq/6mer-features.csv" 
    # # Test
    # #feature_extract(seqs[1:10],writeToFile=True,file=filename)
    # feature_extract(seqs, writeToFile=True, file=filename)
    
    input_file = args.input_sequence
    output_file = args.output_features
    k = args.kmer
    seq_type = args.seq_type
    region = args.region
    
    seqs = read_FASTA_sequences(input_file) 
    if os.path.exists(output_file):
        os.system("rm -f %s" % (output_file))
    
    if seq_type == "dna":
        feature_extract(seqs, writeToFile=True, file=output_file, k = k)
    elif seq_type == "rna":
        feature_extract(seqs, seq_type=seq_type, region=region, writeToFile=True, file=output_file, k = k)
    os.system('pigz -p 7 %s' % (output_file))
