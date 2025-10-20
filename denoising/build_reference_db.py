from Bio import SeqIO
from pathlib import Path
from collections import Counter
import gzip
import sys 

# === CONFIGURATION ===
if len(sys.argv) != 3:
    sys.exit(f"Usage: {sys.argv[0]} <fastq_files> <outfile>")

fastq_files = sys.argv[1]  
output_fasta = sys.argv[2]

# === FUNCTION TO READ FASTQ ===
def read_fastq_sequences(file_path):
    """Yield sequences from a FASTQ or gzipped FASTQ file."""
    if file_path.suffix == ".gz":
        handle = gzip.open(file_path, "rt")
    else:
        handle = open(file_path, "r")
    with handle as f:
        for record in SeqIO.parse(f, "fastq"):
            yield str(record.seq)

# === MAIN FUNCTION ===
def collapse_fastq_to_fasta(input_dir, output_fasta):
    input_dir = Path(input_dir)
    fastq_files = sorted(input_dir.glob("*.fq*"))
    print(f"Found {len(fastq_files)} FASTQ files")

    seq_counter = Counter()

    for fq_file in fastq_files:
        print(f"Processing: {fq_file.name}")
        for seq in read_fastq_sequences(fq_file):
            seq_counter[seq] += 1

    print(f"Total unique sequences: {len(seq_counter):,}")

    with open(output_fasta, "w") as out_fa:
        for i, (seq, count) in enumerate(seq_counter.items(), 1):
            out_fa.write(f">Uniq{i};size={count};\n{seq}\n")

    print(f"✅ Wrote {len(seq_counter):,} unique sequences to {output_fasta}")

# === RUN ===
if __name__ == "__main__":
    collapse_fastq_to_fasta(fastq_files, output_fasta)
