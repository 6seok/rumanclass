import multiprocessing as mp

# Function to divide the FASTA file into blocks
def chunk_fasta_file(fasta_file):
    with open(fasta_file, 'r') as f:
        chunk = []
        for line in f:
            if line.startswith(">") and chunk:
                yield chunk
                chunk = [line]
            else:
                chunk.append(line)
        if chunk:
            yield chunk

# Function to process blocks of the FASTA file
def parse_fasta_block(block, sequence_dict):
    hash = None
    for line in block:
        if line.startswith(">"):
            hash = line[1:].strip()  # Remove '>' and strip
        else:
            sequence_dict[hash] = line.strip()

# Function to process chunks of the TSV file
def process_tsv_chunk(chunk, sequence_dict):
    result = []
    for row in chunk:
        row = row.strip().split('\t')  # Split each row by tab
        hash = row[0]
        sequence = sequence_dict.get(hash, hash)  # Find the sequence or use the hash
        row[0] = sequence
        result.append('\t'.join(row))  # Save the modified row
    return result

# Function to divide the TSV file into chunks
def chunk_tsv_file(tsv_file, chunk_size=1000):
    with open(tsv_file, 'r') as f:
        chunk = []
        for line in f:
            if len(chunk) < chunk_size:
                chunk.append(line)
            else:
                yield chunk
                chunk = [line]
        if chunk:
            yield chunk

# Parse FASTA file in parallel using multiple cores
def parse_fasta_in_parallel(fasta_file, manager, num_workers):
    sequence_dict = manager.dict()  # Create a shared dictionary

    # Process FASTA file blocks in parallel
    with mp.Pool(num_workers) as pool:
        chunked_blocks = chunk_fasta_file(fasta_file)
        pool.starmap(parse_fasta_block, [(block, sequence_dict) for block in chunked_blocks])

    return sequence_dict

# Process TSV file in parallel using multiple cores
def process_tsv_in_parallel(tsv_file, sequence_dict, num_workers, output_file):
    with mp.Pool(num_workers) as pool:
        chunked_tsv = chunk_tsv_file(tsv_file)
        results = pool.starmap(process_tsv_chunk, [(chunk, sequence_dict) for chunk in chunked_tsv])

    with open(output_file, 'w') as out_file:
        for result_chunk in results:
            out_file.write('\n'.join(result_chunk) + '\n')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process FASTA and TSV files in parallel.")
    parser.add_argument("--input-seq", type=str, required=True, help="Path to the input FASTA file.")
    parser.add_argument("--input-table", type=str, required=True, help="Path to the input TSV file.")
    parser.add_argument("--output", type=str, required=True, help="Path to the output TSV file.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use (default: 4).")

    args = parser.parse_args()

    fasta_file = args.input_seq
    tsv_file = args.input_table
    output_file = args.output
    num_workers = args.threads

    # Create a Manager object for shared memory
    manager = mp.Manager()

    # Parse FASTA file in parallel
    sequence_dict = parse_fasta_in_parallel(fasta_file, manager, num_workers)

    # Process TSV file in parallel
    process_tsv_in_parallel(tsv_file, sequence_dict, num_workers, output_file)
