import os
import polars as pl
import argparse

def main():
    parser = argparse.ArgumentParser(description='Concatenate TSV files into a parquet file')
    parser.add_argument('--dir', required=True, help='Directory containing TSV files')
    parser.add_argument('--output', default='output.parquet', help='Output parquet file name')
    
    args = parser.parse_args()
    
    # Get all TSV files in directory
    tsv_files = [
        os.path.join(args.dir, f) 
        for f in os.listdir(args.dir) 
        if f.endswith('.tsv')
    ]
    
    if not tsv_files:
        print(f"No TSV files found in {args.dir}")
        return
    
    # Read and concatenate
    df = pl.concat([
        pl.read_csv(file, separator="\t")
        for file in tsv_files
    ])
    
    # Save as parquet
    df.write_parquet(args.output)
    print(f"Concatenated {len(tsv_files)} files into {args.output}")

if __name__ == "__main__":
    main()