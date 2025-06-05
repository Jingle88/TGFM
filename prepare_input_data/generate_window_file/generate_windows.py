import pandas as pd
import os
import sys

def generate_windows(bim_file, output_file, chrom, window_size=3000000, step_size=1000000):
    """Generate overlapping genomic windows for a specific chromosome from a PLINK bim file."""
    if not os.path.exists(bim_file):
        print(f"Error: {bim_file} not found.")
        return
    
    # Read bim file
    df = pd.read_csv(bim_file, sep='\s+', header=None, usecols=[0, 3], names=["chr", "pos"])
    
    # Filter for the specific chromosome
    df = df[df["chr"] == chrom]
    if df.empty:
        print(f"Warning: No data found for chromosome {chrom}.")
        return
    
    start_pos = max(1, df["pos"].min() - 1000000)  # Start 1MB before first variant
    last_variant_pos = df["pos"].max()  # Last variant position
    
    windows = []
    pos = start_pos
    while pos <= last_variant_pos:  # Keep generating windows until the last variant is covered
        window_start = pos
        window_end = pos + window_size
        window_name = f"{chrom}:{window_start}-{window_end}"
        
        windows.append([window_name, chrom, window_start, window_end])
        
        # If the last variant is within the current window, break
        if last_variant_pos <= window_end:
            break
        
        pos += step_size
    
    # Save the windows to a file
    window_df = pd.DataFrame(windows, columns=["window_name", "chr", "window_start", "window_end"])
    window_df.to_csv(output_file, sep="\t", index=False)
    print(f"Generated: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_windows.py <bim_file> <output_file> <chromosome>")
        sys.exit(1)
    
    bim_file = sys.argv[1]
    output_file = sys.argv[2]
    chrom = int(sys.argv[3])  # Convert chromosome to integer
    generate_windows(bim_file, output_file, chrom)
