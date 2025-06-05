import os
import pandas as pd
import sys

def process_chromosome(chromosome, window_input_dir, npy_dir, variant_info_dir, output_dir):
    """Generate final window file for a given chromosome."""
    # Construct the path to the input window file
    window_input_file = os.path.join(window_input_dir, f"windows_chr{chromosome}_corrected_flip.tsv")
    final_output = os.path.join(output_dir, f"window_file_chr{chromosome}_corrected_flip.tsv")

    if not os.path.exists(window_input_file):
        print(f"Error: {window_input_file} not found. Skipping chromosome {chromosome}.")
        return

    # Process the window file
    df = pd.read_csv(window_input_file, sep="\t")

    df["LD_matrix"] = df["window_name"].apply(lambda x: os.path.join(npy_dir, f"{x}_ld_corrected_flip.npy"))
    df["variant_info_file"] = df["window_name"].apply(lambda x: os.path.join(variant_info_dir, f"{x}_variant_info_corrected_flip.txt"))

    # Save the final output
    df.to_csv(final_output, sep="\t", index=False)
    print(f"Processed chromosome {chromosome}: {final_output} created.")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python generate_final_windows_file.py <chromosome> <window_input_dir> <npy_dir> <variant_info_dir> <output_dir>")
        sys.exit(1)

    process_chromosome(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
