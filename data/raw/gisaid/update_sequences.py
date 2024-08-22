import subprocess
import glob
import hashlib
import pandas as pd

checksum_file = "checksum.txt"

def concatenate_fasta_files():
    # Create the output filename
    fasta_out = "H5N1_concatenated.fasta"

    # Get all .fasta and .fa files in the current directory
    fasta_files = [file for file in glob.glob("*.fasta") + glob.glob("*.fa") if file != fasta_out]
    
    if not fasta_files:
        print("No FASTA files found in the current directory.")
        return
    
    # Sort the files to ensure consistent order
    fasta_files.sort()
    
    # Use cat command to concatenate all files
    cat_command = f"cat {' '.join(fasta_files)} > {fasta_out}"
    
    try:
        # Execute the cat command
        subprocess.run(cat_command, shell=True, check=True)
        print(f"Successfully concatenated {len(fasta_files)} FASTA files into {fasta_out}")

        # Calculate the MD5 checksum of the output file
        with open(fasta_out, "rb") as file:
            file_content = file.read()
            fasta_checksum = hashlib.md5(file_content).hexdigest()

        print(f"MD5 checksum of updated fasta file is:\n\t{fasta_checksum}")

        with open(checksum_file, "w") as file:
            file.write(f"FASTA: {fasta_checksum}\n")
        
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while concatenating files: {e}")

def concatenate_excel_files():
    # Create the output filename
    metadata_out = "H5N1_concatenated.csv"

    # Get all Excel files in the current directory
    excel_files = glob.glob("*.xlsx") + glob.glob("*.xls")
    
    if not excel_files:
        print("No Excel files found in the current directory.")
        return
    
    # Sort the files to ensure consistent order
    excel_files.sort()
    
    # Create an empty list to store DataFrames
    dfs = []
    
    # Read each Excel file and append to the list
    for file in excel_files:
        try:
            df = pd.read_excel(file)
            dfs.append(df)
            print(f"Read {file} successfully.")
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    if not dfs:
        print("No data was read from the Excel files.")
        return
    
    # Concatenate all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Save the combined DataFrame to CSV
    combined_df.to_csv(metadata_out, index=False)
    
    print(f"Successfully concatenated {len(excel_files)} Excel files into {metadata_out}")

    # Calculate the MD5 checksum of the output file
    with open(metadata_out, "rb") as file:
        file_content = file.read()
        metadata_checksum = hashlib.md5(file_content).hexdigest()

    print(f"MD5 checksum of updated metadata file is:\n\t{metadata_checksum}")

    with open(checksum_file, "a") as file:
        file.write(f"Metadata: {metadata_checksum}")
    

if __name__ == "__main__":
    concatenate_fasta_files()
    concatenate_excel_files()