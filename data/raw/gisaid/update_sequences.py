import subprocess
import glob
import os
import hashlib

def concatenate_fasta_files():
    # Create the output filename
    output_file = "H5N1_concatenated.fasta"

    # Get all .fasta and .fa files in the current directory
    fasta_files = [file for file in glob.glob("*.fasta") + glob.glob("*.fa") if file != output_file]
    
    if not fasta_files:
        print("No FASTA files found in the current directory.")
        return
    
    # Sort the files to ensure consistent order
    fasta_files.sort()
    
    # Use cat command to concatenate all files
    cat_command = f"cat {' '.join(fasta_files)} > {output_file}"
    
    try:
        # Execute the cat command
        subprocess.run(cat_command, shell=True, check=True)
        print(f"Successfully concatenated {len(fasta_files)} FASTA files into {output_file}")
        
        # Get the size of the output file
        file_size = os.path.getsize(output_file)
        print(f"Size of {output_file}: {file_size} bytes")

        # Calculate the MD5 checksum of the output file
        with open(output_file, "rb") as file:
            file_content = file.read()
            md5_checksum = hashlib.md5(file_content).hexdigest()

        print(f"MD5 checksum of latest update is:\n\t{md5_checksum}")

        checksum_file = f"checksum.txt"

        with open(checksum_file, "w") as file:
            file.write(md5_checksum)

        print(f"MD5 checksum saved in {checksum_file}")
        
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while concatenating files: {e}")

if __name__ == "__main__":
    concatenate_fasta_files()