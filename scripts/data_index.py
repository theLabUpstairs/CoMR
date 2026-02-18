import random
import string
from Bio import SeqIO

input_file = snakemake.input[0]  # Input FASTA file
output_file_1 = snakemake.output.output1  # Output FASTA file
output_file_2 = snakemake.output.output2  # Mapping file

def generate_random_id(length=16, existing_ids=None):
    """Generate a unique random string of upper case letters and digits."""
    if existing_ids is None:
        existing_ids = set()
    while True:
        new_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))
        if new_id not in existing_ids:
            existing_ids.add(new_id)
            return new_id

def sanitize_header(header):
    """Replace spaces and commas in headers with underscores."""
    return header.replace(" ", "_").replace(",", "_")  # Replace spaces and commas with underscores

def index_fasta(input_file, output_file, mapping_file):
    """Read a FASTA file, assign unique new random headers, write to a new file, and save the mapping."""
    existing_ids = set()  # Keep track of generated IDs to ensure uniqueness
    
    # Open files and write the header explicitly
    with open(output_file, 'w') as output_handle, open(mapping_file, 'w') as mapping_handle:
        print("Writing header to mapping file...")
        mapping_handle.write("Index,RealSeqName\n")
        
        # Iterate over the FASTA file
        for record in SeqIO.parse(input_file, "fasta"):
            new_id = generate_random_id(existing_ids=existing_ids)
            sanitized_description = sanitize_header(record.description)  # Sanitize the entire description
            # Write the mapping of new_id to the sanitized original sequence name
            mapping_handle.write(f"{new_id},{sanitized_description}\n")
            # Update the record ID for the FASTA output
            record.id = new_id
            record.description = ''  # Clear description to only have the new ID
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")

        print("Finished writing mapping file and output FASTA.")

# Main script execution
with open(snakemake.log[0], "w") as log_file:
    sys.stdout = log_file
    sys.stderr = log_file
    print("Starting FASTA indexing...", flush=True)
    index_fasta(input_file, output_file_1, output_file_2)
    print("FASTA indexing done.", flush=True)
