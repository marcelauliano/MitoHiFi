def modify_circularity_in_header(input_file, output_file):
    with open(input_file, 'r') as input_fasta, open(output_file, 'w') as output_fasta:
        for line in input_fasta:
            if line.startswith('>'):  # Check if the line is a header
                header = line.strip()[1:]  # Remove '>' and any leading/trailing whitespace
                seq_id = header.split()[0] # 
                is_circular = "suggestCircular=yes" in header  # Check if the suggestion is "yes"
                is_linear = "suggestCircular=no" in header  # Check if the suggestion is "no"
                if is_circular:
                    new_header = f">{seq_id}c"
                elif is_linear:
                    new_header = f">{seq_id}l"
                else:
                    new_header = f">{header}"
                output_fasta.write(new_header + "\n")
            else:
                output_fasta.write(line.strip() + "\n")
