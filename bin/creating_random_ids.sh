#!/bin/bash
input_file=/home/glbrc.org/millican/ref_db/trait-files/pgpt-ontology_headers.txt
output_file=/home/glbrc.org/millican/ref_db/trait-files/pgpt-ontology_headers_cleaned.txt
final_file=/home/glbrc.org/millican/ref_db/trait-files/pgpt-ontology_uid_headers.txt

awk -F ">" '{print substr($2, 1)}' $input_file | awk '{sub(/-/, "\t"); print}' | awk '{sub(/-/, "\t"); print}' | awk '{sub(/ /, "\t"); print}' | awk '{gsub(/\//, "\t"); print}' > $output_file

# count the number of lines in the input file
num_lines=$(cat $output_file | wc -l)

# generate the same number of unique random strings
strings=$(python -c "import random,string;print('\n'.join(set([''.join(random.choice(string.ascii_letters) + ''.join(random.choices(string.ascii_letters + string.digits, k=7))) for _ in range($num_lines)])))")

# write the strings to a temporary file
tmp_file=$(mktemp)
echo "$strings" > $tmp_file

# join the random strings with the output file
paste $tmp_file $output_file > $final_file

# clean up the temporary file
rm $tmp_file