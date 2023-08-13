convert_fortran_to_json() {
    local input_file=$1
    local output_file=$2
    local current_section=""
    local json_output="{\n"

    # Read the input file line by line
    while IFS= read -r line; do
        # Skip empty lines and lines starting with !
       [[ -z "$line" || ${line} =~ ^[[:space:]]*! ]] && continue

        # Check if the line starts a new section
        if [[ "$line" =~ ^[[:space:]]*\& ]]; then
            current_section="${line//[& ]/}"
            continue
        fi

        # If the line contains "/END", it signifies the end of a section
        if [[ "$line" =~ \/END ]]; then
            current_section=""
            continue
        fi

        # If we're in a section, split the variable assignments and process each one
        if [[ -n "$current_section" ]]; then
            # Remove trailing commas
            line=${line%,}
            IFS=',' read -ra assignments <<< "$line"
            for assignment in "${assignments[@]}"; do
                # Separate variable and value
                key=$(echo "$assignment" | cut -d'=' -f1 | xargs)
                value=$(echo "$assignment" | cut -d'=' -f2 | xargs)

                # Append to the JSON string
                json_output+="\"$current_section.$key\": \"$value\",\n"
            done
        fi
    done < "$input_file"

    # Remove the last comma and close the JSON
    json_output=${json_output%,}
    json_output+="\n}"

    # Write to output file
    echo -e "$json_output" > "$output_file"
}

# Call the function with the input and output file paths
convert_fortran_to_json "INPUT_ORG" "output.json"
