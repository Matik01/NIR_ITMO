#!/bin/bash

# Define paths
smi_file="experiment.smi"
workspace_dir="../workspace"

# Extract CHEMBL IDs from experiment.smi
smi_ids=$(awk '{print $2}' "$smi_file" | sort)

# Extract CHEMBL IDs from workspace directories
workspace_ids=$(ls "$workspace_dir" | awk -F@ '{print $1}' | sort)

# Find missing CHEMBL IDs
missing_ids=$(comm -23 <(echo "$smi_ids") <(echo "$workspace_ids"))

# Output missing CHEMBL IDs
echo "Missing CHEMBL IDs:"
echo "$missing_ids"
