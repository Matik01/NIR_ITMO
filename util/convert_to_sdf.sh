#!/bin/bash
#RUN FROM RECEPTOR_GROUP DIRECTORY
source_dir="."

for package_dir in "${source_dir}"/*_package; do
  if [[ -d "$package_dir" ]]; then
    receptor_id=$(basename "$package_dir" | sed 's/_package//')

    find "${package_dir}/data" -type f -name "*.smi" -exec cat {} + > temp.smi

    obabel -ismi temp.smi -osdf --gen2D | obabel -isdf --gen3D -osdf > "${receptor_id}_lig.sdf"

    mv "${receptor_id}_lig.sdf" "${source_dir}/"
  fi
done

rm -f temp.smi

echo "Complete!"

