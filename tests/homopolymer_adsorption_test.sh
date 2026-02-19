#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

binary="${repo_root}/bin/namics"
template_file="${repo_root}/tests/homopolymer_adsorption.in"
output_dir="${repo_root}/output"
reference_dir="${repo_root}/tests/reference"
runtime_input="${output_dir}/homopolymer_adsorption.in"
runtime_input_save_memory="${output_dir}/homopolymer_adsorption_save_memory.in"
chi_values=(0 -2 -4 -6)

if [[ ! -x "${binary}" ]]; then
  echo "ERROR: built binary not found or not executable: ${binary}" >&2
  echo "Build it first with: make" >&2
  exit 1
fi

if [[ ! -f "${template_file}" ]]; then
  echo "ERROR: input template file not found: ${template_file}" >&2
  exit 1
fi

# 1) Ensure output folder exists.
mkdir -p "${output_dir}"

for chi in "${chi_values[@]}"; do
  reference_file="${reference_dir}/homopolymer_adsorption.chi_${chi}.pro.ref"
  output_file="${output_dir}/homopolymer_adsorption.pro"
  output_file_save_memory="${output_dir}/homopolymer_adsorption_save_memory.pro"

  if [[ ! -f "${reference_file}" ]]; then
    echo "ERROR: reference file not found: ${reference_file}" >&2
    exit 1
  fi

  # 2) Render input and run built namics.
  sed "s/{chi_Si}/${chi}/g" "${template_file}" > "${runtime_input}"
  rm -f "${output_file}"
  "${binary}" "${runtime_input}" > /dev/null

  # 3) Verify a tabulated output exists.
  if [[ ! -s "${output_file}" ]]; then
    echo "ERROR: expected output file was not created for chi_Si=${chi}: ${output_file}" >&2
    exit 1
  fi

  if ! awk 'NF < 2 { exit 1 } END { if (NR < 2) exit 1 }' "${output_file}"; then
    echo "ERROR: output file is not tabulated for chi_Si=${chi}: ${output_file}" >&2
    exit 1
  fi

  # 4) Compare run output with the stored reference snapshot.
  if ! diff -u "${reference_file}" "${output_file}"; then
    echo "ERROR: output differs from reference for chi_Si=${chi}: ${reference_file}" >&2
    exit 1
  fi

  # 5) Repeat with save_memory enabled and verify exact same result.
  sed \
    -e "s/{chi_Si}/${chi}/g" \
    -e 's#^//mol : pol : save_memory : true#mol : pol : save_memory : true#' \
    "${template_file}" > "${runtime_input_save_memory}"
  rm -f "${output_file_save_memory}"
  "${binary}" "${runtime_input_save_memory}" > /dev/null

  if [[ ! -s "${output_file_save_memory}" ]]; then
    echo "ERROR: expected save_memory output was not created for chi_Si=${chi}: ${output_file_save_memory}" >&2
    exit 1
  fi

  if ! awk 'NF < 2 { exit 1 } END { if (NR < 2) exit 1 }' "${output_file_save_memory}"; then
    echo "ERROR: save_memory output is not tabulated for chi_Si=${chi}: ${output_file_save_memory}" >&2
    exit 1
  fi

  if ! diff -u "${reference_file}" "${output_file_save_memory}"; then
    echo "ERROR: save_memory output differs from reference for chi_Si=${chi}: ${reference_file}" >&2
    exit 1
  fi

  if ! diff -u "${output_file}" "${output_file_save_memory}"; then
    echo "ERROR: save_memory output differs from baseline output for chi_Si=${chi}" >&2
    exit 1
  fi
done

echo "PASS: homopolymer adsorption test (chi_Si=0,-2,-4,-6, with/without save_memory)"
