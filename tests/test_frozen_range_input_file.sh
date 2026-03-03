#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

binary="${repo_root}/bin/namics"
input_file="${repo_root}/tests/frozen_range_input_file.in"
reference_file="${repo_root}/tests/reference/frozen_range_input_file.pro.ref"
output_dir="${repo_root}/output"
output_file_pro="${output_dir}/frozen_range_input_file.pro"
output_file_json="${output_dir}/frozen_range_input_file.json"
compare_script="${repo_root}/tests/compare_profile_content.py"
phi_tolerance="1e-6"
xy_tolerance="1e-12"

if [[ ! -x "${binary}" ]]; then
  echo "ERROR: built binary not found or not executable: ${binary}" >&2
  echo "Build it first with: make" >&2
  exit 1
fi

if [[ ! -f "${input_file}" ]]; then
  echo "ERROR: input file not found: ${input_file}" >&2
  exit 1
fi

if [[ ! -f "${reference_file}" ]]; then
  echo "ERROR: reference file not found: ${reference_file}" >&2
  exit 1
fi

if [[ ! -f "${compare_script}" ]]; then
  echo "ERROR: compare script not found: ${compare_script}" >&2
  exit 1
fi

mkdir -p "${output_dir}"
rm -f "${output_file_pro}" "${output_file_json}"

"${binary}" "${input_file}" > /dev/null

output_file=""
if [[ -s "${output_file_json}" ]]; then
  output_file="${output_file_json}"
elif [[ -s "${output_file_pro}" ]]; then
  output_file="${output_file_pro}"
fi

if [[ -z "${output_file}" ]]; then
  echo "ERROR: expected output file was not created: ${output_file_json} or ${output_file_pro}" >&2
  exit 1
fi

if ! python3 "${compare_script}" --left "${reference_file}" --right "${output_file}" --coord-tol "${xy_tolerance}" --value-tol "${phi_tolerance}"; then
  echo "ERROR: output differs from reference (within tolerance): ${reference_file}" >&2
  exit 1
fi

echo "PASS: frozen range input file regression test"
