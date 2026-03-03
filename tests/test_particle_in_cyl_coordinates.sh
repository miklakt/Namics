#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

binary="${repo_root}/bin/namics"
input_file="${repo_root}/tests/particle_in_cyl_coordinates.in"
reference_file="${repo_root}/tests/reference/particle_in_cyl_coordinates.json.ref"
output_dir="${repo_root}/output"
output_file_json="${output_dir}/particle_in_cyl_coordinates.json"
compare_script="${repo_root}/tests/compare_profile_content.py"
phi_tolerance="1e-6"
xy_tolerance="1e-12"

cleanup() {
  rm -f "${output_file_json}"
}
trap cleanup EXIT

if [[ ! -x "${binary}" ]]; then
  echo "ERROR: built binary not found or not executable: ${binary}" >&2
  echo "Build it first with: make" >&2
  exit 1
fi

for required_file in "${input_file}" "${reference_file}" "${compare_script}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "ERROR: required file not found: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p "${output_dir}"
rm -f "${output_file_json}"

"${binary}" "${input_file}" > /dev/null

if [[ ! -s "${output_file_json}" ]]; then
  echo "ERROR: expected JSON output file was not created: ${output_file_json}" >&2
  exit 1
fi

if ! python3 "${compare_script}" --left "${reference_file}" --right "${output_file_json}" --coord-tol "${xy_tolerance}" --value-tol "${phi_tolerance}"; then
  echo "ERROR: output differs from reference: ${reference_file}" >&2
  exit 1
fi

echo "PASS: particle_in_cyl_coordinates regression test"
