#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

binary="${repo_root}/bin/namics"
generate_input="${repo_root}/tests/micelle_guess_generate.in"
use_input="${repo_root}/tests/micelle_guess_use.in"
reference_file="${repo_root}/tests/reference/micelle_guess_use.pro.ref"
output_dir="${repo_root}/output"
guess_file="${output_dir}/micelle_2.outi"
output_file_json="${output_dir}/micelle_guess_use.json"
compare_script="${repo_root}/tests/compare_profile_content.py"
x_tolerance="1e-12"
phi_tolerance="1e-9"

cleanup() {
  rm -f "${guess_file}" "${output_file_json}"
}
trap cleanup EXIT

if [[ ! -x "${binary}" ]]; then
  echo "ERROR: built binary not found or not executable: ${binary}" >&2
  echo "Build it first with: make" >&2
  exit 1
fi

for required_file in "${generate_input}" "${use_input}" "${reference_file}" "${compare_script}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "ERROR: required file not found: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p "${output_dir}"
rm -f "${guess_file}" "${output_file_json}"

"${binary}" "${generate_input}" > /dev/null

if [[ ! -s "${guess_file}" ]]; then
  echo "ERROR: expected guess file was not created: ${guess_file}" >&2
  exit 1
fi

"${binary}" "${use_input}" > /dev/null

if [[ ! -s "${output_file_json}" ]]; then
  echo "ERROR: expected JSON output file was not created: ${output_file_json}" >&2
  exit 1
fi

if ! python3 "${compare_script}" --left "${reference_file}" --right "${output_file_json}" --coord-tol "${x_tolerance}" --value-tol "${phi_tolerance}"; then
  echo "ERROR: output differs from reference: ${reference_file}" >&2
  exit 1
fi

echo "PASS: micelle self-assembly regression test"
