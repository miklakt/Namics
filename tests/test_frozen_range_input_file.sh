#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

binary="${repo_root}/bin/namics"
input_file="${repo_root}/tests/frozen_range_input_file.in"
reference_file="${repo_root}/tests/reference/frozen_range_input_file.pro.ref"
output_dir="${repo_root}/output"
output_file="${output_dir}/frozen_range_input_file.pro"
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

mkdir -p "${output_dir}"
rm -f "${output_file}"

"${binary}" "${input_file}" > /dev/null

if [[ ! -s "${output_file}" ]]; then
  echo "ERROR: expected output file was not created: ${output_file}" >&2
  exit 1
fi

if ! awk 'NF < 2 { exit 1 } END { if (NR < 2) exit 1 }' "${output_file}"; then
  echo "ERROR: output file is not tabulated: ${output_file}" >&2
  exit 1
fi

if ! awk -v reference="${reference_file}" -v phi_tol="${phi_tolerance}" -v xy_tol="${xy_tolerance}" '
function absval(x) {
  return x < 0 ? -x : x
}
BEGIN {
  n_ref = 0
  while ((getline line < reference) > 0) {
    ref_lines[++n_ref] = line
  }
  close(reference)
}
{
  if (NR > n_ref) {
    printf "ERROR: output has more lines than reference (extra line %d)\n", NR > "/dev/stderr"
    exit 1
  }

  n_out = split($0, out, "\t")
  n_expected = split(ref_lines[NR], expected, "\t")
  if (n_out < 5 || n_expected < 5) {
    printf "ERROR: malformed tabulated row at line %d\n", NR > "/dev/stderr"
    exit 1
  }

  if (NR == 1) {
    if (out[1] != expected[1] || out[2] != expected[2] || out[3] != expected[3] || out[4] != expected[4] || out[5] != expected[5]) {
      printf "ERROR: header mismatch at line 1\n" > "/dev/stderr"
      exit 1
    }
    next
  }

  if (absval((out[1] + 0) - (expected[1] + 0)) > xy_tol || absval((out[2] + 0) - (expected[2] + 0)) > xy_tol) {
    printf "ERROR: coordinate mismatch at line %d (x,y)\n", NR > "/dev/stderr"
    exit 1
  }

  if (absval((out[3] + 0) - (expected[3] + 0)) > phi_tol ||
      absval((out[4] + 0) - (expected[4] + 0)) > phi_tol ||
      absval((out[5] + 0) - (expected[5] + 0)) > phi_tol) {
    printf "ERROR: phi mismatch at line %d (tol=%s)\n", NR, phi_tol > "/dev/stderr"
    exit 1
  }
}
END {
  if (NR != n_ref) {
    printf "ERROR: output has fewer lines than reference (output=%d, reference=%d)\n", NR, n_ref > "/dev/stderr"
    exit 1
  }
}
' "${output_file}"; then
  echo "ERROR: output differs from reference (within tolerance): ${reference_file}" >&2
  exit 1
fi

echo "PASS: frozen range input file regression test"
