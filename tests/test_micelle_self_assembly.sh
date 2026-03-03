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
output_file="${output_dir}/micelle_guess_use.pro"
x_tolerance="1e-12"
phi_tolerance="1e-9"

compare_tabulated() {
  local reference="$1"
  local output="$2"
  local x_tol="$3"
  local value_tol="$4"

  awk -v reference="${reference}" -v x_tol="${x_tol}" -v value_tol="${value_tol}" '
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

  if (NR == 1) {
    if ($0 != ref_lines[NR]) {
      printf "ERROR: header mismatch at line 1\n" > "/dev/stderr"
      exit 1
    }
    next
  }

  if (n_out != n_expected || n_out < 4) {
    printf "ERROR: malformed or inconsistent row at line %d\n", NR > "/dev/stderr"
    exit 1
  }

  if (absval((out[1] + 0) - (expected[1] + 0)) > x_tol) {
    printf "ERROR: x mismatch at line %d\n", NR > "/dev/stderr"
    exit 1
  }

  for (i = 2; i <= n_out; i++) {
    if (absval((out[i] + 0) - (expected[i] + 0)) > value_tol) {
      printf "ERROR: value mismatch at line %d, col %d (tol=%s)\n", NR, i, value_tol > "/dev/stderr"
      exit 1
    }
  }
}
END {
  if (NR != n_ref) {
    printf "ERROR: output has fewer lines than reference (output=%d, reference=%d)\n", NR, n_ref > "/dev/stderr"
    exit 1
  }
}
' "${output}"
}

if [[ ! -x "${binary}" ]]; then
  echo "ERROR: built binary not found or not executable: ${binary}" >&2
  echo "Build it first with: make" >&2
  exit 1
fi

for required_file in "${generate_input}" "${use_input}" "${reference_file}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "ERROR: required file not found: ${required_file}" >&2
    exit 1
  fi
done

mkdir -p "${output_dir}"
rm -f "${guess_file}" "${output_file}"

"${binary}" "${generate_input}" > /dev/null

if [[ ! -s "${guess_file}" ]]; then
  echo "ERROR: expected guess file was not created: ${guess_file}" >&2
  exit 1
fi

"${binary}" "${use_input}" > /dev/null

if [[ ! -s "${output_file}" ]]; then
  echo "ERROR: expected output file was not created: ${output_file}" >&2
  exit 1
fi

if ! awk 'NF < 2 { exit 1 } END { if (NR < 2) exit 1 }' "${output_file}"; then
  echo "ERROR: output file is not tabulated: ${output_file}" >&2
  exit 1
fi

if ! compare_tabulated "${reference_file}" "${output_file}" "${x_tolerance}" "${phi_tolerance}"; then
  echo "ERROR: output differs from reference: ${reference_file}" >&2
  exit 1
fi

echo "PASS: micelle self-assembly regression test"
