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
benchmark_dir="${repo_root}/tests/benchmarks"
benchmark_history="${benchmark_dir}/homopolymer_adsorption_test_benchmark.csv"
benchmark_threshold_pct="15"
run_id="$(date -u +%Y%m%dT%H%M%S%3NZ)"
run_log="${benchmark_dir}/homopolymer_adsorption_test_${run_id}.log"
chi_values=(0 -2 -4 -6)

now_ms() {
  echo $(( $(date +%s%N) / 1000000 ))
}

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
mkdir -p "${benchmark_dir}"

benchmark_start_ms="$(now_ms)"
solver_runtime_ms=0

{
  echo "run_id=${run_id}"
  echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "threshold_pct=${benchmark_threshold_pct}"
  echo "cases=chi_Si(0,-2,-4,-6) x save_memory(off,on)"
} > "${run_log}"

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
  run_start_ms="$(now_ms)"
  "${binary}" "${runtime_input}" > /dev/null
  run_elapsed_ms=$(( $(now_ms) - run_start_ms ))
  solver_runtime_ms=$(( solver_runtime_ms + run_elapsed_ms ))
  echo "chi_Si=${chi},mode=baseline,elapsed_ms=${run_elapsed_ms}" >> "${run_log}"

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
  run_start_ms="$(now_ms)"
  "${binary}" "${runtime_input_save_memory}" > /dev/null
  run_elapsed_ms=$(( $(now_ms) - run_start_ms ))
  solver_runtime_ms=$(( solver_runtime_ms + run_elapsed_ms ))
  echo "chi_Si=${chi},mode=save_memory,elapsed_ms=${run_elapsed_ms}" >> "${run_log}"

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

benchmark_wall_ms=$(( $(now_ms) - benchmark_start_ms ))
previous_solver_runtime_ms="$(awk -F, 'NR > 1 { last=$2 } END { if (last != "") print last }' "${benchmark_history}" 2>/dev/null || true)"
delta_pct="NA"
performance_state="baseline"

if [[ -n "${previous_solver_runtime_ms}" ]] && [[ "${previous_solver_runtime_ms}" =~ ^[0-9]+$ ]] && (( previous_solver_runtime_ms > 0 )); then
  delta_pct="$(awk -v current="${solver_runtime_ms}" -v previous="${previous_solver_runtime_ms}" \
    'BEGIN { printf "%.2f", ((current - previous) * 100.0) / previous }')"
  is_significant="$(awk -v delta="${delta_pct}" -v threshold="${benchmark_threshold_pct}" \
    'BEGIN { abs_delta = (delta < 0 ? -delta : delta); print (abs_delta >= threshold ? 1 : 0) }')"
  if [[ "${is_significant}" -eq 1 ]]; then
    if (( solver_runtime_ms > previous_solver_runtime_ms )); then
      performance_state="slower"
    else
      performance_state="faster"
    fi
  else
    performance_state="stable"
  fi
fi

if [[ ! -f "${benchmark_history}" ]]; then
  echo "run_id,solver_runtime_ms,wall_runtime_ms,threshold_pct,previous_solver_runtime_ms,delta_pct,performance_state" > "${benchmark_history}"
fi
echo "${run_id},${solver_runtime_ms},${benchmark_wall_ms},${benchmark_threshold_pct},${previous_solver_runtime_ms:-NA},${delta_pct},${performance_state}" >> "${benchmark_history}"

{
  echo "solver_runtime_ms=${solver_runtime_ms}"
  echo "wall_runtime_ms=${benchmark_wall_ms}"
  echo "previous_solver_runtime_ms=${previous_solver_runtime_ms:-NA}"
  echo "delta_pct=${delta_pct}"
  echo "performance_state=${performance_state}"
  echo "history_file=${benchmark_history}"
} >> "${run_log}"

echo "Benchmark: solver=${solver_runtime_ms} ms, wall=${benchmark_wall_ms} ms"
echo "Benchmark log: ${run_log}"
echo "Benchmark history: ${benchmark_history}"

echo "PASS: homopolymer adsorption test (chi_Si=0,-2,-4,-6, with/without save_memory)"
