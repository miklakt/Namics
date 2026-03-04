#!/usr/bin/env bash

# Render INPUT to OUTPUT while forcing the solver method line to METHOD.
# If the line is missing, insert it before the first "start" block.
ensure_solver_method_line() {
  local input_file="$1"
  local output_file="$2"
  local method="$3"

  awk -v method="${method}" '
    {
      lines[NR] = $0
      if ($0 ~ /^[[:space:]]*newton[[:space:]]*:[[:space:]]*isaac[[:space:]]*:[[:space:]]*method[[:space:]]*:/) {
        has_method = 1
      }
      if (first_start == 0 && $0 ~ /^[[:space:]]*start[[:space:]]*$/) {
        first_start = NR
      }
    }
    END {
      for (i = 1; i <= NR; i++) {
        if (lines[i] ~ /^[[:space:]]*newton[[:space:]]*:[[:space:]]*isaac[[:space:]]*:[[:space:]]*method[[:space:]]*:/) {
          print "newton : isaac : method : " method
        } else {
          if (!has_method && first_start > 0 && i == first_start) {
            print "newton : isaac : method : " method
          }
          print lines[i]
        }
      }
      if (!has_method && first_start == 0) {
        print "newton : isaac : method : " method
      }
    }
  ' "${input_file}" > "${output_file}"
}
