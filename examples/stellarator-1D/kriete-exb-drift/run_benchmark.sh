#!/usr/bin/env bash
# Run the binormal ExB drift benchmark (Kriete et al. NF 2023) and compare
# with the analytic solution.
#
# Usage: ./run_benchmark.sh /path/to/hermes-3 [output_dir]
#
# Runs:
#   - gamma scan on the parallel-length grid (ny = 200)
#   - resolution scan at gamma = 0.25 (ny = 100, 200, 400, 800)
#   - gamma = 0.25 on an island-poloidal grid (y = s_theta, g^22 = Theta^2)
#   - gamma = 0.25 in Kriete's geometry with source-free legs (y_xpoint = 0.4)
#
# Parameters are written into each run's BOUT.inp (not passed on the command
# line) so that runs can be restarted consistently.

set -euo pipefail

HERMES=${1:?"Usage: $0 /path/to/hermes-3 [output_dir]"}
OUT=${2:-runs}
HERE=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$OUT"

# Write BOUT.inp for a run: run_input NAME GAMMA NY Y_XPOINT [theta]
run_input() {
  local name=$1 gamma=$2 ny=$3 yx=$4 grid=${5:-par}
  mkdir -p "$OUT/$name"
  sed -e "s/^gamma = .*/gamma = $gamma/" \
      -e "s/^ny = .*/ny = $ny/" \
      -e "s/^y_xpoint = .*/y_xpoint = $yx/" \
      "$HERE/BOUT.inp" > "$OUT/$name/BOUT.inp"
  if [ "$grid" = "theta" ]; then
    # Island-poloidal grid: the 1D coordinate is y itself and Theta enters
    # through the metric: Grad_par = Theta d/dy, Div_par(F) = d/dy(Theta F)
    sed -i -e "s|^length = .*|length = 2 * y_target\ng22 = Theta^2|" \
           -e "s|^ypol = .*|ypol = ypos - length / 2|" "$OUT/$name/BOUT.inp"
  fi
}

run() {
  echo "Running $1"
  (cd "$OUT/$1" && "$HERMES" -d . -q -q > run.log 2>&1)
}

for gamma in 0.0 0.25 -0.25 0.4 0.5; do
  run_input "gamma_${gamma}" $gamma 200 0.5
  run "gamma_${gamma}"
done

for ny in 100 400 800; do
  run_input "gamma_0.25_ny${ny}" 0.25 $ny 0.5
  run "gamma_0.25_ny${ny}"
done

run_input gamma_0.25_theta_grid 0.25 200 0.5 theta
run gamma_0.25_theta_grid

run_input gamma_0.25_legs 0.25 400 0.4
run gamma_0.25_legs

python3 "$HERE/compare.py" "$OUT"/gamma_0.0 "$OUT"/gamma_0.25 "$OUT"/gamma_-0.25 \
  "$OUT"/gamma_0.4 "$OUT"/gamma_0.5 --plot "$OUT/kriete_gamma_scan.png"
python3 "$HERE/compare.py" "$OUT"/gamma_0.25_ny100 "$OUT"/gamma_0.25 \
  "$OUT"/gamma_0.25_ny400 "$OUT"/gamma_0.25_ny800 "$OUT"/gamma_0.25_theta_grid \
  "$OUT"/gamma_0.25_legs
