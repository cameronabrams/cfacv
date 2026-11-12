#!/usr/bin/env bash
#
# String method in collective variables: RMSD of the string, frame by frame,
# from a reference string, for monitoring convergence.
#
# Reads the history file(s) that tcl/stringmethod.tcl writes (lines
# "SM iter <n> replica <i> : reparam z <z_1> ... <z_ncv>") and writes
# rmsd.dat, one line per frame:
#
#     <iteration> <RMSD>
#
# where RMSD = sqrt( (1/NI) sum_i |z_i - z0_i|^2 ) over the NI images, and the
# reference z0 is the average of the last <avg_last> frames (default 1: the
# final string). In dual-image mode only the first NI replicas, the images
# themselves, are used.
#
# Usage (from the run directory):
#
#     measurermsd.sh [-smc alad_stringmethod.conf] [-c cv.inp]
#                    [-nA avg_last] [-o rmsd.dat] [-smh history] [history ...]
#
# With no history file named, uses output/0/alad_sm.job*.0.history, which
# collects the histories of a run and its restarts in order.
#
# Plain bash and POSIX awk; no bc. Supersedes measurermsd_JMG.sh and the
# awk-based rmsd.sh formerly in sbin/.
#
# (c) 2016-2026 Cameron F Abrams, Drexel University; 2018 Jasmine M Gardner

set -euo pipefail

STRING_METHOD_CONFIG=alad_stringmethod.conf
CV_INP=cv.inp
AVG_LAST=1
OUT=rmsd.dat
HISTORIES=()

while [ $# -gt 0 ]; do
    case "$1" in
        -smc|--string_method_config) STRING_METHOD_CONFIG="$2"; shift 2 ;;
        -c|--cv_inp)                 CV_INP="$2"; shift 2 ;;
        -smh|--string_method_history) HISTORIES+=("$2"); shift 2 ;;
        -nA|--avg_last)              AVG_LAST="$2"; shift 2 ;;
        -o|--output)                 OUT="$2"; shift 2 ;;
        -h|--help)                   sed -n '2,27p' "$0"; exit 0 ;;
        -*) echo "ERROR: unknown option $1" >&2; exit 1 ;;
        *)  HISTORIES+=("$1"); shift ;;
    esac
done

if [ ${#HISTORIES[@]} -eq 0 ]; then
    shopt -s nullglob
    HISTORIES=(output/0/alad_sm.job*.0.history)
    shopt -u nullglob
fi
for f in "$STRING_METHOD_CONFIG" "$CV_INP" ${HISTORIES[@]+"${HISTORIES[@]}"}; do
    [ -f "$f" ] || { echo "ERROR: file $f not found." >&2; exit 1; }
done
[ ${#HISTORIES[@]} -gt 0 ] || { echo "ERROR: no history file found; name one." >&2; exit 1; }

N_SYSTEMS=$(awk '$1=="set" && $2=="num_replicas" {print $3}' "$STRING_METHOD_CONFIG")
DUAL=$(awk '$1=="set" && $2=="SMPARAMS(dual)" {print $3}' "$STRING_METHOD_CONFIG")
NCV=$(grep -cv '^[[:space:]]*\(#\|$\)' "$CV_INP")
NI=$N_SYSTEMS
[ "${DUAL:-0}" = "1" ] && NI=$((N_SYSTEMS / 2))

echo "measurermsd.sh: $NI images, $NCV CVs, reference = mean of last $AVG_LAST frame(s)"
echo "measurermsd.sh: history ${HISTORIES[*]}"

awk -v NI="$NI" -v NCV="$NCV" -v NA="$AVG_LAST" -v OUT="$OUT" '
    # "SM iter <n> replica <i> : reparam z v1 ... vNCV"; the run ends with
    # "SM FINALiter <n> ...", the final string
    $1 == "SM" && ($2 == "iter" || $2 == "FINALiter") && $6 == ":" && $7 == "reparam" && $8 == "z" {
        it = $3; img = $5
        if (img >= NI) next
        if (!(it in seen)) { seen[it] = 1; order[++nf] = it }
        for (k = 1; k <= NCV; k++) z[it, img, k] = $(8 + k)
        cnt[it]++
    }
    END {
        # use only frames that list every image
        n = 0
        for (f = 1; f <= nf; f++) if (cnt[order[f]] == NI) fr[++n] = order[f]
        if (n == 0) { print "ERROR: no complete frames" > "/dev/stderr"; exit 1 }
        na = (NA > n) ? n : NA
        for (i = 0; i < NI; i++) for (k = 1; k <= NCV; k++) {
            s = 0
            for (f = n - na + 1; f <= n; f++) s += z[fr[f], i, k]
            z0[i, k] = s / na
        }
        print "# iter RMSD" > OUT
        for (f = 1; f <= n; f++) {
            d = 0
            for (i = 0; i < NI; i++) for (k = 1; k <= NCV; k++) {
                e = z[fr[f], i, k] - z0[i, k]; d += e * e
            }
            printf "%d %.8f\n", fr[f], sqrt(d / NI) > OUT
        }
        printf "measurermsd.sh: %d frames written to %s\n", n, OUT
    }
' "${HISTORIES[@]}"
