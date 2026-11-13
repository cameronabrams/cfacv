#!/usr/bin/env bash
#
# String method in collective variables: free energy along the converged
# string, from the mean force measured at each image.
#
# For each image i of the string, this script
#   1. takes the image's CV values z_i from the string-method history (the
#      average of the last <avg_last> frames; default 1, the final string),
#   2. writes fep_<i>.conf from the NAMD template (%IMG% -> i) and the
#      restraint file it names from the restraint template (%K% -> k,
#      %CV0%, %CV1%, ... -> z_i),
#   3. runs restrained MD on it with NAMD, writing fep_<i>.log, and
#   4. extracts the mean force from the log with bin/ForcesFromLog,
# then integrates the mean force along the string with bin/FEPFromForces,
# writing fep.dat (columns: image, arc length, free energy | z, f, tangent).
#
# Usage (from the run directory):
#
#     measurefep.sh [-np 8] [-namd namd3] [-cr charmrun] [-k 100.0]
#                   [-nA avg_last] [-c cv.inp] [-rX restr_fepX.inp]
#                   [-nX fepX.conf] [-smc alad_stringmethod.conf]
#                   [-smh history] [-cbd cfacv_basedir] [-nodihed]
#                   [-n] [-force]
#
#   -np N      processors for each NAMD run (default 8)
#   -namd EXE  NAMD executable (default namd3; also -namd2)
#   -cr EXE    launch NAMD through this charmrun (default: none; NAMD is
#              run as "EXE +pN"; with -cr, as "charmrun +pN EXE")
#   -k K       spring constant of the restraints, substituted for %K%
#   -nA N      average the image positions over the last N frames
#   -smh FILE  string-method history; may be repeated (default:
#              output/0/alad_sm.job*.0.history, the run and its restarts)
#   -cbd DIR   the cfacv tree (default: $CFACV_BASEDIR, else the tree this
#              script is in)
#   -nodihed   the CVs are not dihedral angles (default: they are, and
#              ForcesFromLog wraps their differences into [-pi, pi])
#   -n         write the configuration and restraint files, run nothing
#   -force     rerun images whose fep_<i>.log already ended normally
#
# An image whose fep_<i>.log already ends with NAMD's "End of program" is not
# rerun, so an interrupted measurement resumes where it stopped.
#
# (c) 2016-2026 Cameron F Abrams, Drexel University

set -euo pipefail

CV_INP=cv.inp
RESTR_TEMPLATE=restr_fepX.inp
NAMD_TEMPLATE=fepX.conf
STRING_METHOD_CONFIG=alad_stringmethod.conf
HISTORIES=()
AVG_LAST=1
NP=8
K=100.0
NAMD=namd3
CHARMRUN=
DIHED=-dihed
DRY_RUN=0
FORCE=0
BASEDIR=${CFACV_BASEDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}

die() { echo "measurefep.sh: ERROR: $*" >&2; exit 1; }

while [ $# -gt 0 ]; do
    case "$1" in
        -c|--cv_inp)                  CV_INP="$2"; shift 2 ;;
        -rX|--restr_inp_template)     RESTR_TEMPLATE="$2"; shift 2 ;;
        -nX|--namd_config_template)   NAMD_TEMPLATE="$2"; shift 2 ;;
        -smc|--string_method_config)  STRING_METHOD_CONFIG="$2"; shift 2 ;;
        -smh|--string_method_history) HISTORIES+=("$2"); shift 2 ;;
        -nA|--avg_last)               AVG_LAST="$2"; shift 2 ;;
        -cbd|--cfacv_basedir)         BASEDIR="$2"; shift 2 ;;
        -namd|-namd2|--namd)          NAMD="$2"; shift 2 ;;
        -cr|--charmrun)               CHARMRUN="$2"; shift 2 ;;
        -np|--number_of_processors)   NP="$2"; shift 2 ;;
        -k|--spring_constant)         K="$2"; shift 2 ;;
        -nodihed)                     DIHED=; shift ;;
        -n|--dry-run)                 DRY_RUN=1; shift ;;
        -force)                       FORCE=1; shift ;;
        -h|--help)                    sed -n '2,44p' "$0"; exit 0 ;;
        *) die "unknown option $1 (try -h)" ;;
    esac
done

if [ ${#HISTORIES[@]} -eq 0 ]; then
    shopt -s nullglob
    HISTORIES=(output/0/alad_sm.job*.0.history)
    shopt -u nullglob
fi
[ ${#HISTORIES[@]} -gt 0 ] || die "no history file found; name one with -smh."
for f in "$CV_INP" "$RESTR_TEMPLATE" "$NAMD_TEMPLATE" "$STRING_METHOD_CONFIG" "${HISTORIES[@]}"; do
    [ -f "$f" ] || die "file $f not found."
done
for t in ForcesFromLog FEPFromForces; do
    [ -x "$BASEDIR/bin/$t" ] || die "$BASEDIR/bin/$t not found; build it (cd src; make $t) or give the cfacv tree with -cbd."
done
if [ $DRY_RUN -eq 0 ]; then
    command -v "$NAMD" >/dev/null || die "NAMD executable $NAMD not found; give it with -namd."
    [ -z "$CHARMRUN" ] || command -v "$CHARMRUN" >/dev/null || die "charmrun $CHARMRUN not found."
fi

N_SYSTEMS=$(awk '$1=="set" && $2=="num_replicas" {print $3}' "$STRING_METHOD_CONFIG")
DUAL=$(awk '$1=="set" && $2=="SMPARAMS(dual)" {print $3}' "$STRING_METHOD_CONFIG")
NCV=$(grep -cv '^[[:space:]]*\(#\|$\)' "$CV_INP")
[ -n "$N_SYSTEMS" ] || die "no 'set num_replicas' in $STRING_METHOD_CONFIG."
NI=$N_SYSTEMS
[ "${DUAL:-0}" = "1" ] && NI=$((N_SYSTEMS / 2))

echo "measurefep.sh: $NI images, $NCV CVs, k = $K, image positions = mean of last $AVG_LAST frame(s)"
echo "measurefep.sh: history ${HISTORIES[*]}"

# One line per image: <image> <z_1> ... <z_NCV>. Frames are read as in
# measurermsd.sh: "SM iter|FINALiter <n> replica <i> : reparam z ...", only
# the first NI replicas, only frames that list every image.
IMAGE_Z=$(awk -v NI="$NI" -v NCV="$NCV" -v NA="$AVG_LAST" '
    $1 == "SM" && ($2 == "iter" || $2 == "FINALiter") && $6 == ":" && $7 == "reparam" && $8 == "z" {
        it = $3; img = $5
        if (img >= NI) next
        if (!(it in seen)) { seen[it] = 1; order[++nf] = it }
        for (k = 1; k <= NCV; k++) z[it, img, k] = $(8 + k)
        cnt[it]++
    }
    END {
        n = 0
        for (f = 1; f <= nf; f++) if (cnt[order[f]] == NI) fr[++n] = order[f]
        if (n == 0) { print "measurefep.sh: ERROR: no complete frames in the history" > "/dev/stderr"; exit 1 }
        na = (NA > n) ? n : NA
        for (i = 0; i < NI; i++) {
            line = i
            for (k = 1; k <= NCV; k++) {
                s = 0
                for (f = n - na + 1; f <= n; f++) s += z[fr[f], i, k]
                line = line " " sprintf("%.17g", s / na)
            }
            print line
        }
    }
' "${HISTORIES[@]}") || die "could not read image positions from the history."

# fepX.conf finds the cfacv tree through CFACV_BASEDIR
export CFACV_BASEDIR=$BASEDIR

if [ -n "$CHARMRUN" ]; then LAUNCH=("$CHARMRUN" "+p$NP" "$NAMD"); else LAUNCH=("$NAMD" "+p$NP"); fi

IMAGE_FORCES=image_forces.dat
if [ $DRY_RUN -eq 0 ]; then
    {
        echo "# image mean forces from restrained MD with k = $K"
        echo "# MD engine: ${LAUNCH[*]}"
        echo "# ---z----   -----f-----"
    } > "$IMAGE_FORCES"
fi

while read -r img zs; do
    sed "s/%IMG%/$img/g" "$NAMD_TEMPLATE" > "fep_${img}.conf"
    rfn=$(awk '$1=="set" && $2=="restrINP" {print $3}' "fep_${img}.conf")
    [ -n "$rfn" ] || die "fep_${img}.conf sets no restrINP."
    expr="s/%K%/$K/g"
    i=0
    for z in $zs; do expr="$expr; s/%CV${i}%/$z/g"; i=$((i + 1)); done
    sed "$expr" "$RESTR_TEMPLATE" > "$rfn"
    if grep -q '%CV[0-9]*%\|%K%' "$rfn"; then die "$rfn still has unsubstituted placeholders."; fi
    echo "Image $img : z $zs -> fep_${img}.conf, $rfn"
    [ $DRY_RUN -eq 1 ] && continue

    log="fep_${img}.log"
    if [ $FORCE -eq 0 ] && [ -f "$log" ] && grep -q "End of program" "$log"; then
        echo "  $log already complete; not rerun (use -force to rerun)"
    else
        echo "  running: ${LAUNCH[*]} fep_${img}.conf > $log"
        "${LAUNCH[@]}" "fep_${img}.conf" > "$log" 2>&1 || die "NAMD failed on image $img; see $log."
        grep -q "End of program" "$log" || die "NAMD did not finish image $img; see $log."
    fi
    "$BASEDIR/bin/ForcesFromLog" -k "$K" -f "$log" -dim "$NCV" -ofn "fep_fra_${img}.dat" $DIHED >> "$IMAGE_FORCES" \
        || die "ForcesFromLog failed on $log."
done <<< "$IMAGE_Z"

[ $DRY_RUN -eq 1 ] && { echo "measurefep.sh: dry run; no NAMD runs."; exit 0; }

"$BASEDIR/bin/FEPFromForces" -i "$IMAGE_FORCES" -dim "$NCV" -o fep.dat
echo "measurefep.sh: fep.dat created."
