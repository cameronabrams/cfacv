#!/bin/bash
#
# String method in collective variables via NAMD/replica
#
# Measurement of string RMSD
#
# Required files:
#   - cv.inp used by cfacv.tcl in the TAMD sweeps
#   - restr_initX.tmp - template restr.inp file for restrained MD via cfacv
#   - initX.conf - template NAMD configuration file for the restrained MD
#   - stringmethod.conf - string method configuration file used in the SMCV simulation
#
# (c) 2016-2018 Cameron F Abrams, Drexel University
#
# default values
cv_inp=cv.inp
string_method_history=output/0/alad_sm.job0.0.history
STRING_METHOD_CONFIG=alad_stringmethod.conf
CFACV_BASEDIR=${HOME}/research/cfacv

avg_last=1  # number of iterations to average over to generate anchor points
DIHED="-dihed" # set to "-dihed" if CV's are dihedral angles

while [[ $# -gt 1 ]]
do
key="$1"
case $key in
    -c|--cv_inp)
    cv_inp="$2"
    shift # past argument
    ;;
    -smc|--string_method_config)
    STRING_METHOD_CONFIG="$2"
    shift # past argument
    ;;
    -smh|--string_method_history)
    string_method_history="$2"
    shift # past argument
    ;;
    -nA|--avg_last)
    avg_last="$2"
    shift # past argument
    ;;
    -cbd|--cfacv_basedir)
    CFACV_BASEDIR="$2"
    shift
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

for f in $cv_inp $STRING_METHOD_CONFIG ; do
  if [ ! -f $f ]; then
    echo "ERROR: file $f not found."
    exit
  fi
done

for d in $CFACV_BASEDIR; do
  if [ ! -d $d ]; then
    echo "ERROR: directory $d not found."
    exit
  fi
done

N_SYSTEMS=`grep num_replicas $STRING_METHOD_CONFIG|awk '{print $3}'`
DUAL=`grep "SMPARAMS(dual)" $STRING_METHOD_CONFIG|awk '{print $3}'`
SYSTEMS_PER_IMAGE=1
if [ "$DUAL" -eq "1" ] ; then
  SYSTEMS_PER_IMAGE=2
fi

NI=`echo "$N_SYSTEMS / $SYSTEMS_PER_IMAGE" | bc`
NCV=`grep -v ^\# $cv_inp | wc -l | awk '{print $1}'`
NCVI=`echo "$NCV - 1" | bc`
CVL=`grep -v ^\# $cv_inp | awk '{print $1}'`
NII=`echo "$NI - 1" | bc`


cat ${string_method_history}*.history > tmp.hist

NI_avg_last=`echo "$NI*$avg_last" | bc`
echo $NI_avg_last
echo $avg_last
tac tmp.hist | grep reparam | awk -v NI=$NI_avg_last 'BEGIN {delimit="\n"} NR<=NI {newa=a delimit; a=newa $0; next} 1; END {print a}' | tac  > tmp1.hist

rm tmp.hist

awk -v NI=$NI_avg_last -v CV=$NCV -v avg_num=$avg_last '
BEGIN {
	printf "# iter RMSD \n"; }
{
	if (NR<=(NI-NI/avg_num)) {
		for (i=9; i<9+CV; i++) {
			Z0[$5][i-9]+=$i; 
	}} else if (NR<=NI) {
		for (i=9; i<9+CV; i++) {
			Z0[$5][i-9]+=$i; 
			Z0[$5][i-9]/=avg_num;
	}} else {
		for (i=9; i<9+CV; i++) {
			D[$3]+=($i-Z0[$5][i-9])*($i-Z0[$5][i-9]);
	}} 
	
	if ($5==(NI/avg_num-1) && NR>NI) {
		printf("%i %.8f \n", $3, sqrt(D[$3]/NI*avg_num));
}}' tmp1.hist > rmsd.dat

rm tmp1.hist
