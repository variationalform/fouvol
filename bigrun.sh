#!/bin/bash    

# """
# Coordinate the running of fouvol.py in a bash shell. See for example,
# https://github.com/variationalform/fouvol
# 
# Copyright (c) 2020, Simon Shaw
# (https://github.com/variationalform, https://www.brunel.ac.uk/people/simon-shaw).
# The moral right of the author has been asserted.
# 
# These codes are free software; you can redistribute them and/or
# modify them under the terms of the GNU General Public License Version 3 - the terms
# of which should accompany this source code.
# """

# Examples...
# sudo apt-get install python-numpy python-matplotlib python-scipy dvipng
# time ./bigrun.sh -A "-v 0 -s 2 -a -0.5 -T 10 "  -J 8 17  -L 0 0 -C | tee -a runout.txt

# prevent formation of core file (takes ages), doesn't seem to work!
ulimit -c 0 # maybe do this in docker

# set up the defaults
STDARGS="-v 0 -s 3 -a -0.5 -m 3 --T1 0.1 "
MINJ=0
MAXJ=4
MINL=0
MAXL=8
DL=2      # needs command line - will affect legacy scripts. MUST mirror in plotter.py, compare.py
VERBOSE=0
NUMPROCS=1
COMPARE=0

# parse the command line
while [[ $# -gt 0 ]]
do
  OPT="$1"
  shift
  
  case $OPT in
  
    -A|--stdargs)
      STDARGS="$1"
      STDARGS=" "$STDARGS" "
      shift # past argument
    ;;
    
    -J|--maxj)
      MINJ="$1"
      shift # past argument
      MAXJ="$1"
      shift # past argument
    ;;
    
    -L|--maxl)
      MINL="$1"
      shift # past argument
      MAXL="$1"
      shift # past argument
    ;;
    
    -C|--compare)
      COMPARE=1
    ;;
    
    -v|--verbose)
      VERBOSE=1
    ;;
    
    -n|--np)
      NUMPROCS="$1"
      shift # past argument
    ;;
    
    -p|--PDF)
      enscript --color=1 --margins 10c::: -C -f Courier8 -E -M A4 -p ./fouker.ps fouker.py
      ps2pdf ./fouker.ps; rm fouker.ps
      enscript --color=1 --margins 10c::: -C -f Courier8 -E -M A4 -p ./fouvol.ps fouvol.py
      ps2pdf ./fouvol.ps; rm fouvol.ps
      enscript --color=1 --margins 10c::: -C -f Courier8 -E -M A4 -p ./bigrun.ps bigrun.sh
      ps2pdf ./bigrun.ps; rm bigrun.ps
      enscript --color=1 --margins 10c::: -C -f Courier8 -E -M A4 -p ./listrun.ps listrun.sh
      ps2pdf ./listrun.ps; rm listrun.ps
      enscript --color=1 --margins 10c::: -C -f Courier8 -E -M A4 -p ./plotter.ps plotter.py
      ps2pdf ./plotter.ps; rm plotter.ps
      enscript --color=1 --margins 10c::: -C -f Courier8 -E -M A4 -p ./compare.ps compare.py
      ps2pdf ./compare.ps; rm compare.ps
      exit 0
    ;;
    
    -b|--backup)
      DATESTAMP=`date +%Y_%m_%d_%H-%M-%S`
      echo "backing up to $DATESTAMP-{fouvol,fouker,plotter,compare}.py and $DATESTAMP-{bigrun,listrun}.sh"
      cp fouker.py ./backup/$DATESTAMP-fouker.py
      cp fouvol.py ./backup/$DATESTAMP-fouvol.py
      cp plotter.py ./backup/$DATESTAMP-plotter.py
      cp bigrun.sh ./backup/$DATESTAMP-bigrun.sh
      cp listrun.sh ./backup/$DATESTAMP-listrun.sh
      cp compare.py ./backup/$DATESTAMP-compare.py
      exit 0
    ;;
    
    -h|--help)
      echo "Help text: e.g. time ./bigrun.sh -J 8 -v "
      echo "    -p (to create PDF's)"
      echo "    -b (to backup fouvol.py and bigrun.sh to backup directory)"
      echo "    -C set the comparison path as an environment variable in compare.sh"
      echo "    -v (verbose)"
      echo "    -h"
      echo "    -n np or --np np (number of processors for mpirun)"
      echo "    -A \"-v 0 -s 3 -a -0.5 -m 3 --T1 0.1 \""
      echo "    -L 4 8"
      echo "        This -L option corresponds, for example, to ..."
      for((l=0; l<=5; l+=1)); do
        echo "          -L $l => 2^$l = `awk "BEGIN {print int(0.5+2^($l))}"` "
      done
      echo "    -J 1 19"
      echo "        This -J option corresponds, for example, to ..."
      for((j=4; j<=24; j+=2)); do
        echo "          -J $j => 2^$j = `awk "BEGIN {print int(0.5+2^($j))}"` "
      done
      exit 0
    ;;
    
    *)
      echo "option $OPT unknown"
      exit 1 # for echo $?
    ;;
  esac
  
done

# report back if desired
if [ ${VERBOSE} -ne 0 ]; then
  echo MINJ  = "${MINJ}"
  echo MAXJ  = "${MAXJ}"
  echo MINL  = "${MINL}"
  echo MAXL  = "${MAXL}"
  echo   DL  = "${DL}"
fi

#STDARGS="-v 0 -s 3 -a -0.5 -m 1 --T1 0.2 "

# run a sequence of tests.
date +%Y_%m_%d_%H-%M-%S
rm -f output.txt
rm -f results.txt
rm -f errors.txt
declare -A errors
declare -A timings
#for((l=${MINL}; l<=${MAXL}; l++)); do
for((l=${MINL}; l<=${MAXL}; l+=${DL})); do
  for((j=${MINJ}; j<=${MAXJ}; j++)); do
    # build the argument list
    Nt=`awk "BEGIN {print int(0.5+2^($j))}"`
    L=`awk "BEGIN {print int(0.5+2^($l))}"`
    ARGLIST=$STDARGS"-L $L --Nt $Nt "
    # run the solver
    echo; echo `date`
    echo "Running ./fouvol.py $ARGLIST " | tee -a output.txt
    # start a timer
    THEN=$SECONDS
    ./fouvol.py $ARGLIST | tee -a results.txt
    # start a timer
    NOW=$SECONDS
    DURATION=$(( NOW - THEN ))
    err1=`tail -1 results.txt | tr -d ' ' | cut -d"=" -f2`
    errors[$j,$l]=$err1
    timings[$j,$l]=$DURATION
    awk -v Nt="$Nt" -v err1="$err1" -v err2="$err2" -v dur="$DURATION" \
        'BEGIN {printf "%8d & %10.3e & %10.3e & %10.3e & %d\n", Nt, err1, err2, err2/err1, dur}' | tee -a errors.txt
    err2=$err1
  done
  printf "\n"
done
echo 
echo 
echo '   Nt        err1         err2        err2/err1   duration(s) '
cat errors.txt
echo 

# print out the error table (tex and raw)
printf "%% errors for last arglist as: $ARGLIST %%\n"   >   errortable.raw
printf "%% errors for last arglist as: $ARGLIST %%\n" | tee errortable.tex
printf "\\\\begin{tabular}{|l|"    >>    errortable.raw
printf "\\\\begin{tabular}{|l|" | tee -a errortable.tex
#for((l=${MINL}; l<=${MAXL}; l++)); do
for((l=${MINL}; l<=${MAXL}; l+=${DL})); do
  printf "l"    >>    errortable.raw
  printf "l" | tee -a errortable.tex
done
printf "|}\\hline \n"    >>    errortable.raw
printf "|}\\hline \n" | tee -a errortable.tex
for((j=${MINJ}; j<=${MAXJ}; j++)); do
  if [ $j -eq ${MINJ} ]; then
    printf "\$N_t\\downarrow\\ \\\\vert\\ L\\\\rightarrow\$ " $Nt    >>    errortable.raw
    printf "\$N_t\\downarrow\\ \\\\vert\\ L\\\\rightarrow\$ " $Nt | tee -a errortable.tex
#    printf "\$N_t\\ {}_{\\\\rotatebox{90}{\,=}}\\ \\\\vert\ L=\$ " $Nt    >>    errortable.raw
#    printf "\$N_t\\ {}_{\\\\rotatebox{90}{\,=}}\\ \\\\vert\ L=\$ " $Nt | tee -a errortable.tex
#    for((l=${MINL}; l<=${MAXL}; l++)); do
    for((l=${MINL}; l<=${MAXL}; l+=${DL})); do
      # Assume that if we're comparing then there is no Fourier series
      if [ ${COMPARE} -ne 0 ]; then
        L=0
      else
        L=`awk "BEGIN {print int(0.5+2^($l))}"`
      fi
#      L=`awk "BEGIN {print int(0.5+2^($l))}"`
      printf " &    %3d     " $L    >>    errortable.raw
      printf " &   $%3d$    " $L | tee -a errortable.tex
    done
    printf "  \\\\\\\\\\hline\n"    >>    errortable.raw
    printf "  \\\\\\\\\\hline\n" | tee -a errortable.tex
  fi  
  Nt=`awk "BEGIN {print int(0.5+2^($j))}"`
  printf " %-8d  " $Nt    >>    errortable.raw
  printf "$%-8d$ " $Nt | tee -a errortable.tex
#  for((l=${MINL}; l<=${MAXL}; l++)); do
  for((l=${MINL}; l<=${MAXL}; l+=${DL})); do
    # only include the first dollar, sed the other one
    printf " &  %10.3e " ${errors[$j,$l]}                                    >>    errortable.raw
    printf " & $%10.3e " ${errors[$j,$l]} | sed 's/e/(/' | sed 's/$/)\$/' | tee -a errortable.tex
  done
  printf "  \\\\\\\\\n"    >>    errortable.raw
  printf "  \\\\\\\\\n" | tee -a errortable.tex
done
printf "\\hline\\\\end{tabular}%%\n"    >>    errortable.raw
printf "\\hline\\\\end{tabular}%%\n" | tee -a errortable.tex

# print out the timings table (tex and raw)
printf "%% timings for last arglist as: $ARGLIST %%\n"   >   timestable.raw
printf "%% timings for last arglist as: $ARGLIST %%\n" | tee timestable.tex
printf "\\\\begin{tabular}{|l|"    >>    timestable.raw
printf "\\\\begin{tabular}{|l|" | tee -a timestable.tex
#for((l=${MINL}; l<=${MAXL}; l++)); do
for((l=${MINL}; l<=${MAXL}; l+=${DL})); do
  printf "l"    >>    timestable.raw
  printf "l" | tee -a timestable.tex
done
printf "|}\\hline \n"    >>    timestable.raw
printf "|}\\hline \n" | tee -a timestable.tex
#printf "|}\\hline\\\\\\\\\n" | tee -a timestable.tex
for((j=${MINJ}; j<=${MAXJ}; j++)); do
  if [ $j -eq ${MINJ} ]; then
    printf "\$N_t\\downarrow\\ \\\\vert\ L\\\\rightarrow\$ " $Nt    >>    timestable.raw
    printf "\$N_t\\downarrow\\ \\\\vert\ L\\\\rightarrow\$ " $Nt | tee -a timestable.tex
#    printf "\$N_t\\ {}_{\\\\rotatebox{90}{\,=}}\\ \\\\vert\ L=\$ " $Nt    >>    timestable.raw
#    printf "\$N_t\\ {}_{\\\\rotatebox{90}{\,=}}\\ \\\\vert\ L=\$ " $Nt | tee -a timestable.tex
#    for((l=${MINL}; l<=${MAXL}; l++)); do
    for((l=${MINL}; l<=${MAXL}; l+=${DL})); do
      # Assume that if we're comparing then there is no Fourier series
      if [ ${COMPARE} -ne 0 ]; then
        L=0
      else
        L=`awk "BEGIN {print int(0.5+2^($l))}"`
      fi
#      L=`awk "BEGIN {print int(0.5+2^($l))}"`
      printf " &   $%3d$    " $L    >>    timestable.raw
      printf " &    %3d     " $L | tee -a timestable.tex
    done
    printf "  \\\\\\\\\\hline\n"    >>    timestable.raw
    printf "  \\\\\\\\\\hline\n" | tee -a timestable.tex
  fi  
  Nt=`awk "BEGIN {print int(0.5+2^($j))}"`
  printf " %-8d  " $Nt    >>    timestable.raw
  printf "$%-8d$ " $Nt | tee -a timestable.tex
#  for((l=${MINL}; l<=${MAXL}; l++)); do
  for((l=${MINL}; l<=${MAXL}; l+=${DL})); do
    # only include the first dollar, sed the other one
    printf " &  %10.3e " ${timings[$j,$l]}                                     >>    timestable.raw
    printf " & $%10.3e " ${timings[$j,$l]}  | sed 's/e/(/' | sed 's/$/)\$/' | tee -a timestable.tex
  done
  printf "  \\\\\\\\\n"    >>    timestable.raw
  printf "  \\\\\\\\\n" | tee -a timestable.tex
done
printf "\\hline\\\\end{tabular}%%\n"    >>    timestable.raw
printf "\\hline\\\\end{tabular}%%\n" | tee -a timestable.tex

# create the plots
./plotter.py

echo "Using results directory: ./results"
if [ ! -d "./results" ]; then
  mkdir ./results
fi
DIRNAME=results-`echo $STDARGS | tr -d ' '`
echo "Moving results to ./results/$DIRNAME (clearing earlier versions)"
if [ -d "./results/$DIRNAME" ]; then
  rm -rf ./results/$DIRNAME
fi
sleep 5

# update a file called compare.sh to capture output directories for compare.py
if [ ${COMPARE} -ne 0 ]; then
  echo "REFPATH=\"./results/$DIRNAME\"" >> compare.sh
else
  echo "./compare.py -X 2 -c ./results/$DIRNAME -r \$REFPATH" >> compare.sh
fi

mkdir ./results/$DIRNAME
mv errortable.raw errortable.tex errortable.dat ./results/$DIRNAME
mv timestable.raw timestable.tex timestable.dat output.txt ./results/$DIRNAME
#mv errors.eps errors.png errors.jpg ./results/$DIRNAME
mv errors.eps errors.png ./results/$DIRNAME
#mv timings.eps timings.png timings.jpg ./results/$DIRNAME
mv timings.eps timings.png ./results/$DIRNAME
#mv varphirep.eps varphirep.png varphirep.jpg ./results/$DIRNAME
mv varphirep.eps varphirep.png ./results/$DIRNAME
#mv Ut.eps Ut.png Ut.jpg ./results/$DIRNAME
mv Ut.eps Ut.png ./results/$DIRNAME
echo "Contents of ./results"
ls ./results/$DIRNAME
echo "Recommending..."
echo "   mv runout.txt ./results/$DIRNAME"
echo "BE SURE TO MOVE/COPY contents of ./results TO THEIR FINAL DESTINATION"

date +%Y_%m_%d_%H-%M-%S
# fortune | cowsay

exit



