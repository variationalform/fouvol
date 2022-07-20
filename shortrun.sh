#!/bin/bash    

# ./shortrun.sh && ./compare.sh

# BUG: -v 0 is necessary (why?) other values cause problems


read -n 1 -s -r -p 'Removing previous content in ./results - CNTRL-C if not OK ... '
if [ -e ./results ]; then
  echo 'executing rm -rf ./results/*' 
  rm -rf ./results/*
fi

if [ -e ./AAA_finished.txt ]; then
  echo 'executing rm -rf ./results/*' 
  rm -f ./AAA_finished.txt
fi

# initialise runout.txt
date +%Y_%m_%d_%H-%M-%S | tee runout.txt
THEN=$SECONDS
uname -a
chmod u+x shortrun.sh listrun.sh plotter.py bigrun.sh fouvol.py compare.py
rm -f compare.sh

# important: -X for compare will not work on headless server!

# -P is required for gfx for U(t), varphi(t) and the proxy

# JARGS="-J 5 16"
# LARGS="-L 3 11"
JARGS="-J 5 7"
LARGS="-L 3 7"

TARG=" -T 10"
# NOTE: compare wont work without the basic s1 solve
time ./bigrun.sh -A "-v 0 -s 1 -a -0.5       $TARG "        $JARGS -L 0 0 -C | tee -a runout.txt

T1ARG=" --T1 0.5"
echo '#!/bin/bash' > compare_1.sh
# NOTE: compare wont work without the basic s1 solve
# time ./bigrun.sh -A "-v 0 -s 1 -a -0.5       $TARG "        $JARGS -L 0 0 -C | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 1  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 5  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 10 $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 15 $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
chmod u+x compare_1.sh
echo "Running ./compare_1.sh..."
./compare_1.sh
cat ./compare_1.sh

T1ARG=" --T1 0.1"
echo '#!/bin/bash' > compare_2.sh
# NOTE: compare wont work without the basic s1 solve
# time ./bigrun.sh -A "-v 0 -s 1 -a -0.5        $TARG "        $JARGS -L 0 0 -C | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 1   $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 5   $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 10  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 15  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
chmod u+x compare_2.sh
echo "Running ./compare_2.sh..."
./compare_2.sh
cat ./compare_2.sh

T1ARG=" --T1 0.05"
echo '#!/bin/bash' > compare_3.sh
# NOTE: compare wont work without the basic s1 solve
#time ./bigrun.sh -A "-v 0 -s 1 -a -0.5        $TARG "        $JARGS -L 0 0 -C | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 1   $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 5   $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 10  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 15  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
chmod u+x compare_3.sh
echo "Running ./compare_3.sh..."
./compare_3.sh
cat ./compare_3.sh

T1ARG=" --T1 0.01"
echo '#!/bin/bash' > compare_4.sh
# NOTE: compare wont work without the basic s1 solve
#time ./bigrun.sh -A "-v 0 -s 1 -a -0.5        $TARG "        $JARGS -L 0 0 -C | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 1   $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 5   $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 10  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
time ./bigrun.sh -A "-v 0 -s 3 -a -0.5 -m 15  $TARG $T1ARG " $JARGS $LARGS    | tee -a runout.txt
chmod u+x compare_4.sh
echo "Running ./compare_4.sh..."
./compare_4.sh
cat ./compare_4.sh

echo "Finished"
date +%Y_%m_%d_%H-%M-%S | tee -a runout.txt
NOW=$SECONDS
DURATION=$(( NOW - THEN ))
echo "Total Run Time for shortrun (s) = $DURATION" 
echo "Total Run Time for shortrun (m) = `awk "BEGIN {print ($DURATION)/60}"`" 
echo "Total Run Time for shortrun (h) = `awk "BEGIN {print ($DURATION)/3600}"`" 
echo "Total Run Time for shortrun (d) = `awk "BEGIN {print ($DURATION)/3600/12}"`" 

echo "Current directory:"
echo `pwd`
echo "Paste the reference path to the -r options in this file, compare.sh: "
cat compare.sh
chmod u+x ./compare.sh
echo "Execute ./compare.sh "
touch ./AAA_finished.txt
