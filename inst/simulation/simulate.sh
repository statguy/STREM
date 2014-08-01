#!/bin/bash

if [ -z "$1" ]
then
  echo "Missing parameter 'exec path'".
  exit
fi

if [ -z "$2" ]
then
  echo "Missing parameter 'scenario'".
  exit
fi

exec_path=$1
scenario=$2
count_days=59

python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R notest "$scenario"
reset
python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest "$scenario" "$count_days"
reset
python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest "$scenario"
reset
python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest "$scenario"
reset
R --vanilla --args notest "$scenario" 1 < ~/git/Winter-Track-Counts/inst/simulation/validate.R
