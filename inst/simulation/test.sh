#!/bin/bash
# ./test.sh ~/git/RParallelScreen/ A

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

python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/simulate.R test "$scenario"
reset
python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R test "$scenario"
reset
python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R test "$scenario"
reset
python "$exec_path"/parallel_r.py -t 1:50 -n 60 -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R test "$scenario"
reset
R --vanilla --args test "$scenario" 1 < ~/git/Winter-Track-Counts/inst/simulation/validate.R
