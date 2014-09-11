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
if [ -z "$3" ]
then
  echo "Missing parameter 'iterations'".
  exit
fi
if [ -z "$4" ]
then
  echo "Missing parameter 'max_nodes'".
  exit
fi

exec_path=$1
scenario=$2
iterations=$3
max_nodes=$4
count_days=1

python "$exec_path"/parallel_r.py -t "$iterations" -n "$max_nodes" -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/count_intersections.R notest "$scenario" "$count_days"
reset
python "$exec_path"/parallel_r.py -t "$iterations" -n "$max_nodes" -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/estimate.R notest "$scenario"
reset
python "$exec_path"/parallel_r.py -t "$iterations" -n "$max_nodes" -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest "$scenario"
reset
python "$exec_path"/parallel_r.py -t "$iterations" -n "$max_nodes" -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/validate.R notest "$scenario"
reset

# ./population_size.sh ~/git/RParallelScreen/ combinedA 1:9 10
