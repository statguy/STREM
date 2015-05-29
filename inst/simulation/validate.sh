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
if [ -z "$5" ]
then
  echo "Missing parameter 'model'".
  exit
fi
if [ -z "$6" ]
then
  echo "Missing parameter 'free_mem'".
  exit
fi

exec_path=$1
scenario=$2
iterations=$3
max_nodes=$4
model=$5
free_mem=$6

if [ ! "$exec_path"/git_uptodate.sh ]
then
  exit
fi

if [ "$model" != "FMPModel" ]
then
  python "$exec_path"/parallel_r.py -p 15 -t "$iterations" -n "$max_nodes" -l 20.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/validate.R notest "$scenario" "$model"
else
  python "$exec_path"/parallel_r.py -p 15 -t "$iterations" -n "$max_nodes" -l 20.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/validate_boot_ci.R notest "$scenario" "$model"
fi

# ./validate.sh ~/git/RParallelScreen/ A 1 1 FMPModel 0
# ./validate.sh ~/git/RParallelScreen/ A 1:50 60 FMPModel 0
