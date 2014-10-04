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

exec_path=$1
scenario=$2
iterations=$3
max_nodes=$4
model=$5

python "$exec_path"/parallel_r.py -t "$iterations" -n "$max_nodes" -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/population_size.R notest "$scenario" "$model"

if [ "$model" != "FMPModel" ]
then
  python "$exec_path"/parallel_r.py -t "$iterations" -n "$max_nodes" -l 10.0 -b ~/tmp/blacklist.txt -v ~/git/Winter-Track-Counts/inst/simulation/validate.R notest "$scenario" "$model"
fi

# ./population_size.sh ~/git/RParallelScreen/ A 1:50 60 FMPModel 1
# ./population_size.sh ~/git/RParallelScreen/ A 1:50 60 SmoothModel-nbinomial-ar1
# ./population_size.sh ~/git/RParallelScreen/ A 1:50 60 SmoothModel-nbinomial-matern-ar1

# ./population_size.sh ~/git/RParallelScreen/ Acombined 1:10 11 FMPModel 1
# ./population_size.sh ~/git/RParallelScreen/ Acombined 1:10 11 SmoothModel-nbinomial-ar1
# ./population_size.sh ~/git/RParallelScreen/ Acombined 1:10 11 SmoothModel-nbinomial-matern-ar1

# ./population_size.sh ~/git/RParallelScreen/ A10days 1:50 60 FMPModel 10
# ./population_size.sh ~/git/RParallelScreen/ A10days 1:50 60 SmoothModel-nbinomial-ar1
# ./population_size.sh ~/git/RParallelScreen/ A10days 1:50 60 SmoothModel-nbinomial-matern-ar1
