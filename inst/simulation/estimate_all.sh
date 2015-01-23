#!/bin/bash                                                                                                                            

scenarios=(A B C D E F)
scenarios2000=(Acombined Bcombined Ccombined Dcombined Ecombined Fcombined)
scenarios10days=(A10days B10days C10days D10days E10days F10days)
models=(FMPModel SmoothModel-nbinomial-ar1 SmoothModel-nbinomial-matern-ar1)
nodes=(1:50 1:50 1:50)
max_nodes=(60 60 60)
days=(1 1 10)

test="no-test"
if [ "$1" == "test" ]
then
  test="test"
fi

function estimate_single {
  echo $test $1 $2 $3 $4 $5
  if [ "$test" != "test" ]
  then
    ./estimate.sh ~/git/RParallelScreen/ $1 $2 $3 $4 $5
  fi
}

function population_size_scenario {
  _scenarios=$1[@]
  scenarios=("${!_scenarios}")
  _models=$2[@]
  models=("${!_models}")
  nodes=$3
  max_nodes=$4
  numdays=$5

  for model in "${models[@]}"
  do
    for scenario in "${scenarios[@]}"
    do
      #estimate_single $scenario $nodes $max_nodes $model $numdays
      # Count intersections separately:
      estimate_single $scenario $nodes $max_nodes $model
    done
  done
}

population_size_scenario scenarios models "${nodes[0]}" "${max_nodes[0]}" "${days[1]}"
#population_size_scenario scenarios2000 models "${nodes[1]}" "${max_nodes[1]}" "${days[2]}"
#population_size_scenario scenarios10days models "${nodes[2]}" "${max_nodes[2]}" "${days[3]}"

