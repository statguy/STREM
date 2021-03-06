#!/bin/bash                                                                                                                            

scenarios=(A B C D E F)
scenarios2000=(Acombined Bcombined Ccombined Dcombined Ecombined Fcombined)
scenarios10days=(A10days B10days C10days D10days E10days F10days)
models=(FMPModel SmoothModelMean-nbinomial-ar1 SmoothModel-nbinomial-ar1 SmoothModel-nbinomial-matern-ar1 SmoothModelMean-nbinomial-ar1-priors1)
nodes=(1:50 1:50 1:50)
max_nodes=(60 60 60)
#free_mem=(10000 10000 10000 10000 28000 28000)
free_mem=(10000 10000 10000 10000 10000 10000) # Assuming habitat uses have been determined before

test="no-test"
if [ "$1" == "test" ]
then
  test="test"
fi

function validate_single {
  echo $test $1 $2 $3 $4 $5
  if [ "$test" != "test" ]
  then
    ./validate.sh ~/git/RParallelScreen/ $1 $2 $3 $4 $5
  fi
}

function population_size_scenario {
  _scenarios=$1[@]
  scenarios=(${!_scenarios})
  _models=$2[@]
  models=(${!_models})
  nodes=$3
  max_nodes=$4
  _free_mem=$5[@]
  free_mem=(${!_free_mem})

  for model in ${models[@]}
  do
    for i in ${!scenarios[@]}
    do
      echo validate_single ${scenarios[$i]} $nodes $max_nodes $model ${free_mem[$i]}
      validate_single ${scenarios[$i]} $nodes $max_nodes $model ${free_mem[$i]}
    done
  done
}

population_size_scenario scenarios models ${nodes[0]} ${max_nodes[0]} free_mem
population_size_scenario scenarios2000 models ${nodes[1]} ${max_nodes[1]} free_mem
population_size_scenario scenarios10days models ${nodes[2]} ${max_nodes[2]} free_mem
