#!/bin/bash

./population_size.sh ~/git/RParallelScreen/ E 1:50 60 FMPModel 1
./population_size.sh ~/git/RParallelScreen/ E 1:50 60 SmoothModel-nbinomial-ar1
./population_size.sh ~/git/RParallelScreen/ E 1:50 60 SmoothModel-nbinomial-matern-ar1
./population_size.sh ~/git/RParallelScreen/ F 1:50 60 FMPModel 1
./population_size.sh ~/git/RParallelScreen/ F 1:50 60 SmoothModel-nbinomial-ar1
./population_size.sh ~/git/RParallelScreen/ F 1:50 60 SmoothModel-nbinomial-matern-ar1

./population_size.sh ~/git/RParallelScreen/ Ecombined 1:10 11 FMPModel 1
./population_size.sh ~/git/RParallelScreen/ Ecombined 1:10 11 SmoothModel-nbinomial-ar1
./population_size.sh ~/git/RParallelScreen/ Ecombined 1:10 11 SmoothModel-nbinomial-matern-ar1
./population_size.sh ~/git/RParallelScreen/ Fcombined 1:10 11 FMPModel 1
./population_size.sh ~/git/RParallelScreen/ Fcombined 1:10 11 SmoothModel-nbinomial-ar1
./population_size.sh ~/git/RParallelScreen/ Fcombined 1:10 11 SmoothModel-nbinomial-matern-ar1

./population_size.sh ~/git/RParallelScreen/ E10days 1:50 60 FMPModel 10
./population_size.sh ~/git/RParallelScreen/ E10days 1:50 60 SmoothModel-nbinomial-ar1
./population_size.sh ~/git/RParallelScreen/ E10days 1:50 60 SmoothModel-nbinomial-matern-ar1
./population_size.sh ~/git/RParallelScreen/ F10days 1:50 60 FMPModel 10
./population_size.sh ~/git/RParallelScreen/ F10days 1:50 60 SmoothModel-nbinomial-ar1
./population_size.sh ~/git/RParallelScreen/ F10days 1:50 60 SmoothModel-nbinomial-matern-ar1
