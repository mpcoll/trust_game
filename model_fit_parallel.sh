#!/bin/bash

# Author: michel-pierre.coll
# Date: 2020-07-21 08:19:07
# Description: Run multiple model fitting for the Irritability cpp function in
# parrallel. This version creates an output for each participant.
# You can run as many in parrallel as there are threads on your system.
# It takes about 2 hours/participant (or 2 hours for a a parrallel batch).
# RAM usage is about 0.8 Gb/paricipant

# Command:
# Program, data in, data out, number of pairs, pair to start at, pait to end at
#   ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &


Change counter values to run more/different participant
i=0
for i in {11..20}  # Will run part 0-20 in parrallel
do
   ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &
done

 wait  # Waits for the first loop to end

i=0
for i in {21..40} # Will run part 21-40 in parrallel
do
    ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &
done

wait

i=0
for i in {41..60}
do
    ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &
done

wait

i=0
for i in {61..80}
do
    ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &
done

wait

i=0
for i in {81..100}
do
    ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &
done

wait

i=0
for i in {101..120}
do
    ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &
done


wait

i=0
for i in {101..127}
do
     ./code/model_fit/fit_ipomcp /code/data/trust_data.bin /code/model_fit/outputs/trust_params_pair$i.bin 128 $i $((i+ 1)) &
done

wait