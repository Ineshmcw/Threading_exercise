#!/bin/bash

# Runner: Single running


# clear
# FILE_NAME=$1 
# echo "running $FILE_NAME.c file"
# gcc -fopenmp -o $FILE_NAME $FILE_NAME.c
# ./$FILE_NAME


# Runner: Mulitple running and average


# if [ -z "$1" ]; then
#   echo "Usage: $0 <file_name_without_extension>"
#   exit 1
# fi

# FILE_NAME=$1
# TOTAL_TIME=0

# echo "Compiling and running $FILE_NAME.c file"
# gcc -fopenmp -o $FILE_NAME $FILE_NAME.c

# if [ $? -ne 0 ]; then
#   echo "Compilation failed"
#   exit 1
# fi

# for i in {1..5}; do
#   OUTPUT=$(./$FILE_NAME)
# #   echo "Run $i: $OUTPUT"

#   TIME=$(echo "$OUTPUT" | grep -oP '[0-9]+\.[0-9]+(?= seconds)')
#   TOTAL_TIME=$(echo "$TOTAL_TIME + $TIME" | bc)
# done

# AVERAGE_TIME=$(echo "scale=6; $TOTAL_TIME / 5" | bc)
# echo "Average time over 5 runs: $AVERAGE_TIME seconds"



# Runner: single running with math library

FILE_NAME=$1 
echo "Running $FILE_NAME.c file"

gcc -fopenmp -o $FILE_NAME $FILE_NAME.c -lm

if [ $? -eq 0 ]; then
    ./$FILE_NAME
else
    echo "Compilation failed."
fi
