#!/bin/bash

LOG_FILE="benchmark_log.txt"

echo "Benchmark Log" > $LOG_FILE
echo "======================" >> $LOG_FILE
echo "Start Time: $(date)" >> $LOG_FILE
echo "======================" >> $LOG_FILE

echo "Compiling C++ programs..." >> $LOG_FILE
g++ -std=c++11 -o single_thread single_thread.cpp >> $LOG_FILE 2>&1
g++ -std=c++11 -o multi_thread multi_thread.cpp >> $LOG_FILE 2>&1

echo -e "\nRunning single-threaded program..." >> $LOG_FILE
START_TIME=$(date +%s.%N)
./single_thread >> $LOG_FILE
END_TIME=$(date +%s.%N)
SINGLE_THREAD_TIME=$(echo "$END_TIME - $START_TIME" | bc)
echo "Single-thread execution time: $SINGLE_THREAD_TIME seconds" >> $LOG_FILE

echo -e "\nRunning multi-threaded program..." >> $LOG_FILE
START_TIME=$(date +%s.%N)
./multi_thread >> $LOG_FILE
END_TIME=$(date +%s.%N)
MULTI_THREAD_TIME=$(echo "$END_TIME - $START_TIME" | bc)
echo "Multi-thread execution time: $MULTI_THREAD_TIME seconds" >> $LOG_FILE

# Calculate the speed improvement in percentage
if (( $(echo "$SINGLE_THREAD_TIME > 0" | bc -l) )); then
    SPEED_IMPROVEMENT=$(echo "scale=2; 100 * ($SINGLE_THREAD_TIME - $MULTI_THREAD_TIME) / $SINGLE_THREAD_TIME" | bc)
    echo "Speed improvement: $SPEED_IMPROVEMENT%" >> $LOG_FILE
else
    echo "Error: Single-thread execution time is zero or invalid." >> $LOG_FILE
fi

# Note the current time of execution
echo -e "\nEnd Time: $(date)" >> $LOG_FILE
echo "======================" >> $LOG_FILE

# Print the log to the console for review
cat $LOG_FILE
