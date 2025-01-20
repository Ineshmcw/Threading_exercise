#!/bin/bash

LOG_FILE="benchmark_log.txt"
NUM_ITERATIONS=5  # Number of iterations to run for more stable results
#!/bin/bash

folders=("output_files_modified" "output_files_original" 
         "./tinyxml2_modified/build" "./tinyxml2_original/build")

# Check and create folders if they don't exist
for folder in "${folders[@]}"; do
  if [ ! -d "$folder" ]; then
    mkdir -p "$folder"
    echo "Created $folder folder"
  fi
done

# Build the projects
for project in "tinyxml2_modified" "tinyxml2_original"; do
  cd "./$project/build"
  cmake ..
  make
  cd ../../
done

echo "Benchmark Log" > $LOG_FILE
echo "======================" >> $LOG_FILE
echo "Start Time: $(date)" >> $LOG_FILE
echo "======================" >> $LOG_FILE

# # Clear system cache before starting the benchmarks
# echo "Clearing system cache..."
# sudo sync; sudo echo 3 > /proc/sys/vm/drop_caches
# sleep 2  # Allow system to stabilize after clearing cache

# Set CPU governor to performance (optional, helps if CPU is throttling)
# echo "Setting CPU governor to performance"
# sudo cpupower frequency-set --governor performance

# Compile C++ programs (before changes and after changes)
echo "Compiling C++ programs..." >> $LOG_FILE
g++ -std=c++17 -I./tinyxml2_original -L./tinyxml2_original/build ./original_multithreading_test.cpp -pthread -l tinyxml2 -o ./original_multithreading_testexecutable >> $LOG_FILE 2>&1
g++ -std=c++17 -I./tinyxml2_modified -L./tinyxml2_modified/build ./modified_multithreading_test.cpp -pthread -l:libtinyxml2.a -o ./modified_multithreading_test_executable >> $LOG_FILE 2>&1

# Run the benchmark for the "before changes" version multiple times
echo -e "\nRunning multi-threaded program (before changes)" >> $LOG_FILE
BEFORE_TIME_TOTAL=0
for (( i=1; i<=NUM_ITERATIONS; i++ ))
do
    echo "Iteration $i..." >> $LOG_FILE
    sleep 2  # Allow system to stabilize after clearing cache
    START_TIME=$(date +%s.%N)
    ./original_multithreading_testexecutable >> $LOG_FILE
    END_TIME=$(date +%s.%N)
    ITERATION_TIME=$(echo "$END_TIME - $START_TIME" | bc)
    BEFORE_TIME_TOTAL=$(echo "$BEFORE_TIME_TOTAL + $ITERATION_TIME" | bc)
done

# Calculate average time for "before changes"
BEFORE_TIME_AVG=$(echo "$BEFORE_TIME_TOTAL / $NUM_ITERATIONS" | bc -l)dgrfg 
echo "Average execution time (before changes): $BEFORE_TIME_AVG seconds" >> $LOG_FILE

# Run the benchmark for the "after changes" version multiple times
echo -e "\nRunning multi-threaded program (after changes)" >> $LOG_FILE
AFTER_TIME_TOTAL=0
for (( i=1; i<=NUM_ITERATIONS; i++ ))
do
    echo "Iteration $i..." >> $LOG_FILE
    START_TIME=$(date +%s.%N)
    ./modified_multithreading_test_executable >> $LOG_FILE
    END_TIME=$(date +%s.%N)
    ITERATION_TIME=$(echo "$END_TIME - $START_TIME" | bc)
    AFTER_TIME_TOTAL=$(echo "$AFTER_TIME_TOTAL + $ITERATION_TIME" | bc)
done

# Calculate average time for "after changes"
AFTER_TIME_AVG=$(echo "$AFTER_TIME_TOTAL / $NUM_ITERATIONS" | bc -l)
echo "Average execution time (after changes): $AFTER_TIME_AVG seconds" >> $LOG_FILE

# Calculate the speed improvement
if (( $(echo "$BEFORE_TIME_AVG > 0" | bc -l) )); then
    SPEED_IMPROVEMENT=$(echo "scale=2; 100 * ($BEFORE_TIME_AVG - $AFTER_TIME_AVG) / $BEFORE_TIME_AVG" | bc)
    echo "Speed improvement: $SPEED_IMPROVEMENT%" >> $LOG_FILE
else
    echo "Error: Before execution time is zero or invalid." >> $LOG_FILE
fi

# Note the current time of execution
echo -e "\nEnd Time: $(date)" >> $LOG_FILE
echo "======================" >> $LOG_FILE

# Print the log to the console for review
cat $LOG_FILE
