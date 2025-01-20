#include <iostream>
#include <vector>
#include <thread>
#include <chrono>  // For measuring time

// Function to calculate sum of squares for a given range
long long sum_of_squares_range(int start, int end) {
    long long total = 0;
    for (int i = start; i <= end; ++i) {
        total += static_cast<long long>(i) * i;
    }
    return total;
}

// Function to be run by each thread
void thread_task(int start, int end, long long& result) {
    result = sum_of_squares_range(start, end);
}

int main() {
    int start = 1;
    int end = 10000000;
    int num_threads = 4;  // Number of threads
    int range_size = (end - start + 1) / num_threads;  // Split work into equal ranges
    
    std::vector<std::thread> threads(num_threads);
    std::vector<long long> results(num_threads, 0);

    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();

    // Create threads and assign each a range of numbers
    for (int i = 0; i < num_threads; ++i) {
        int range_start = start + i * range_size;
        int range_end = (i == num_threads - 1) ? end : range_start + range_size - 1;  // Ensure the last thread takes the remainder

        threads[i] = std::thread(thread_task, range_start, range_end, std::ref(results[i]));
    }

    // Join threads (wait for them to finish)
    for (auto& t : threads) {
        t.join();
    }

    // Combine the results from all threads
    long long total_sum = 0;
    for (const auto& res : results) {
        total_sum += res;
    }

    // End timer
    auto end_time = std::chrono::high_resolution_clock::now();

    // Calculate elapsed time
    std::chrono::duration<double> elapsed_time = end_time - start_time;

    std::cout << "Sum of squares from " << start << " to " << end << ": " << total_sum << std::endl;
    std::cout << "Time taken: " << elapsed_time.count() << " seconds" << std::endl;

    return 0;
}
