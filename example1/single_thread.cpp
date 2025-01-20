#include <iostream>
#include <chrono> 

long long sum_of_squares(int start, int end) {
    long long total = 0;
    for (int i = start; i <= end; ++i) {
        total += static_cast<long long>(i) * i;
    }
    return total;
}

int main() {
    int start = 1;
    int end = 10000000;

    auto start_time = std::chrono::high_resolution_clock::now();

    long long result = sum_of_squares(start, end);

    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time = end_time - start_time;

    std::cout << "Sum of squares from " << start << " to " << end << ": " << result << std::endl;
    std::cout << "Time taken: " << elapsed_time.count() << " seconds" << std::endl;

    return 0;
}
