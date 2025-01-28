#include <iostream>
#include <vector>
#include <chrono>
using namespace std;

void mulMat(vector<vector<int>>& m1, vector<vector<int>>& m2, 
            vector<vector<int>>& res) {
    int r1 = m1.size();
    int c1 = m1[0].size();
    int r2 = m2.size();
    int c2 = m2[0].size();

    if (c1 != r2) {
        cout << "Invalid Input" << endl;
        exit(EXIT_FAILURE);
    }

    // Resize result matrix to fit the result dimensions
    res.resize(r1, vector<int>(c2, 0)); 
  
    for (int i = 0; i < r1; i++) {
        for (int j = 0; j < c2; j++) {
            for (int k = 0; k < c1; k++) {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

// Driver code
int main() {
    int N = 64;
    vector<vector<int>> m1(N, vector<int>(N, 1));
    vector<vector<int>> m2(N, vector<int>(N, 1));
    vector<vector<int>> res;

    auto start = chrono::high_resolution_clock::now();
    mulMat(m1, m2, res);
    auto end = chrono::high_resolution_clock::now();

    cout << "Multiplication of given two matrices (first 5x5) is:\n";
    for (int i = 0; i < 5 && i < res.size(); i++) {
        for (int j = 0; j < 5 && j < res[i].size(); j++) {
            cout << res[i][j] << "\t";
        }
        cout << endl;
    }

    chrono::duration<double> execTime = end - start;
    cout << "Execution time: " << execTime.count() << " seconds" << endl;

    return 0;
}

