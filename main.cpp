#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>

using namespace std;

// Функция для генерации случайной матрицы размером n на n
vector<vector<double>> generate_random_matrix(int n) {
    vector<vector<double>> matrix(n, vector<double>(n));
    
    // Инициализация генератора случайных чисел
    srand(time(0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 2.0 * (static_cast<double>(rand()) / RAND_MAX) - 1.0;
        }
    }

    return matrix;
}

// Функция для LU-разложения
void LU_decomposition(vector<vector<double>>& matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            matrix[j][i] /= matrix[i][i];
            for (int k = i + 1; k < n; k++) {
                matrix[j][k] -= matrix[j][i] * matrix[i][k];
            }
        }
    }
}

// Функция для замера времени
double measure_time(vector<vector<double>>& matrix, int n) {
    auto start = chrono::high_resolution_clock::now();
    LU_decomposition(matrix, n);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    return elapsed.count();
}

int main() {

    vector<int> sizes = {4, 8, 16, 32, 64, 128, 256, 512, 1024}; // размеры матриц
    int sizes_size = sizes.size();

    for (int one_size = 0; one_size < sizes_size; one_size++) {
        vector<vector<double>> matrix = generate_random_matrix(sizes[one_size]);
        double time = measure_time(matrix, sizes[one_size]);
        cout << "Running native. Matrix size: " << sizes[one_size] << " Time: " << time << endl;
    }
    return 0;
}