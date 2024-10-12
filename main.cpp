#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>

extern "C" {
    void dgetrf_(int*, int*, double*, int*, int*, int*);
}

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

// Вызов LU через LAPACK
void LU_decomposition_lapack(vector<double>& matrix, int n) {
    int* ipiv = new int[n + 1];
    int info;
    dgetrf_(&n, &n, matrix.data(), &n, ipiv, &info);
    delete[] ipiv;
}

// Все преобразования вынесли за функции замеры времени, чтобы чисто замерить время LU разложений
// Функция для замера времени LAPACK
double measure_time_lapack(vector<double>& matrix, int n) {
    auto start = chrono::high_resolution_clock::now();
    LU_decomposition_lapack(matrix, n);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    return elapsed.count();
}

// Функция для преобразования матрицы n x n в одномерный массив (LAPACK ест в таком виже)
vector<double> flatten_matrix(const vector<vector<double>>& matrix, int n) {
    vector<double> flat_matrix(n * n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            flat_matrix[i * n + j] = matrix[i][j];
        }
    }

    return flat_matrix;
}


int main() {

    vector<int> sizes = {4, 8, 16, 32, 64, 128, 256, 512, 1024}; // размеры матриц
    int sizes_size = sizes.size();

    for (int one_size = 0; one_size < sizes_size; one_size++) {
        vector<double> matrix = flatten_matrix(generate_random_matrix(sizes[one_size]), sizes[one_size]);
        double time = measure_time_lapack(matrix, sizes[one_size]);
        cout << "Running LAPACK. Matrix size: " << sizes[one_size] << " Time: " << time << endl;
    }
    return 0;
}