#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iomanip>

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


void print_table(const vector<int>& vec1, const vector<double>& vec2, const vector<double>& vec3, 
                const vector<double>& vec4, const vector<double>& vec5, const vector<double>& vec6,
                const vector<double>& vec7, const vector<string>& headers) {
    const int width = 15;

    for (const auto& header : headers) {
        cout << setw(width) << header;
    }
    cout << endl;

    cout << scientific << setprecision(4);

    size_t length = vec1.size();

    for (size_t i = 0; i < length; ++i) {
        cout << setw(width) << vec1[i]
                  << setw(width) << vec2[i]
                  << setw(width) << vec3[i]
                  << setw(width) << vec4[i]
                  << setw(width) << vec5[i]
                  << setw(width) << vec6[i]
                  << setw(width) << vec7[i]
                  << endl;
    }
}

int main() {
    /*
        На вход = вектор из размеров матрицы n x n, которая генерируется случайными числами от -1 до 1
        Каждый размер матрицы провегоняется meas_cnt раз, замеряется min, mean, max времени работы
        Сначала прогоняется обычное нативное LU разложение через цикл
        Далее прогоняется через LAPACK
    */
    
    vector<int> sizes = {4, 8, 16, 32, 64, 128, 256, 512, 1024}; // размеры матриц
    int sizes_size = sizes.size();
    int meas_cnt = 10; // сколько измерений за итерацию делать
    vector<double> output_min_times(sizes_size);
    vector<double> output_mean_times(sizes_size);
    vector<double> output_max_times(sizes_size);
    vector<double> output_min_times_lpck(sizes_size);
    vector<double> output_mean_times_lpck(sizes_size);
    vector<double> output_max_times_lpck(sizes_size);
    vector<string> headers = {"Size","Native_min", "LAPACK_min", "Native_mean", "LAPACK_mean", "Native_max", "LAPACK_max"};

    for (int one_size = 0; one_size < sizes_size; one_size++) {
        vector<double> times_vector(meas_cnt);
        for (int iteration = 0; iteration < meas_cnt; iteration++)  {
            cout << "Running native. Matrix size: " << sizes[one_size] << " Iteration: " << iteration + 1 << endl;
            vector<vector<double>> matrix = generate_random_matrix(sizes[one_size]);
            double time = measure_time(matrix, sizes[one_size]);
            times_vector[iteration] = time;
        }
        output_min_times[one_size] = *min_element(times_vector.begin(), times_vector.end());
        output_max_times[one_size] = *max_element(times_vector.begin(), times_vector.end());
        output_mean_times[one_size] = (accumulate(times_vector.begin(), times_vector.end(), 0.0)) / static_cast<double>(times_vector.size());
    }

    for (int one_size = 0; one_size < sizes_size; one_size++) {
        vector<double> times_vector(meas_cnt);
        for (int iteration = 0; iteration < meas_cnt; iteration++)  {
            cout << "Running LAPACK. Matrix size: " << sizes[one_size] << " Iteration: " << iteration + 1 << endl;
            vector<double> matrix = flatten_matrix(generate_random_matrix(sizes[one_size]), sizes[one_size]);
            double time = measure_time_lapack(matrix, sizes[one_size]);
            times_vector[iteration] = time;
        }
        output_min_times_lpck[one_size] = *min_element(times_vector.begin(), times_vector.end());
        output_max_times_lpck[one_size] = *max_element(times_vector.begin(), times_vector.end());
        output_mean_times_lpck[one_size] =  (accumulate(times_vector.begin(), times_vector.end(), 0.0)) / static_cast<double>(times_vector.size());
    }

    print_table(
        sizes, output_min_times, output_min_times_lpck, 
        output_mean_times, output_mean_times_lpck, output_max_times,
        output_max_times_lpck, headers
    );

    return 0;
}
