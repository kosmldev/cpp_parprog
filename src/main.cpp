#include <iostream>
#include <vector>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <cmath>
#include <string>
#include <cblas.h>

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
    int* ipiv = new int[n];
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


void print_table(const bool& run_only_lapack, vector<int>& vec1, const vector<double>& vec2, const vector<double>& vec3, 
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
        if ( run_only_lapack ){
            cout << setw(width) << vec1[i]
                    << setw(width) << vec3[i]
                    << setw(width) << vec5[i]
                    << setw(width) << vec7[i]
                    << endl;
        } else {
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
}

int main() {
    /*
        На вход = вектор из размеров матрицы n x n, которая генерируется случайными числами от -1 до 1
        Каждый размер матрицы провегоняется meas_cnt раз, замеряется min, mean, max времени работы
        Сначала прогоняется обычное нативное LU разложение через цикл
        Далее прогоняется через LAPACK
    */
    
    bool run_only_lapack = stoi(getenv("RUN_ONLY_LAPACK"));
    if ( run_only_lapack ) {
        cout << "Running only LAPACK" << endl;
    }

    int num_treads_for_lapack = stoi(getenv("SET_THREADS"));
    openblas_set_num_threads(num_treads_for_lapack);
    int num_threads = openblas_get_num_threads();
    cout << "OpenBLAS uses " << num_threads << " threads" << endl;

    int sizes_length = stoi(getenv("SET_SIZES")); // размеры матриц
    cout << "Sizes of matrix = " << sizes_length << endl;
    vector<int> sizes(sizes_length); // размеры матриц
    for (int one_size = 0; one_size < sizes_length; one_size++) {
        sizes[one_size] = 4 * pow(2, one_size);
    }
    int sizes_size = sizes.size();
    int meas_cnt = stoi(getenv("SET_ITERATIONS")); // сколько измерений за итерацию делать
    cout << "Iteration per matrix size = " << meas_cnt << endl;
    vector<double> output_min_times(sizes_size);
    vector<double> output_mean_times(sizes_size);
    vector<double> output_max_times(sizes_size);
    vector<double> output_min_times_lpck(sizes_size);
    vector<double> output_mean_times_lpck(sizes_size);
    vector<double> output_max_times_lpck(sizes_size);
    vector<string> headers;

    if ( !run_only_lapack ) {
        headers = {"Size","Native_min", "LAPACK_min", "Native_mean", "LAPACK_mean", "Native_max", "LAPACK_max"};
    } else {
        headers = {"Size", "LAPACK_min", "LAPACK_mean", "LAPACK_max"};
    }

    if ( !run_only_lapack ) {
        
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
        run_only_lapack,
        sizes, output_min_times, output_min_times_lpck, 
        output_mean_times, output_mean_times_lpck, output_max_times,
        output_max_times_lpck, headers
    );

    return 0;
}
