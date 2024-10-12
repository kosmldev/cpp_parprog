#include <iostream>
#include <vector>
#include <ctime>

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

int main() {
    vector<int> sizes = {4}; // размеры матриц
    int sizes_size = sizes.size();
    
    for (int one_size = 0; one_size < sizes_size; one_size++) {
        vector<vector<double>> matrix = generate_random_matrix(sizes[one_size]);
        for (int i = 0; i < sizes[one_size]; i++){
            
            for (int j = 0; j < sizes[one_size]; j++){
                cout << matrix[i][j] << "\t";
            }
            cout << endl
        }
    }
    return 0;
}