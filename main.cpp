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
    return 0;
}