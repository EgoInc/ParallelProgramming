#include <iostream>
#include <random>
#include <iomanip>
#include <omp.h>
#include <fstream>

using namespace std;

enum MatrixType {
    IDENTITY,
    RANDOM,
    ZERO
};

enum errorCode {
    ZERO_DET
};

double** initMatrix(size_t n, MatrixType type = RANDOM) {
    double** m = new double* [n];
    for (int i = 0; i < n; i++) {
        m[i] = new double[n];
        for (int j = 0; j < n; j++) {
            switch (type) {
            case RANDOM:
                m[i][j] = (double)(rand() % 10);
                break;
            case IDENTITY:
                m[i][j] = (i == j) ? 1 : 0;
                break;
            default:
                m[i][j] = 0;
                break;
            }
        }
    }
    return m;
}

double** copyMatrix(double** m, size_t n) {
    double** m1 = new double* [n];
    for (int i = 0; i < n; i++) {
        m1[i] = new double[n];
        for (int j = 0; j < n; j++) {
            m1[i][j] = m[i][j];
        }
    }
    return m1;
}

void printMatrix(double** m, size_t n) {
    cout << "{\n";
    for (int i = 0; i < n; i++) {
        cout << "    {";
        for (int j = 0; j < n; j++) {
            cout << ((m[i][j] >= 0.0 && m[i][j] < 10.0) ? " " : "") << fixed << setprecision(3) << m[i][j] << (j != n - 1 ? ", " : "");
        }
        cout << "}" << (i != n - 1 ? "," : "") << endl;
    }
    cout << '}';
}

void swapRows(double** m, int row1, int row2, size_t n) {
    double* tmp = m[row1];
    m[row1] = m[row2];
    m[row2] = tmp;
}

double** inverseMatrixOpenMP(double** m, size_t n, int threadsNum = 8) {
    double** I = initMatrix(n, IDENTITY);

    omp_set_num_threads(threadsNum);
    int i, j, k;

#pragma omp parallel for shared(m, I) private(i, j, k)
    for (j = 0; j < n; j++) {
        // swap lines if a_jj == 0
        if (m[j][j] == 0) {
            bool flag = true;
            for (k = j; k < n; k++) {
                if (m[k][j] != 0) {
                    swapRows(m, j, k, n);
                    flag = false;
                    break;
                }
            }
            if (flag) {
                throw ZERO_DET;
            }
        }
        // делим строку на a_jj
        double div = m[j][j];
        for (int i = 0; i < n; i++) {
            m[j][i] /= div;
            I[j][i] /= div;
        }
        // вот здесь надо будет распараллелить: из всех строк, кроме этой, вычитаем эту с коэф-ом
        for (int i = 0; i < n; i++) {
            if (i == j)
                continue;
            double c = m[i][j];

            for (int k = 0; k < n; k++) {
                m[i][k] -= c * m[j][k];
                I[i][k] -= c * I[j][k];
            }
        }
    }
    return I;
}


void InverseMatrix(double** inputMatrix, int matrix_size) {
    
    for (int threads = 2; threads <= 8; threads = threads * 2) {
        double** T = copyMatrix(inputMatrix, matrix_size);

        double timerOpenMp;

        double** N;



        cout << "Matrix size: " << matrix_size << ", threads: " << threads << endl;

        try {
            timerOpenMp = omp_get_wtime();
            N = inverseMatrixOpenMP(T, matrix_size, threads);
            timerOpenMp = omp_get_wtime() - timerOpenMp;
            if (matrix_size <= 15) {
                cout << "\n\nResult:\n\n";
                printMatrix(N, matrix_size);
            }
            cout << "OpenMP execution time: " << timerOpenMp << endl;
        }
        catch (errorCode err) {
            if (err == ZERO_DET) {
                cout << "Inverse matrix doesn't exist\n";
            }
        }
    }

}

double** createMatrix(int matrix_size) {
    double** InputMatrix = initMatrix(matrix_size, RANDOM);
    if (matrix_size <= 15) {
        cout << "\nInput matrix (generated randomly):\n\n";
        printMatrix(InputMatrix, matrix_size);
    }
    return InputMatrix;
}
int main(int argc) {

    cout << "Task 10\n";
    cout << "Construction of the inverse matrix\n";
    
    cout << "\n________________________________Test 1_________________________________\n";
    double** InputMatrix = createMatrix(3);
    cout << "Let's start with 5x5 matrix\n";
        InverseMatrix(InputMatrix, 3);
    
    cout << "\n_____________________Test 3_________________________\n";
    cout << "Matrices of size larger than 15 are not output to the table, but let's check the 500x500 matrix\n";
    cout << "=======500x500 \n";
     InputMatrix = createMatrix(500);
    InverseMatrix(InputMatrix, 500);
   

  
    cout << "\n_____________________Test 4_________________________\n";
    cout << "=======1000x1000 \n";
    InputMatrix = createMatrix(1000);
    InverseMatrix(InputMatrix, 1000);
    

    cout << "\n_____________________Test 5_________________________\n";
    InputMatrix = createMatrix(2000);
    InverseMatrix(InputMatrix, 2000);
    

    cout << "\n_____________________Test 5_________________________\n";
    InputMatrix = createMatrix(4000);
    InverseMatrix(InputMatrix, 4000);
    
    return 0;
}