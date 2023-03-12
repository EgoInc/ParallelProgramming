#include <iostream>
#include <omp.h>
#include <random>
#include <iomanip>
#include <fstream>
#include "windows.h"

using namespace std;
//__________________________ÇÀÄÀÍÈÅ 1_______________________
int task_1()
{
#pragma omp parallel num_threads(8)
    {
        int thread_id = omp_get_thread_num();
        int local_number_of_threads = omp_get_num_threads();
        printf("Hello world | thread id: %d | threads: %d\n", thread_id, local_number_of_threads);
    }
    return 0;
}

//__________________________ÇÀÄÀÍÈÅ 2_______________________
void task_2() {
    const int N = 16000;
    int a[N];
    int b[N];

    for (int i = 0; i < N; i++)
    {
        a[i] = i;
    }

    double start = omp_get_wtime();
#pragma omp parallel num_threads(8)
    {
#pragma omp for schedule(static, 8)
        for (int i = 1; i < N - 1; i++)
        {
            b[i] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
        }
    }
    double end = omp_get_wtime();
    printf("Time with static %f\n", end - start);

    start = omp_get_wtime();
#pragma omp parallel num_threads(8)
    {
#pragma omp for schedule(dynamic, 600)
        for (int i = 1; i < N - 1; i++)
        {
            b[i] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
        }
    }
    end = omp_get_wtime();
    printf("Time with dynamic %f\n", end - start);

    start = omp_get_wtime();
#pragma omp parallel num_threads(8)
    {
#pragma omp for schedule(guided, 100)
        for (int i = 1; i < N - 1; i++)
        {
            b[i] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
        }
    }
    end = omp_get_wtime();
    printf("Time with guided %f\n", end - start);

    start = omp_get_wtime();
#pragma omp parallel num_threads(8)
    {
#pragma omp for schedule(runtime)
        for (int i = 1; i < N - 1; i++)
        {
            b[i] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
        }
    }
    end = omp_get_wtime();
    printf("Time with runtime %f\n", end - start);
}

//__________________________ÇÀÄÀÍÈÅ 3. Âåğñèÿ 1_______________________
int task3_1(){
#pragma omp parallel num_threads(8)
    {

        for (int i = 7; i >= 0; i--)
        {
#pragma omp barrier
            {
                if (i == omp_get_thread_num())
                {
#pragma omp critical
                    printf("Hello world | thread id: %d | threads: %d\n", omp_get_thread_num(), 8);
                }
            }
        }
    }
    return 0;
}

//__________________________ÇÀÄÀÍÈÅ 3. Âåğñèÿ 2_______________________
int task3_2() {
#pragma omp parallel num_threads(8)
    {
        Sleep(100 * (8 - omp_get_thread_num()));
        int thread_id = omp_get_thread_num();
        int local_number_of_threads = omp_get_num_threads();
        printf("Hello world | thread id: %d | threads: %d\n", thread_id, local_number_of_threads);
    }
    return 0;
}

//__________________________ÇÀÄÀÍÈÅ 3. Âåğñèÿ 3_______________________
int task3_3() {
    int counter = 7;
#pragma omp parallel num_threads(8) shared(counter)
    {
        int thread_id = omp_get_thread_num();
        while (thread_id != counter) {
            Sleep(1000);
        }
        counter--;

        int local_number_of_threads = omp_get_num_threads();
        printf("Hello world | thread id: %d | threads: %d\n", thread_id, local_number_of_threads);
    }
    return 0;
}

//__________________________ÇÀÄÀÍÈÅ 4_______________________

#define matrixSize 2
double** initMatrix(size_t n) {
    double** m = new double* [n];
    for (int i = 0; i < n; i++) {
        m[i] = new double[n];
        for (int j = 0; j < n; j++) {
                m[i][j] = (double)(rand() % 10);
        }
    }
    return m;
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

void matrixMultiplication(double** A, double** B, int number_of_threads) {
    double time_spent = 0.0;
    clock_t begin = clock();
    double** C;
    C = new double* [matrixSize];
    #pragma omp parallel for num_threads(number_of_threads)
    for (int i = 0; i < matrixSize; i++)
    {
        C[i] = new double[matrixSize];
        for (int j = 0; j < matrixSize; j++)
        {
            C[i][j] = 0;
            for (int g = 0; g < matrixSize; g++) {
                C[i][j] += A[i][g] * B[g][j];
            }
        }
    }
    clock_t end = clock();
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Number of threads: %d \nThe time: %f  seconds\n", number_of_threads, time_spent);
    //printMatrix(C, matrixSize);
    for (int i = 0; i < matrixSize; ++i) {
        delete[] C[i];
    }
    delete[] C;
}

void task_4() {
    double** A = initMatrix(matrixSize);
    //printf("Matrix A: \n");
    //printMatrix(A, matrixSize);
    double** B = initMatrix(matrixSize);
    //printf("Matrix B: \n");
    //printMatrix(B, matrixSize);
    matrixMultiplication(A, B, matrixSize, 2);
    matrixMultiplication(A, B, matrixSize, 4);
    matrixMultiplication(A, B, matrixSize, 8);
}

void main() {
    printf("TASK 1");
    //task_1();
    printf("TASK 2");
    //task_2();
    printf("TASK 3v1");
    //task3_1();
    printf("TASK 3v2");
    //task3_2();
    printf("TASK 3v3");
    //task3_3();
    printf("TASK 4");
    task_4();
}