#include <iostream>
#include <omp.h>

#include "windows.h"

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

void main() {
    printf("TASK 1");
    task_1();
    printf("TASK 3v1");
    task3_1();
    printf("TASK 3v2");
    task3_2();
    printf("TASK 3v3");
    task3_3();
}