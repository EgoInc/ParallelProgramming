#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <conio.h>
#include <time.h>
#include <omp.h>
#include <windows.h>
#include <vector>
#include <math.h>
#include <iomanip>
#include <cstdlib>
/*
using namespace std;

void printa(int n,int a[][100])            //����� ������� "a" � �������
{
    printf("������ a:\n");
    for (int i = 0; i < n; i++)
    {
        printf("a[%d][] - ", i);
        for (int j = 0; j < n; j++)
            printf("%d  ", a[i][j]);
        printf("\n");
    }
}

void Task4MatrixInverse() //����� �����
{
    cout << endl << endl << "Solution for SLAU Gauss" << endl;
    int i = 0, j = 0, k = 0;
    int n;
    int temp;
    double det = 1;

    const double EPS = 1E-9;
    int timein, timeout, timeres = 0;
    cout << "Enter n: ";
    cin >> n;

    //n=2;
    srand(1);
    vector <vector<double>> Matrix(n, vector<double>(n));
    cout << "original matrix" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Matrix[i][j] = rand() % 10 + 1;
            cout<< Matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << "starts to reverse" << endl;
    //Matrix[0][0] = 1;// = -1
    //Matrix[0][1] = 1;
    //Matrix[1][0] = 1;
    //Matrix[1][1] = 0;

    vector <vector<double>> E(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    timein = GetTickCount();

    // �������� �� ������� ������� (������ �� ���������)
// ������ ����. �� ������ ����� ���������� ������ ���
// � �������� ������� ������������ � ������� �����������
    for (int k = 0; k < n; ++k)
    {
        // ���� ������� �� ������� ��������� � ��������
        // ������ - ����, �� ���� ������, ��� �������
        // ���� �� ������� �� �������, � ������ ������
        // �������
        if (fabs(Matrix[k][k]) < 1e-8)
        {
            // ����, ��������� � ���, ��� ��� ��������� ����� �����
            bool changed = false;

            // ��� �� �������, ������������� ���� ��������
            for (int i = k + 1; i < n; ++i)
            {
                // ���� ����� ������, ��� � ��� �� �������
                // ������� ��������� �������
                if (fabs(Matrix[i][k]) > 1e-8)
                {
                    // ������ ��������� � �������� ������ �������
                    // ��� � �������� �������, ��� � � ���������
                    swap(Matrix[k], Matrix[i]);
                    swap(E[k], E[i]);

                    // ������� ���� - �������� � ������������ ������ �����
                    changed = true;

                    break;
                }
            }

            // ���� ����� ����� ��������� �� ��� - ������� �� ����� ����
            // ��������
            if (!changed)
            {
                cout << "Matrix cannot inverse!!!" << endl;
                // �������� � ������� ���������
            }
        }
        // ���������� �������� - ������������ �������
        double div = Matrix[k][k];
        // ��� �������� �������� ������ ����� �� ������������
        // ������� ��� � �������� �������, ��� � � ���������
        for (int j = 0; j < n; ++j)
        {
            Matrix[k][j] /= div;
            E[k][j] /= div;
        }
        // ��� �� �������, ������� ����������� ���� ��������
        for (int i = k + 1; i < n; ++i)
        {
            // ���������� ��������� - ������� ��������� ������,
            // ������������� ��� ������������ ��������� ��������
            // ������
            double multi = Matrix[i][k];

            // �������� �� ��������� ������ ��������, ����������
            // �� ���������� ����� ��������� ��� � ��������,
            // ��� � � ��������� �������
            for (int j = 0; j < n; ++j)
            {
                Matrix[i][j] -= multi * Matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }

    // �������� �� �������� ����������� �������, ����������
    // �� ������ ����, ����� �����
    // �� ������ ����� ���������� �������� ���, � �� ��������
    // ������� ������������ ����������� ���������, � �� ��������� -
    // ��������
    for (int k = n - 1; k > 0; --k)
    {
        // ��� �� �������, ������� ����������� ���� ��������
        for (int i = k - 1; i + 1 > 0; --i)
        {
            // ���������� ��������� - ������� ��������� ������,
            // ������������� ��� ������������ ��������� ��������
            // ������
            double multi = Matrix[i][k];

            // �������� �� ��������� ������ ��������, ����������
            // �� ���������� ����� ��������� ��� � ��������,
            // ��� � � ��������� �������
            for (int j = 0; j < n; ++j)
            {
                Matrix[i][j] -= multi * Matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }

    /*for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
            cout<<E[i][j];
        cout<<endl;
    }*/
    cout << "inverse matrix" << endl;;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << E[i][j] << " ";
        }
        cout << endl;
    }
    timeout = GetTickCount();
    timeres = timeout - timein;
    cout << "Our time: " << timeres << endl;
}

void Task4MatrixInverseParallel()
{
    int i = 0, j = 0, k = 0;
    int n;
    int temp;
    double det = 1;

    const double EPS = 1E-9;
    int timein, timeout, timeres = 0;
    cout << "Enter n: ";
    cin >> n;

    //n=2;
    srand(1);
    cout << "original matrix" << endl;
    vector <vector<double>> Matrix(n, vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            Matrix[i][j] = rand() % 10 +1;
            cout << Matrix[i][j] << " ";
        }
        cout << endl;
    }
    //Matrix[0][0] = 1;// = -1
    //Matrix[0][1] = 1;
    //Matrix[1][0] = 1;
    //Matrix[1][1] = 0;

    vector <vector<double>> E(n, vector<double>(n));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    timein = GetTickCount();

    // �������� �� ������� ������� (������ �� ���������)
// ������ ����. �� ������ ����� ���������� ������ ���
// � �������� ������� ������������ � ������� �����������
    for (int k = 0; k < n; ++k)
    {
        // ���� ������� �� ������� ��������� � ��������
        // ������ - ����, �� ���� ������, ��� �������
        // ���� �� ������� �� �������, � ������ ������
        // �������
        if (fabs(Matrix[k][k]) < 1e-8)
        {
            // ����, ��������� � ���, ��� ��� ��������� ����� �����
            bool changed = false;

            // ��� �� �������, ������������� ���� ��������
            for (int i = k + 1; i < n; ++i)
            {
                // ���� ����� ������, ��� � ��� �� �������
                // ������� ��������� �������
                if (fabs(Matrix[i][k]) > 1e-8)
                {
                    // ������ ��������� � �������� ������ �������
                    // ��� � �������� �������, ��� � � ���������
                    swap(Matrix[k], Matrix[i]);
                    swap(E[k], E[i]);

                    // ������� ���� - �������� � ������������ ������ �����
                    changed = true;

                    break;
                }
            }

            // ���� ����� ����� ��������� �� ��� - ������� �� ����� ����
            // ��������
            if (!changed)
            {
                cout << "Matrix cannot inverse!!!" << endl;
                // �������� � ������� ���������
            }
        }

        // ���������� �������� - ������������ �������
        double div = Matrix[k][k];

        // ��� �������� �������� ������ ����� �� ������������
        // ������� ��� � �������� �������, ��� � � ���������
        for (int j = 0; j < n; ++j)
        {
            Matrix[k][j] /= div;
            E[k][j] /= div;
        }

        // ��� �� �������, ������� ����������� ���� ��������
#pragma omp parallel for
        for (int i = k + 1; i < n; ++i)
        {
            // ���������� ��������� - ������� ��������� ������,
            // ������������� ��� ������������ ��������� ��������
            // ������
            double multi = Matrix[i][k];

            // �������� �� ��������� ������ ��������, ����������
            // �� ���������� ����� ��������� ��� � ��������,
            // ��� � � ��������� �������
            for (int j = 0; j < n; ++j)
            {
                Matrix[i][j] -= multi * Matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }

    // �������� �� �������� ����������� �������, ����������
    // �� ������ ����, ����� �����
    // �� ������ ����� ���������� �������� ���, � �� ��������
    // ������� ������������ ����������� ���������, � �� ��������� -
    // ��������
#pragma omp parallel for
    for (int k = n - 1; k > 0; --k)
    {
        // ��� �� �������, ������� ����������� ���� ��������
        for (int i = k - 1; i + 1 > 0; --i)
        {
            // ���������� ��������� - ������� ��������� ������,
            // ������������� ��� ������������ ��������� ��������
            // ������
            double multi = Matrix[i][k];

            // �������� �� ��������� ������ ��������, ����������
            // �� ���������� ����� ��������� ��� � ��������,
            // ��� � � ��������� �������
            for (int j = 0; j < n; ++j)
            {
                Matrix[i][j] -= multi * Matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }
    cout << "inverse matrix" << endl;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
            cout<<E[i][j] << " ";
        cout<<endl;
    }

    timeout = GetTickCount();
    timeres = timeout - timein;
    cout << "Our time: " << timeres << endl;
}




int main()
{


    Task4MatrixInverse(); Task4MatrixInverseParallel();

    system("pause");

    return 0;
}
*/