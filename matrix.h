#pragma once
#include <iostream>

using namespace std;
//int getGlobeNum(int i, int j, int k, int n_x, int n_y, int n_z) { return k * n_x * n_y + j * n_x + i; }


struct matrix {
	int* ig; // массив указателей начала строк (столбцов) в массивах ggl и ggu
	int* jg; // массив номеров строк(столбцов) заданного треугольника матрицы
	double* gg; // массив внедиагональных элементов нижнего треугольника матрицы
	double* di; // массив диагональных элементов матрицы
	int n; // размерность матрицы
	int n_gg; // размер массивов gg и jg
	double* ggl;
	double* ggu;
	bool* isUntouch;
	matrix(matrix& Rec)
	{
		ig = Rec.ig;
		jg = Rec.jg;
		n = Rec.n;
		n_gg = Rec.n_gg;
		ggl = new double[n_gg];
		ggu = new double[n_gg];
		di = new double[n];
	}
	matrix(int n_x, int n_y, int n_z) // Построение профиля
	{
		n = (n_x) * (n_y) * (n_z);
		n_gg = (n_x - 1) + 5 * (n_y - 1) + (n_x - 2) * (n_y - 1) * 4; // Добавил 1 потому что вроде надо, а вообще лучше проверить
		n_gg += (n_z - 1) * (22 + (n_x + n_y - 4) * 17 + (n_x - 2) * (n_y - 2) * 13);
		//di = new double[2 * n];
		di = new double[n];
		isUntouch = new bool[n];
		ig = new int[n + 1];
		//gg = new double[2 * n_gg];
		gg = new double[n_gg];
		jg = new int[n_gg];
		int ni_ig = 0;
		ig[0] = ni_ig;
		ig[1] = ni_ig;
		int i_ig = 2;
		int i_jg = 0;
		int i, j, k;
		ni_ig++;
		/*Отдельное заполнение профиля для плоскости z_0*/
		for (i = 1; i < n_x; i++, ni_ig++, i_ig++, i_jg++)
		{
			ig[i_ig] = ni_ig;
			jg[i_jg] = i - 1;
		}
		ni_ig--;

		for (j = 1; j < n_y; j++)
		{
			/*Заполнение для x_0*/
			ni_ig += 2;
			ig[i_ig] = ni_ig;
			i_ig++;
			jg[i_jg] = getGlobeNumb(0, j - 1, 0, n_x, n_y, n_z);
			jg[i_jg + 1] = getGlobeNumb(1, j - 1, 0, n_x, n_y, n_z);
			i_jg += 2;
			/*Заполнение для x_i*/
			for (i = 1; i < n_x - 1; i++)
			{
				ni_ig += 4;
				ig[i_ig] = ni_ig;
				i_ig++;
				jg[i_jg] = getGlobeNumb(i - 1, j - 1, 0, n_x, n_y, n_z);
				jg[i_jg + 1] = getGlobeNumb(i, j - 1, 0, n_x, n_y, n_z);
				jg[i_jg + 2] = getGlobeNumb(i + 1, j - 1, 0, n_x, n_y, n_z);
				jg[i_jg + 3] = getGlobeNumb(i - 1, j, 0, n_x, n_y, n_z);
				i_jg += 4;
			}
			ni_ig += 3;
			ig[i_ig] = ni_ig;
			i_ig++;
			jg[i_jg] = getGlobeNumb(i - 1, j - 1, 0, n_x, n_y, n_z);
			jg[i_jg + 1] = getGlobeNumb(i, j - 1, 0, n_x, n_y, n_z);
			jg[i_jg + 2] = getGlobeNumb(i - 1, j, 0, n_x, n_y, n_z);
			i_jg += 3;
		}

		/*Заполнение для плоскостей z_i*/
		for (k = 1; k < n_z; k++)
		{
			ni_ig += 4;
			ig[i_ig] = ni_ig;
			i_ig++;
			jg[i_jg] = getGlobeNumb(0, 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 1] = getGlobeNumb(1, 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 2] = getGlobeNumb(0, 1, k - 1, n_x, n_y, n_z);
			jg[i_jg + 3] = getGlobeNumb(1, 1, k - 1, n_x, n_y, n_z);
			i_jg += 4;

			for (i = 1; i < n_x - 1; i++)
			{
				ni_ig += 7;
				ig[i_ig] = ni_ig;
				i_ig++;

				jg[i_jg + 0] = getGlobeNumb(i - 1, 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 1] = getGlobeNumb(i + 0, 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 2] = getGlobeNumb(i + 1, 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 3] = getGlobeNumb(i - 1, 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 4] = getGlobeNumb(i + 0, 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 5] = getGlobeNumb(i + 1, 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 6] = getGlobeNumb(i - 1, 0, k, n_x, n_y, n_z);
				i_jg += 7;
			}

			ni_ig += 5;
			ig[i_ig] = ni_ig;
			i_ig++;
			jg[i_jg + 0] = getGlobeNumb(i - 1, 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 1] = getGlobeNumb(i - 0, 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 2] = getGlobeNumb(i - 1, 1, k - 1, n_x, n_y, n_z);
			jg[i_jg + 3] = getGlobeNumb(i - 0, 1, k - 1, n_x, n_y, n_z);
			jg[i_jg + 4] = getGlobeNumb(i - 1, 0, k + 0, n_x, n_y, n_z);
			i_jg += 5;

			for (j = 1; j < n_y - 1; j++)
			{
				ni_ig += 8;
				ig[i_ig] = ni_ig;
				i_ig++;
				jg[i_jg + 0] = getGlobeNumb(0, j - 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 1] = getGlobeNumb(1, j - 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 2] = getGlobeNumb(0, j + 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 3] = getGlobeNumb(1, j + 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 4] = getGlobeNumb(0, j + 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 5] = getGlobeNumb(1, j + 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 6] = getGlobeNumb(0, j - 1, k + 0, n_x, n_y, n_z);
				jg[i_jg + 7] = getGlobeNumb(1, j - 1, k + 0, n_x, n_y, n_z);
				i_jg += 8;

				for (i = 1; i < n_x - 1; i++)
				{
					ni_ig += 13;
					ig[i_ig] = ni_ig;
					i_ig++;
					jg[i_jg + 0] = getGlobeNumb(i - 1, j - 1, k - 1, n_x, n_y, n_z);
					jg[i_jg + 1] = getGlobeNumb(i + 0, j - 1, k - 1, n_x, n_y, n_z);
					jg[i_jg + 2] = getGlobeNumb(i + 1, j - 1, k - 1, n_x, n_y, n_z);
					jg[i_jg + 3] = getGlobeNumb(i - 1, j + 0, k - 1, n_x, n_y, n_z);
					jg[i_jg + 4] = getGlobeNumb(i + 0, j + 0, k - 1, n_x, n_y, n_z);
					jg[i_jg + 5] = getGlobeNumb(i + 1, j + 0, k - 1, n_x, n_y, n_z);
					jg[i_jg + 6] = getGlobeNumb(i - 1, j + 1, k - 1, n_x, n_y, n_z);
					jg[i_jg + 7] = getGlobeNumb(i + 0, j + 1, k - 1, n_x, n_y, n_z);
					jg[i_jg + 8] = getGlobeNumb(i + 1, j + 1, k - 1, n_x, n_y, n_z);
					jg[i_jg + 9] = getGlobeNumb(i - 1, j - 1, k + 0, n_x, n_y, n_z);
					jg[i_jg + 10] = getGlobeNumb(i + 0, j - 1, k + 0, n_x, n_y, n_z);
					jg[i_jg + 11] = getGlobeNumb(i + 1, j - 1, k + 0, n_x, n_y, n_z);
					jg[i_jg + 12] = getGlobeNumb(i - 1, j + 0, k + 0, n_x, n_y, n_z);
					i_jg += 13;
				}

				ni_ig += 9;
				ig[i_ig] = ni_ig;
				i_ig++;
				jg[i_jg + 0] = getGlobeNumb(i - 1, j - 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 1] = getGlobeNumb(i + 0, j - 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 2] = getGlobeNumb(i - 1, j + 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 3] = getGlobeNumb(i + 0, j + 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 4] = getGlobeNumb(i - 1, j + 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 5] = getGlobeNumb(i + 0, j + 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 6] = getGlobeNumb(i - 1, j - 1, k + 0, n_x, n_y, n_z);
				jg[i_jg + 7] = getGlobeNumb(i + 0, j - 1, k + 0, n_x, n_y, n_z);
				jg[i_jg + 8] = getGlobeNumb(i - 1, j + 0, k + 0, n_x, n_y, n_z);
				i_jg += 9;
			}

			/*y - final*/
			ni_ig += 6;
			ig[i_ig] = ni_ig;
			i_ig++;
			jg[i_jg] = getGlobeNumb(0, j - 1, k - 1, n_x, n_y, n_z);
			jg[i_jg + 1] = getGlobeNumb(1, j - 1, k - 1, n_x, n_y, n_z);
			jg[i_jg + 2] = getGlobeNumb(0, j + 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 3] = getGlobeNumb(1, j + 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 4] = getGlobeNumb(0, j - 1, k + 0, n_x, n_y, n_z);
			jg[i_jg + 5] = getGlobeNumb(1, j - 1, k + 0, n_x, n_y, n_z);
			i_jg += 6;

			for (i = 1; i < n_x - 1; i++)
			{
				ni_ig += 10;
				ig[i_ig] = ni_ig;
				i_ig++;

				jg[i_jg + 0] = getGlobeNumb(i - 1, j - 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 1] = getGlobeNumb(i + 0, j - 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 2] = getGlobeNumb(i + 1, j - 1, k - 1, n_x, n_y, n_z);
				jg[i_jg + 3] = getGlobeNumb(i - 1, j + 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 4] = getGlobeNumb(i + 0, j + 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 5] = getGlobeNumb(i + 1, j + 0, k - 1, n_x, n_y, n_z);
				jg[i_jg + 6] = getGlobeNumb(i - 1, j - 1, k + 0, n_x, n_y, n_z);
				jg[i_jg + 7] = getGlobeNumb(i + 0, j - 1, k + 0, n_x, n_y, n_z);
				jg[i_jg + 8] = getGlobeNumb(i + 1, j - 1, k + 0, n_x, n_y, n_z);
				jg[i_jg + 9] = getGlobeNumb(i - 1, j + 0, k + 0, n_x, n_y, n_z);
				i_jg += 10;
			}

			ni_ig += 7;
			ig[i_ig] = ni_ig;
			i_ig++;
			jg[i_jg + 0] = getGlobeNumb(i - 1, j - 1, k - 1, n_x, n_y, n_z);
			jg[i_jg + 1] = getGlobeNumb(i + 0, j - 1, k - 1, n_x, n_y, n_z);
			jg[i_jg + 2] = getGlobeNumb(i - 1, j + 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 3] = getGlobeNumb(i + 0, j + 0, k - 1, n_x, n_y, n_z);
			jg[i_jg + 4] = getGlobeNumb(i - 1, j - 1, k + 0, n_x, n_y, n_z);
			jg[i_jg + 5] = getGlobeNumb(i + 0, j - 1, k + 0, n_x, n_y, n_z);
			jg[i_jg + 6] = getGlobeNumb(i - 1, j + 0, k + 0, n_x, n_y, n_z);
			i_jg += 7;
		}
		ggl = gg;
		ggu = gg;
		cout << "n = " << endl;
		cout << "teor n_gg = " << n_gg << endl;
		cout << "prac n_gg = " << ni_ig << endl;
		cout << "prac i_jg = " << i_jg << endl;

	}
	void printMatrix()
	{
		cout << "Matrix" << endl;
		cout << "\tig" << endl;
		for (int i = 0; i < n + 1; i++)
		{
			cout << ig[i] << " ";
		}
		cout << endl;

		cout << "\tjg" << endl;
		for (int i = 0; i < n_gg; i++)
		{
			cout << jg[i] << " ";
		}
		cout << endl;
		
		cout << "\tdi" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << di[i] << " ";
		}
		cout << endl;

		cout << "\tgg" << endl;
		for (int i = 0; i < n_gg; i++)
		{
			cout << gg[i] << " ";
		}
		cout << endl;
	}
	void init()
	{
		for (int i = 0; i < n; i++)
		{
			di[i] = 0;
			isUntouch[i] = true;
		}

		for (int i = 0; i < n_gg; i++)
			gg[i] = 0;
		
	}
	int getGlobeNumb(int i, int j, int k, int nX, int nY, int nZ) { return k * nX * nY + j * nX + i; }
};
