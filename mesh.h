#pragma once
#include <fstream>
#include <string>
#include <iostream>

using namespace std;

struct borderFirstCond
{
	int fun_num;
	int vertL;
	int vertR;
	int horizL;
	int horizR;
	int heightL;
	int heightR;
};

struct borderSecondCond
{
	int fun_num;
	int vertL;
	int vertR;
	int horizL;
	int horizR;
	int heightL;
	int heightR;
};

struct borderThirdCond
{
	int fun_num;
	int vertL;
	int vertR;
	int horizL;
	int horizR;
	int heightL;
	int heightR;
	double betta;
};

struct border
{
	borderFirstCond* firstCondAr;
	borderSecondCond* secondCondAr;
	borderThirdCond* thirdCondAr;
	int nFirstCond;
	int nSecondCond;
	int nThirdCond;
};

struct Area // ��������� ����������
{
	double lambda;
	double gamma;
	int lXW, rXW; // ������ ����� � ������ ������� ������� �� X
	int lYW, rYW; // ������ ����� � ������ ������� ������� �� Y
	int lZW, rZW; // ������ ����� � ������ ������� ������� �� Z
};

struct SepParam // ��������� ��������� ������� �����
{
	int n; // ���������� �������� ���������
	double q; // ����������� ���������
};

struct SpatioMesh // ���������������� �����
{
	double* XW, *YW, *ZW; // ���������� ������� �����
	int nW; // ����� ���������� ������� �����
	int nXW, nYW, nZW; // ���������� ������� �������� �� ��������������� ����

	int* IXW, *IYW, *IZW;

	double *X, *Y, *Z; // ���������� ���� ����� �����
	int n; // ����� ���������� �����
	int nX, nY, nZ; // ����� ���������� �������� �� �������������� ����
	
	int nArea; // ���������� �����������
	Area* areas; // ����������

	border borders;
	SpatioMesh(string fileNameCord, string fileNameSep, bool writeInCordFile = false)
	{
		ifstream in(fileNameCord);

		/*���������� �������� ��������*/
		in >> nXW >> nYW >> nZW; // ���������� ���������� �������� �������� �� ������������� ����
		
		nW = nXW * nYW * nZW; 

		XW = new double[nW];
		YW = new double[nW];
		ZW = new double[nW];

		int i;
		for (i = 0; i < nW; i++)
		{
			in >> XW[i] >> YW[i] >> ZW[i];
		}
		
		/*���������� ���������� � �����������*/
		in >> nArea;
		areas = new Area[nArea];
		for (i = 0; i < nArea; i++)
		{
			in >> areas[i].lambda >> areas[i].gamma;
			in >> areas[i].lXW >> areas[i].rXW;
			in >> areas[i].lYW >> areas[i].rYW;
			in >> areas[i].lZW >> areas[i].rZW;
		}
		in.close();

		cout << "XW" << endl;
		for (i = 0; i < nW; i++)
			cout << XW[i] << " ";
		cout << endl << endl;

		cout << "YW" << endl;
		for (i = 0; i < nW; i++)
			cout << YW[i] << " ";
		cout << endl << endl;

		cout << "ZW" << endl;
		for (i = 0; i < nW; i++)
			cout << ZW[i] << " ";
		cout << endl << endl;

		/*���������� ���������� ���������*/
		in.open(fileNameSep);
		int nSepX = nXW - 1, 
			nSepY = nYW - 1, 
			nSepZ = nZW - 1;
		
		SepParam* sepX = new SepParam[nSepX],
				* sepY = new SepParam[nSepY],
				* sepZ = new SepParam[nSepZ];

		
		IXW = new int[nXW];
		IYW = new int[nYW];
		IZW = new int[nZW];

		nX = 0;
		for (i = 0; i < nSepX; i++)
		{
			IXW[i] = nX;
			in >> sepX[i].n >> sepX[i].q;
			nX += sepX[i].n;
		}
		IXW[i] = nX;
		nX += 1;

		cout << "IXW" << endl;
		for (i = 0; i < nXW; i++)
			cout << IXW[i] << " ";
		cout << endl << endl;


		nY = 0;
		for (i = 0; i < nSepY; i++)
		{
			IYW[i] = nY;
			in >> sepY[i].n >> sepY[i].q;
			nY += sepY[i].n;
		}
		IYW[i] = nY;
		nY += 1;

		cout << "IYW" << endl;
		for (i = 0; i < nYW; i++)
			cout << IYW[i] << " ";
		cout << endl << endl;

		nZ = 0;
		for (i = 0; i < nSepZ; i++)
		{
			IZW[i] = nZ;
			in >> sepZ[i].n >> sepZ[i].q;
			nZ += sepZ[i].n;
		}
		IZW[i] = nZ;
		nZ += 1;

		cout << "IZW" << endl;
		for (i = 0; i < nZW; i++)
			cout << IZW[i] << " ";
		cout << endl << endl;

		int knotsInLine = nX; // ���������� ����� � �������������� �����
		int knotsInPlane = nY * knotsInLine; // ���������� ����� � ������������ �����
		int knotsN = nZ * knotsInPlane;

		int knotsWInLine = nXW; // ���������� �������� ����� � �������������� �����
		int knotsWInPlane = nYW * knotsWInLine; // ���������� �������� ����� � ������������ �����

		/*
		double heigtHorizPx, heigtHorizPy, heigtHorizPz;
		double heigtVertPx, heigtVertPy, heigtVertPz;
		double heigtPx, heigtPy, heigtPz;

		double heigtdxHorizPx, heigtdyHorizPy, heigtdzHorizPz;
		double heigtdxVertPx, heigdytVertPy, heigtdzVertPz;
		double heigtdxPx, heigtdyPy, heigtdzPz;

		double heigtHorizPdx, heigtHorizPdy, heigtHorizPdz;
		double heigtVertPdx, heigtVertPdy, heigtVertPdz;
		double heigtPdx, heigtPdy, heigtPdz;
		*/

		/*��������� �����*/
		
		int iXW, iYW, iZW;
		/*
		int inX, inY, inZ;
		double iqX, iqY, iqZ;
		int ix, iy, iz;
		*/
		X = new double[knotsN];
		Y = new double[knotsN];
		Z = new double[knotsN];

		int iKnot = 0;
		int iKnotPrev = 0;
		int c, l;

		/*
		double dxHoriz, dyHoriz, dzHoriz;
		double dxHeght, dyHeight, dzHeight;
		double dxVert, dyVert, dzVert;
		*/
		double qHoriz, qVert, qHeight;
		double nHoriz, nVert, nHeight;
		
		double koef;

		int pWCurr, // ������� ����� (�������� ����������)
			pWNextHeight, // ��������� ����� �� ������ (�������� ����������)
			pWNextHoriz, // ��������� �������������� ����� (�������� ����������)
			pWNextVert; // ��������� ������������ ����� (�������� ����������)

		int pCurr, // ������� ����� (����� ����������)
			pNextHeight, // ��������� ����� �� ������ (����� ����������)
			pNextHoriz, // ��������� �������������� ����� (����� ����������)
			pNextVert; // ��������� ������������ ����� (����� ����������)

		double dx, dy, dz;
		double x, y, z;
		int iX = 0, iY = 0, iZ = 0;
		int p;
		for (i = 0, iZW = 0; iZW < nZW; iZW++)
		{
			nHeight = sepZ[iZW].n;
			qHeight = sepZ[iZW].q;
			iZ = IZW[iZW];
			iY = 0;
			for (iYW = 0; iYW < nYW; iYW++)
			{
				nVert = sepY[iYW].n;
				qVert = sepY[iYW].q;
				iY = IYW[iYW];
				iX = 0;
				for (iXW = 0; iXW < nXW; iXW++)
				{
					pWCurr = iZW * knotsWInPlane + iYW * knotsWInLine + iXW;
					pWNextHoriz = pWCurr + 1;
					pWNextVert = pWCurr + knotsWInLine;
					pWNextHeight = pWCurr + knotsWInPlane;

					iX = IXW[iXW];
					pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
					pNextHoriz = pCurr + sepX[iXW].n;
					pNextVert = pCurr + sepY[iYW].n * knotsInLine;
					pNextHeight = pCurr + sepZ[iZW].n * knotsInPlane;

					/*�������� � ����������� �������������� ��������*/
					if (iXW != nXW - 1)
					{
						nHoriz = sepX[iXW].n;
						qHoriz = sepX[iXW].q;
						if (qHoriz == 1)
						{
							dx = (XW[pWNextHoriz] - XW[pWCurr]) / nHoriz;
							dy = (YW[pWNextHoriz] - YW[pWCurr]) / nHoriz;
							dz = (ZW[pWNextHoriz] - ZW[pWCurr]) / nHoriz;
						}
						else
						{
							koef = (1 - qHoriz) / (1 - pow(qHoriz, nHoriz));
							dx = (XW[pWNextHoriz] - XW[pWCurr]) * koef;
							dy = (YW[pWNextHoriz] - YW[pWCurr]) * koef;
							dz = (ZW[pWNextHoriz] - ZW[pWCurr]) * koef;
						}

						x = XW[pWCurr];
						y = YW[pWCurr];
						z = ZW[pWCurr];

						for (p = pCurr; p < pNextHoriz; p++)
						{
							X[p] = x;
							Y[p] = y;
							Z[p] = z;

							x += dx;
							y += dy;
							z += dz;

							dx *= qHoriz;
							dy *= qHoriz;
							dz *= qHoriz;
						}
					}

					/*�������� � ����������� ������������ ��������*/
					if (iYW != nYW - 1)
					{
						if (qVert == 1)
						{
							dx = (XW[pWNextVert] - XW[pWCurr]) / nVert;
							dy = (YW[pWNextVert] - YW[pWCurr]) / nVert;
							dz = (ZW[pWNextVert] - ZW[pWCurr]) / nVert;
						}
						else
						{
							koef = (1 - qVert) / (1 - pow(qVert, nVert));
							dx = (XW[pWNextVert] - XW[pWCurr]) * koef;
							dy = (YW[pWNextVert] - YW[pWCurr]) * koef;
							dz = (ZW[pWNextVert] - ZW[pWCurr]) * koef;
						}

						x = XW[pWCurr];
						y = YW[pWCurr];
						z = ZW[pWCurr];

						for (p = pCurr; p < pNextVert; p += knotsInLine)
						{
							X[p] = x;
							Y[p] = y;
							Z[p] = z;

							x += dx;
							y += dy;
							z += dz;

							dx *= qVert;
							dy *= qVert;
							dz *= qVert;
						}
					}
					

					/*�������� � ����������� �������� �� ������*/
					if (iZW != nZW - 1)
					{
						if (qHeight == 1)
						{
							dx = (XW[pWNextHeight] - XW[pWCurr]) / nHeight;
							dy = (YW[pWNextHeight] - YW[pWCurr]) / nHeight;
							dz = (ZW[pWNextHeight] - ZW[pWCurr]) / nHeight;
						}
						else
						{
							koef = (1 - qHeight) / (1 - pow(qHeight, nHeight));
							dx = (XW[pWNextHeight] - XW[pWCurr]) * koef;
							dy = (YW[pWNextHeight] - YW[pWCurr]) * koef;
							dz = (ZW[pWNextHeight] - ZW[pWCurr]) * koef;
						}

						x = XW[pWCurr];
						y = YW[pWCurr];
						z = ZW[pWCurr];

						for (p = pCurr; p < pNextHeight; p += knotsInPlane)
						{
							X[p] = x;
							Y[p] = y;
							Z[p] = z;

							x += dx;
							y += dy;
							z += dz;

							dx *= qHeight;
							dy *= qHeight;
							dz *= qHeight;
						}
					}
					
				}
			}
		}
		X[knotsN - 1] = XW[nW - 1];
		Y[knotsN - 1] = YW[nW - 1];
		Z[knotsN - 1] = ZW[nW - 1];

		for (i = 0, iZW = 0; iZW < nZW; iZW++)
		{
			iZ = IZW[iZW];
			for (iYW = 0; iYW < nYW - 1; iYW++)
			{
				for (iY = IYW[iYW] + 1; iY < IYW[iYW + 1]; iY++)
				{
					for (iXW = 0; iXW < nXW - 1; iXW++)
					{
						nHoriz = sepX[iXW].n;
						qHoriz = sepX[iXW].q;
						iX = IXW[iXW];
						pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
						pNextHoriz = pCurr + sepX[iXW].n;

						if (qHoriz == 1)
						{
							dx = (X[pNextHoriz] - X[pCurr]) / nHoriz;
							dy = (Y[pNextHoriz] - Y[pCurr]) / nHoriz;
							dz = (Z[pNextHoriz] - Z[pCurr]) / nHoriz;
						}
						else
						{
							koef = (1 - qHoriz) / (1 - pow(qHoriz, nHoriz));
							dx = (X[pNextHoriz] - X[pCurr]) * koef;
							dy = (Y[pNextHoriz] - Y[pCurr]) * koef;
							dz = (Z[pNextHoriz] - Z[pCurr]) * koef;
						}

						x = X[pCurr] + dx;
						y = Y[pCurr] + dy;
						z = Z[pCurr] + dz;

						dx *= qHoriz;
						dy *= qHoriz;
						dz *= qHoriz;

						for (p = pCurr + 1; p < pNextHoriz; p++)
						{
							X[p] = x;
							Y[p] = y;
							Z[p] = z;

							x += dx;
							y += dy;
							z += dz;

							dx *= qHoriz;
							dy *= qHoriz;
							dz *= qHoriz;
						}

					}
				}
			}
		}
		int j;

		for (i = 0, iZW = 0; iZW < nZW - 1; iZW++)
		{
			iZ = IZW[iZW];
			nHeight = sepZ[iZW].n;
			qHeight = sepZ[iZW].q;
			j = 0;
			pCurr = iZ * knotsInPlane;
			pNextHeight = pCurr + nHeight*knotsInPlane;
			
			for (iYW = 0; iYW < nYW - 1; iYW++)
			{
				for (iY = IYW[iYW]; iY < IYW[iYW + 1] + 1; iY++)
				{
					for (iXW = 0; iXW < nXW - 1; iXW++)
					{
						for (iX = IXW[iXW]; iX < IXW[iXW + 1] + 1; iX++, pNextHeight++)
						{
							if (iX == IXW[iXW] && iY == IYW[iYW]) continue;
							pCurr = iZ * knotsInPlane + iY * knotsInLine + iX;
							pNextHeight = pCurr + nHeight*knotsInPlane;
							if (qHeight == 1)
							{
								dx = (X[pNextHeight] - X[pCurr]) / nHeight;
								dy = (Y[pNextHeight] - Y[pCurr]) / nHeight;
								dz = (Z[pNextHeight] - Z[pCurr]) / nHeight;
							}
							else
							{
								koef = (1 - qHeight) / (1 - pow(qHeight, nHeight));
								dx = (X[pNextHeight] - X[pCurr]) * koef;
								dy = (Y[pNextHeight] - Y[pCurr]) * koef;
								dz = (Z[pNextHeight] - Z[pCurr]) * koef;
							}

							x = X[pCurr] + dx;
							y = Y[pCurr] + dy;
							z = Z[pCurr] + dz;

							dx *= qHeight;
							dy *= qHeight;
							dz *= qHeight;

							for (p = pCurr + knotsInPlane; p < pNextHeight; p += knotsInPlane)
							{
								X[p] = x;
								Y[p] = y;
								Z[p] = z;

								x += dx;
								y += dy;
								z += dz;

								dx *= qHeight;
								dy *= qHeight;
								dz *= qHeight;
							}

						}
					}
				}
			}
		}
	}

	int writeMesh(string fileName)
	{
		int i, j, k, p;
		ofstream out;
		out.open(fileName);
		out << nX << " " << nY << " " << nZ << endl;
		for (k = 0, p = 0; k < nZ; k++)
		{
			for (j = 0; j < nY; j++)
			{
				for (i = 0; i < nX; i++, p++)
					out << X[p] << " " << Y[p] << " " << Z[p] << "  ";
				out << endl;
			}
			out << endl;
		}
		return 0;
	}

	int readBorders(string fileName)
	{
		ifstream in;
		in.open(fileName);
		borderFirstCond borderFirstCond_buf;
		borderSecondCond  borderSecondCond_buf;
		borderThirdCond borderThirdCond_buf;

		borderFirstCond* borderFirstCondAr;
		borderSecondCond* borderSecondCondAr;
		borderThirdCond* borderThirdCondAr;

		int nFirstCond;
		in >> nFirstCond;
		borders.nFirstCond = nFirstCond;
		borders.firstCondAr = new borderFirstCond[nFirstCond];
		borderFirstCondAr = borders.firstCondAr;
		//borders = new border[N_bord];
		for (int i = 0; i < nFirstCond; i++)
		{
			in >> borderFirstCond_buf.fun_num;
			in >> borderFirstCond_buf.horizL;
			in >> borderFirstCond_buf.horizR;
			in >> borderFirstCond_buf.vertL;
			in >> borderFirstCond_buf.vertR;
			in >> borderFirstCond_buf.heightL;
			in >> borderFirstCond_buf.heightR;
			borderFirstCondAr[i] = borderFirstCond_buf;
		}

		int nSecondCond;
		in >> nSecondCond;
		borders.nSecondCond = nSecondCond;
		if (nSecondCond == 0)
			borders.secondCondAr = NULL;
		else
			borders.secondCondAr = new borderSecondCond[nFirstCond];
		borderSecondCondAr = borders.secondCondAr;
		//borders = new border[N_bord];
		for (int i = 0; i < nSecondCond; i++)
		{
			//in >> borderSecondCond_buf.fun_num;
			in >> borderSecondCond_buf.fun_num;
			in >> borderSecondCond_buf.horizL;
			in >> borderSecondCond_buf.horizR;
			in >> borderSecondCond_buf.vertL;
			in >> borderSecondCond_buf.vertR;
			in >> borderSecondCond_buf.heightL;
			in >> borderSecondCond_buf.heightR;
			borderSecondCondAr[i] = borderSecondCond_buf;
		}

		int nThirdCond;
		in >> nThirdCond;
		borders.nThirdCond = nThirdCond;
		if (nThirdCond == 0)
			borders.thirdCondAr = NULL;
		else
			borders.thirdCondAr = new borderThirdCond[nFirstCond];
		borderThirdCondAr = borders.thirdCondAr;
		//borders = new border[N_bord];
		for (int i = 0; i < nThirdCond; i++)
		{
			in >> borderThirdCond_buf.fun_num;
			in >> borderThirdCond_buf.horizL;
			in >> borderThirdCond_buf.horizR;
			in >> borderThirdCond_buf.vertL;
			in >> borderThirdCond_buf.vertR;
			in >> borderThirdCond_buf.heightL;
			in >> borderThirdCond_buf.heightR;
			in >> borderThirdCond_buf.betta;
			borderThirdCondAr[i] = borderThirdCond_buf;
		}
		return 0;
	}
private:
	int readSep(string fileNameSep)
	{

	}
};
