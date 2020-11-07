#include "File_string.h"
#include "FCI.h"
#include "CI.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


int main()
{
	clock_t startTime,endTime;
    startTime = clock();//计时开始

	//日志文件，输出一些中间结果与已经完成的模块
	ofstream logfile;
	logfile.open("data/FCI.log");
	logfile.setf(ios::showpoint); //设置为始终输出小数点后的数字，就是说 a = 3，它也输出 3.00000 这样
	logfile.precision(12);
	logfile.setf(ios::fixed); //设置为小数位始终有 6 位，没有这个的话就会像上面那个代码那样固定的不是小数点后面的数字了。

	//输入文件，文件内为几份文件FCIDUMP的名称，方便随时更换文件
	ifstream infile;
	string file_name;
	infile.open("data/FCI.input");
	if(infile.is_open())
        logfile << "Successfully open the input file" << endl;
	else{
        logfile << "fail to open the input file" << endl;
		return 0;
	}
	//获取读取文件所在文件夹的名称，一般是体系名
	getline(infile, file_name);
	infile.close();


	file_name = "data/" + file_name + "/FCIDUMP";
	int nelec = 0;							//电子数
	int nOrb = 0;							//轨道数
	double MS = 0;							//自旋z分量
	double h_nuc = 0;						//核积分项
	vector<double> h;						//单电子积分
	vector<double> g;						//双电子积分
	
	//读取积分文件中的积分值
	bool flag = read_int(file_name, nelec, nOrb, MS, h_nuc, h, g);
	if(flag){
		logfile << "successfully read the intergal file" << endl;
	}
	else
	{
		logfile << "fail to read the intergal file" << endl;
		return 0;
	}

	for (int i = 0; i < nOrb; i++)
	{
		for (int j = 0; j < nOrb; j++)
		{
			for (int k = 0; k < nOrb; k++)
			{
				for (int l = 0; l < nOrb; l++)
				{
					logfile << i + 1 << " " << j + 1 << " " << k + 1 << " " << l + 1 << "\t\t" << g[i * pow(nOrb, 3) + j * pow(nOrb, 2) + k * nOrb + l] << endl;
					//if(fabs(g[i * pow(nOrb, 3) + j * pow(nOrb, 2) + k * nOrb + l])>1e-20)
					//	cout << "bad" << endl;
				}

			}
        }
    }

	cout << g.size() << endl
		 << h.size() << endl;
	logfile <<endl;
	logfile << h << endl;
	logfile << h_nuc << endl;
	logfile <<endl;

	logfile << "the number of electrons: " << nelec << endl
		 << "the number of active orbitals: " << nOrb << endl
		 << "MS: " << MS << endl;
	cout << "the number of electrons: " << nelec << endl
		 << "the number of active orbitals: " << nOrb << endl
		 << "MS: " << MS << endl;


	//创建FCI类，存储分子轨道积分
	FCI FCI(h_nuc, h, g, nOrb);
	//创建组态，CI_Array用于存储组态，Orbital表示暂时表示组态的数组，没有实际用处
	vector<Slater_det> CI_Array;
	vector<int> Orbital(2 * nOrb);
	int nCI = CI_new(nelec, 2 * nOrb, MS, Orbital, CI_Array);


	double S = 0;
	int start = 0;
	vector<CSF> CSF_Array;
	vector<int> Orbital_CSF(nOrb);
	int nCSF = CSF_new(nelec, nOrb, S, Orbital_CSF, CSF_Array, start);


	logfile << "the number of CIs: " << nCI << endl;
	cout << "the number of CIs: " << nCI << endl;
	cout << "the number of CSFs: " << nCSF << endl;

	for (int i = 0; i < nCI; i++)
		logfile << CI_Array[i];

	//构建Hamilton矩阵
	double temp = 0;
	MatrixXd H(nCI, nCI);
	MatrixXd H2(nCSF, nCSF);
	for (int i = 0; i < nCI; i++)
	{
		H(i, i) = FCI.H_ij(CI_Array[i], CI_Array[i]);
		for (int j = i+1; j < nCI; j++)
		{
			temp = FCI.H_ij(CI_Array[i], CI_Array[j]);
			H(i, j) = temp;
			H(j, i) = temp;
		}
	}

	for (int i = 0; i < nCSF; i++)
	{
		H2(i, i) = FCI.H_ij(CSF_Array[i], CSF_Array[i]);
		for (int j = i+1; j < nCSF; j++)
		{
			temp = FCI.H_ij(CSF_Array[i], CSF_Array[j]);
			H2(i, j) = temp;
			H2(j, i) = temp;
		}
	}
	logfile << "Hamilton matrix: " << endl
			<< H << endl;
	logfile << "Hamilton matrix: " << endl
			<< H2 << endl;

	/*
	//构建Hamilton矩阵
	double temp = 0;
	vector<double> H(nCI * nCI);
	for (int i = 0; i < nCI; i++){
		H[i * nCI + i] = FCI.H_ij(CI_Array[i], CI_Array[i]) + h_nuc;
		for (int j = i+1; j < nCI; j++){
			temp = FCI.H_ij(CI_Array[i], CI_Array[j]);
			H[i * nCI + j] = temp;
			H[j * nCI + i] = temp;
		}
	}
	logfile << "Hamilton matrix";
	logfile << H;

	vector<double> F_alpha(nCI * nCI);
	vector<double> F_beta(nCI * nCI);
	vector<double> G1_alpha(nCI * nCI);
	vector<double> G1_beta(nCI * nCI);
	vector<double> G2(nCI * nCI);

	for (int i = 0; i < nCI; i++){
		for (int j = i; j < nCI; j++){
			vector<int> Num(2);//位置不同的电子的数量
    		vector<int> index(2 * 2 * nOrb);//2个组态，2种自旋，最多nOrb个轨道
			find(CI_Array[i], CI_Array[j], Num, index);
			temp = FCI.F_ij(CI_Array[i], CI_Array[j], 0, Num, index);
			F_alpha[i * nCI + j] = temp;
			F_alpha[j * nCI + i] = temp;

			temp = FCI.F_ij(CI_Array[i], CI_Array[j], 1, Num, index);
			F_beta[i * nCI + j] = temp;
			F_beta[j * nCI + i] = temp;
			
			temp = FCI.G1_ij(CI_Array[i], CI_Array[j], 0, Num, index);
			G1_alpha[i * nCI + j] = temp;
			G1_alpha[j * nCI + i] = temp;
			
			temp = FCI.G1_ij(CI_Array[i], CI_Array[j], 1, Num, index);
			G1_beta[i * nCI + j] = temp;
			G1_beta[j * nCI + i] = temp;
			
			temp = FCI.G2_ij(CI_Array[i], CI_Array[j], Num, index);
			G2[i * nCI + j] = temp;
			G2[j * nCI + i] = temp;
		}
	}

	logfile << endl
			<< "F_alpha" << F_alpha;

	logfile << endl
			<< "F_beta" << F_beta;

	logfile << endl
			<< "G1_alpha" << G1_alpha;

	logfile << endl
			<< "G1_beta" << G1_beta;

	logfile << endl
			<< "G2" << G2;
	CI_Array.clear();
	

	//对角化Hamilton矩阵
	
	int nJt = 30;
	vector<double> dbVectors(nCI * nCI);//输入一个nCI x nCI大小的矩阵
	vector<double> dbEigenvalues(nCI);//特征值
	//flag=eigen(H, dbVectors, dbEigenvalues);
	flag=eigenh(H, dbVectors, dbEigenvalues, nJt);
	if(flag){
		logfile << "Successfully diagonalize the Hamilton matrix" << endl;
	}
	else{
		logfile << "fail to diagonalize the Hamilton matrix,,you need more iterations" << endl;
	}
	
*/
	//求本征值
	EigenSolver<MatrixXd> eig(H);
	MatrixXd evecs = eig.eigenvectors().real();
	MatrixXd evals = eig.eigenvalues().real();

	MatrixXd::Index minRow, minCol;
	evals.minCoeff(&minRow, &minCol);


	ofstream outfile;
	outfile.open("data/FCI.output");
	if (outfile.is_open())
	{
		logfile << "Successfully open the output file" << endl;
		outfile.setf(ios::showpoint); //设置为始终输出小数点后的数字，就是说 a = 3，它也输出 3.00000 这样
		outfile.precision(8);
		outfile.setf(ios::fixed); //设置为小数位始终有 6 位，没有这个的话就会像上面那个代码那样固定的不是小数点后面的数字了。
	}
	else
	{
        logfile << "fail to open the output file" << endl;
		return 0;
	}
	endTime = clock();//计时结束
	cout << "Running time:" << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s";

	outfile << "Running time: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl
			<< endl;
	outfile << "nCI: " << nCI << endl
			<< endl;
	//输出最小本征值
	outfile << "Eigenvalues:" << endl
			<< evals(minRow, minCol) + h_nuc << endl
			<< endl;
	outfile << "Eigenvectors: " << endl
			<< evecs.col(minRow) << endl;

	/*
	for (int i = 0; i < nCI; i++)
	{
		for (int j = 0; j < nCI;j++)
		{
			outfile << H[i * nCI + j] << "\t";
		}
		outfile << endl;
	}


	outfile << "Eigenvalues:" << endl;
	for (int i = 0; i < nCI;i++)
		outfile << dbEigenvalues[i] << endl;
	outfile << endl;


	outfile << "Eigenvectors: ";
	outfile << dbVectors;
	outfile << endl;
	*/

	/*
	//寻找本征值中的最小值的位置，然后输出对应本征向量
	vector<double>::iterator min = min_element(dbEigenvalues.begin(), dbEigenvalues.end());
	int index_min = distance(dbEigenvalues.begin(), min);
	outfile << "Energy:";
	outfile << dbEigenvalues[index_min] << endl;


	outfile << "Eigenvectors: ";
	for (int i = 0; i < nCI;i++)
		outfile << dbVectors[i * nCI + index_min] << "\t";
	outfile << endl;
*/
	outfile.close();
	logfile.close();

	cin >> flag;
}