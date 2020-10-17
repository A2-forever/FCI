#include "FCI.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using std::vector;
using std::string;
using std::endl;


Slater_det::Slater_det()                                    //默认构造函数
{
}

Slater_det::Slater_det(const vector<int> &Orbital_ex)       //将轨道信息写入k轨道中
{
    nOrb = Orbital_ex.size() / 2;
    int nalpha = count(Orbital_ex.begin(), Orbital_ex.begin() + nOrb, 1);
    int nbeta = count(Orbital_ex.begin() + nOrb, Orbital_ex.end(), 1);
    nelec = nalpha + nbeta;
    MS = (double(nalpha - nbeta)) / 2.0; //(nalpha + nbeta)  = (nalpha - (nelec - nalpha))

    Orbital.resize(2 * nOrb);
    nelec_occ.resize(2 * nOrb);

    for (int i = 0; i < 2 * nOrb; i++)
        Orbital[i] = Orbital_ex[i];

    this->occ_cal();
}

Slater_det::~Slater_det()                                //默认析构函数
{
}

int Slater_det::Orb(const int &I, const int &sigma)
{
    return Orbital[sigma * nOrb + I];
}

bool Slater_det::occ_cal()                               //计算轨道i前的轨道的电子占据数
{
    bool flag=0;
    for (int sigma = 0; sigma < 2; sigma++)
    {
        nelec_occ[sigma * nOrb + 0] = 0;
        for (int I = 1; I < nOrb; I++)
            nelec_occ[sigma * nOrb + I] = nelec_occ[sigma * nOrb + I - 1] + Orbital[sigma * nOrb + I - 1]; //利用动态规划完成计算，sigma表示电子自旋
    }
    flag=1;
    return flag;
}

int Slater_det::gamma(const int &I,const int &sigma)
{
    //轨道的I_n的gamma值
    return sgn(nelec_occ[sigma * nOrb + I]);
}



//判断分子轨道i与分子轨道j之间占据情况不同的轨道与数量
//sigma表示电子自旋，0为alpha，1为beta
//index用于存储电子不同的位置，Num用于存储不同的占据情况不同的轨道数量
//eg:index[0][I][1]表示k1的I_beta位置占据电子，而k2的该位置未占据
//eg:Num[0]表示占据情况不同的alpha轨道的数量
//返回0表示出现问题
bool find(Slater_det &k1, Slater_det &k2, vector<int> &Num, vector<int> &index)
{
    if(k1.nOrb != k2.nOrb)
        return 0;
    
    int count1=0;
    int count2=0;
    bool flag=1;
    int nOrb = k1.nOrb;

    for (int sigma = 0; sigma < 2; sigma++)
    {
        count1=0;
        count2=0;
        for (int i = 0; i < nOrb; i++)
        {
            if(k1.Orbital[sigma * nOrb + i] != k2.Orbital[sigma * nOrb + i])
            {                                                                       //出现占据情况不同
                if(k1.Orbital[sigma * nOrb + i] == 1)                            //k1该位置占据电子，k2未占据
                {
                    index[0 * 2 * nOrb + sigma * nOrb + count1] = i;
                    count1++;
                }
                else                                                               //k2该位置占据电子，k1未占据
                {
                    index[1 * 2 * nOrb + sigma * nOrb + count2] = i;
                    count2++;
                }
            }
        }
        if (count1 != count2)
            flag = 0;
        Num[sigma] = count1;
    }

    return flag;
}

//构建所有符合条件的CI组态，输入指定电子数，轨道数与自旋z分量，CI组态存于AI_Array中
//ex表示现存的轨道或电子数
//Orbital_ex用于表示暂时表示组态的数组
//CI_Array用于存储组态
int CI_new(int nelec_ex, int nOrb_ex, const double &MS, vector<int> &Orbital_ex, vector<Slater_det> &CI_Array)
{ 
    //cout << nelec_ex <<" "<< nOrb_ex<< endl;
    if (fabs(nelec_ex) < 1e-6)
    {
        int nOrb = Orbital_ex.size() / 2;
        int nalpha = count(Orbital_ex.begin(), Orbital_ex.begin() + nOrb, 1);
        int nbeta = count(Orbital_ex.begin() + nOrb, Orbital_ex.end(), 1);
        double MS_ex = (double(nalpha - nbeta)) / 2; //(nalpha + nbeta) * 2 = (nalpha - (nelec - nalpha)) * 2

        if (fabs(MS_ex - MS) > 1e-6) //判断是否满足自旋条件
        {
            return 0;
        }
        else //满足自旋条件
        {
            Slater_det new_det(Orbital_ex);
            CI_Array.push_back(new_det);
            return 1;
        }
    }

    //本次使用的i用于计数，即标记占据最高能级的电子e_m的几种排列方式
    //现有nelec_ex个电子，n_Orb_ex个分子轨道，那么电子e_m可以从第nelec_ex号轨道排列到第n_Orb_ex号轨道
    //那么在vector中的序号自然为i-1
    int count = 0;
    for (int i = nelec_ex; i <= nOrb_ex; i++)
    {
        Orbital_ex[i - 1] = 1;
        count += CI_new(nelec_ex - 1, i - 1, MS, Orbital_ex, CI_Array);
        Orbital_ex[i - 1] = 0;
    }

    return count;
}


//输出Slater_det类的轨道占据情况
std::ostream &operator<<(std::ostream &os, const Slater_det &k)
{
    string E_sigma[2];
    E_sigma[0]="alpha: ";
    E_sigma[1]="beta:  ";
    os << endl;
    os << "nelec: " << k.nelec << endl;
    os << "nOrb:  " << k.nOrb << endl;
    os << "MS:  " << k.MS << endl;
    for (int sigma = 0; sigma < 2; sigma++)
    {
        os << E_sigma[sigma];
        for (int i = 0; i < k.nOrb; i++)
        {
            os << k.Orbital[sigma * k.nOrb + i] << "  ";
        }
        os << endl;
    }
    
    return os;
}



CI::CI()                                                 //默认构造函数
{
}

CI::CI(const double &h_nuc_ex, const vector<double> &h_ex, const vector<double> &g_ex, const int nOrb_ex)
{
    h_nuc = h_nuc_ex;
    h = h_ex;
    g = g_ex;
    nOrb = nOrb_ex;

    dim =nOrb_ex;
    dim2 = pow(nOrb_ex, 2);
    dim3 = pow(nOrb_ex, 3);

}

CI::~CI()                                                   //默认析构函数
{
}

double CI::get(const int &i = -1, const int &j = -1, const int &k = -1, const int &l = -1)
{
    if(k!=-1)
        return g[i * dim3 + j * dim2 + k * dim + l];
    else if (i != -1){
        return h[i * dim + j];
    }
    else{
        return h_nuc;
    }
}

//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
double CI::H_ij(Slater_det &k1, Slater_det &k2)
{
    vector<int> Num(2);//位置不同的电子的数量
    vector<int> index(2 * 2 * nOrb);//2个组态，2种自旋，最多nOrb个轨道
    bool flag = 1;
    flag = find(k1, k2, Num, index);
    if(!flag)
        return 0;

    if (Num[0] + Num[1] > 2)
        return 0;
    /*
	for (int i = 0; i < 2;i++)
		cout << Num[i] << "\t";
    cout << endl;
    for (int i = 0; i < 2;i++)
    {
        for (int j = 0; j < 2;j++)
        {
            for (int k = 0; k < Num[j];k++)
                cout << index[i * 2 * nOrb + k * 2 + j] << "\t";
            cout << endl;
        }
        cout << endl;
    }
    */

    double f_ij_alpha = this->F_ij(k1, k2, 0, Num, index);      //单电子耦合项的alpha部分
    //cout << f_ij_alpha << endl;
    double f_ij_beta = this->F_ij(k1, k2, 1, Num, index);       //单电子耦合项的beta部分
    //cout << f_ij_beta << endl;

    double g_ij = this->G_ij(k1, k2, Num, index);               //双电子耦合项
    //cout << g_ij << endl;

    //cout << f_ij_alpha + f_ij_beta + g_ij << endl;
    return f_ij_alpha+f_ij_beta+g_ij;

}


//计算分子轨道i与分子轨道j关于hamilton算符的耦合项单电子部分
//sigma表示电子自旋，0为alpha，1为beta
double CI::F_ij(Slater_det &k1, Slater_det &k2, const int sigma, const vector<int> &Num, const vector<int> &index)
{
    if (Num[1 - sigma] != 0) //在进行alpha部分的耦合时，beta部分必须完全一样，反之亦然
        return 0;

    double sum_int=0;
    if (Num[sigma] >= 2)                                        //k1与k2的n自旋轨道相差两个及以上的电子
    {
        sum_int = 0;
    }
    else if(Num[sigma] == 1)                                    //k1与k2的n自旋轨道相差一个电子
    {
        int I = index[1 * 2 * dim + sigma * dim + 0];           //k2的I_sigma分子轨道的n自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + sigma * dim + 0];           //k1的J_sigma分子轨道的n自旋轨道，位置占据电子，k2未占据
        sum_int = k2.gamma(I, sigma) * k1.gamma(J, sigma) * h[I * dim + J];
    }
    else if(Num[sigma] == 0)                                    //k1与k2的n自旋轨道相同 
    {
        for (int P = 0; P < nOrb; P++)
            sum_int += h[P * dim + P] * k1.Orb(P, sigma);
    }

    return sum_int;
}

//计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
//sigma表示电子自旋，0为alpha，1为beta
double CI::G_ij(Slater_det &k1, Slater_det &k2, const vector<int> &Num, const vector<int> &index)
{
    double g1_ij_alpha = this->G1_ij(k1, k2, 0, Num, index);
    double g1_ij_beta = this->G1_ij(k1, k2, 1, Num, index);

    double g2_ij = this->G2_ij(k1, k2, Num, index);

    return g1_ij_alpha + g1_ij_beta + g2_ij;

}
//计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
double CI::G1_ij(Slater_det &k1, Slater_det &k2, const int sigma, const vector<int> &Num, const vector<int> &index)
{
    if (Num[1 - sigma] != 0) //在进行alpha部分的耦合时，beta部分必须完全一样，反之亦然
        return 0;

    double sum_int = 0;
    if(Num[sigma] >= 3)                                         //k1与k2的n自旋轨道相差三个及以上的电子
    {
        sum_int=0;
    }
    else if(Num[sigma] == 2)                                    //k1与k2的n自旋轨道相差两个电子
    {
        int I = index[1 * 2 * dim + sigma * dim + 0];           //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[1 * 2 * dim + sigma * dim + 1];           //k2的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据
        int K = index[0 * 2 * dim + sigma * dim + 0];           //k1的K分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int L = index[0 * 2 * dim + sigma * dim + 1];           //k1的L分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        sum_int = k2.gamma(I, sigma) * k2.gamma(J, sigma) * k1.gamma(K, sigma) * k1.gamma(L, sigma) * (g[I * dim3 + K * dim2 + J * dim + L] - g[I * dim3 + L * dim2 + J * dim + K]);
    }
    else if(Num[sigma] == 1)                                    //k1与k2的n自旋轨道相差一个电子
    {
        int I = index[1 * 2 * dim + sigma * dim + 0];           //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + sigma * dim + 0];           //k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for (int R = 0; R < nOrb; R++)
            sum_int += k1.Orb(R, sigma) * (g[I * dim3 + J * dim2 + R * dim + R] - g[I * dim3 + R * dim2 + R * dim + J]);
        
        sum_int = k2.gamma(I, sigma) * k1.gamma(J, sigma) * sum_int;
    }
    else if(Num[sigma] == 0)                               //k1与k2的sigma自旋轨道相同
    {
        for (int P = 0; P < nOrb; P++)
            for (int R = 0; R < nOrb; R++)
                sum_int += k1.Orb(P, sigma) * k1.Orb(R, sigma) * (g[P * dim3 + P * dim2 + R * dim + R] - g[P * dim3 + R * dim2 + R * dim + P]);

        sum_int *= 0.5;
    }
    

    return sum_int;
    
}

double CI::G2_ij(Slater_det &k1, Slater_det &k2, const vector<int> &Num, const vector<int> &index)//计算分子轨道i与分子轨道j关于hamilton算符的耦合项
{
    double sum_int=0;

    if(Num[0] == 0 && Num[1] == 0)                            //k1与k2的自旋轨道相同
    {
        for (int P = 0; P < nOrb; P++)
            for (int R = 0; R < nOrb; R++)
                sum_int += k1.Orb(P, 0) * k1.Orb(R, 1) * g[P * dim3 + P * dim2 + R * dim + R];
    }
    else if(Num[0] == 1 && Num[1] == 0)                         //k1与k2的alpha自旋轨道相差一个电子,beta自旋轨道相同
    {
        int I = index[1 * 2 * dim + 0 * dim + 0];               //1*nOrb*nOrb+0*nOrb+0,k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + 0 * dim + 0];               //0*nOrb*nOrb+0*nOrb+0,k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for (int R = 0; R < nOrb; R++)
            sum_int += k1.Orb(R, 1) * g[I * dim3 + J * dim2 + R * dim + R];

        sum_int = k2.gamma(I, 0) * k1.gamma(J, 0) * sum_int;
    }
    else if(Num[0] == 0 && Num[1] == 1)                         //k1与k2的beta自旋轨道相差一个电子,alpha自旋轨道相同
    {
        int I = index[1 * 2 * dim + 1 * dim + 0];               //k2的I分子轨道的sigma自旋轨道，位置占据电子，k1未占据
        int J = index[0 * 2 * dim + 1 * dim + 0];               //k1的J分子轨道的sigma自旋轨道，位置占据电子，k2未占据

        for (int P = 0; P < nOrb; P++)
            sum_int += k1.Orb(P, 0) * g[P * dim3 + P * dim2 + I * dim + J];

        sum_int = k2.gamma(I, 1) * k1.gamma(J, 1) * sum_int;
    }
    else if(Num[0] == 1 && Num[1] == 1)                      //k1与k2的alpha与beta自旋轨道均相差一个电子
    {
        int I_alpha = index[1 * 2 * dim + 0 * dim + 0];      //k2的I分子轨道的alpha自旋轨道，位置占据电子，k1未占据
        int J_alpha = index[0 * 2 * dim + 0 * dim + 0];      //k1的J分子轨道的alpha自旋轨道，位置占据电子，k2未占据
        int I_beta = index[1 * 2 * dim + 1 * dim + 0];       //k2的I分子轨道的beta自旋轨道，位置占据电子，k1未占据
        int J_beta = index[0 * 2 * dim + 1 * dim + 0];       //k1的J分子轨道的beta自旋轨道，位置占据电子，k2未占据

        sum_int = k2.gamma(I_alpha, 0) * k1.gamma(J_alpha, 0) * k2.gamma(I_beta, 1) * k1.gamma(J_beta, 1) * g[I_alpha * dim3 + J_alpha * dim2 + I_beta * dim + J_beta];
    }
    else
    {
        sum_int=0;
    }

    return sum_int;
}
