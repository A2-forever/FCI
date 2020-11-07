#include "CI.h"

#include <iostream>
#include <vector>
#include <algorithm>

using std::vector;
using std::string;
using std::cout;
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

int Slater_det::Gamma(const int &I,const int &sigma)
{
    //轨道的I_n的Gamma值
    return sgn(nelec_occ[sigma * nOrb + I]);
}

int Slater_det::gamma(const int &I, const int &sigma)
{
    int n = count(Orbital.begin() + sigma * nOrb,Orbital.begin() + sigma * nOrb + I, 1);
    return n;

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
                if(k1.Orbital[sigma * nOrb + i] == 1)                               //k1该位置占据电子，k2未占据
                {
                    index[0 * 2 * nOrb + sigma * nOrb + count1] = i;
                    count1++;
                }
                else                                                                //k2该位置占据电子，k1未占据
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


CSF::CSF()                                   //默认构造函数
{
}

CSF::CSF(const std::vector<int> &Orbital_ex, double S_ex, double MS_ex)                                   //默认构造函数
{
    nelec = count(Orbital_ex.begin(), Orbital_ex.end(), 1) + count(Orbital_ex.begin(), Orbital_ex.end(), 2) + 2 * count(Orbital_ex.begin(), Orbital_ex.end(), 3);
    nOrb = Orbital_ex.size();
    S = S_ex;
    MS = MS_ex;
    MS2 = 2 * MS;
    Orbital.resize(nOrb);
    for (int i = 0; i < nOrb; i++)
        Orbital[i] = Orbital_ex[i];
    
    this->CSF2Slater();
    

}

CSF::~CSF()                                  //默认析构函数
{
}

bool CSF::CSF2Slater()
{
    
    vector<int> Orbital_Slater(2 * nOrb);
    int n3 = 2 * count(Orbital.begin(), Orbital.end(), 3);//双占电子数

    vector<int> nsigma(2);
    nsigma[0] = (nelec -n3 + MS2)/2;
    nsigma[1] = (nelec -n3 - MS2)/2;

    vector<int> index;


    for(int i = 0; i < nOrb; i++)
    {
        if(Orbital[i] == 3)
        {
            Orbital_Slater[i] = 1;
            Orbital_Slater[nOrb + i] = 1;
        }
        else if(Orbital[i] == 1||Orbital[i] == 2)
        {
            index.push_back(i);
        }
    }
    //for(int i = 0; i < index.size(); i++)
        //cout<<index[i]<<"  ";
    //cout<<endl;

    this->vector2Slater(Orbital_Slater, index, nelec - n3, nsigma[0]);
    n_Slater_CI = Slater_CI.size();
    this->coefficient.resize(n_Slater_CI);
    this->coff_cal();

    return 1;

}

bool CSF::vector2Slater(std::vector<int> &Orbital_Slater, std::vector<int> &index, int nOrb_ex, int nalpha)
{
    if (nalpha < 1e-6)
    {
        int n = index.size();
        int n_up = (n + S)/2;
        int n_down = (n - S)/2;
        for(int i = 0; i < n; i++)
        {
            int p = index[i];
            
            if (Orbital_Slater[p] == 0)
                Orbital_Slater[nOrb + p] = 1;
        }
        Slater_det new_det(Orbital_Slater);
        Slater_CI.push_back(new_det);
        return 1;
        
    }

    for (int i = nalpha; i <= nOrb_ex; i++)
    {
        int p = index[i - 1];
        Orbital_Slater[p] = 1;
        vector2Slater(Orbital_Slater, index, i, nalpha - 1);
        Orbital_Slater[p] = 0;
    }

    return 1;
}

bool CSF::coff_cal()
{
    int a[nOrb] = {0};
    int b[nOrb] = {0};
    int c[nOrb] = {0};
    int d[nOrb] = {0};


    int p = Orbital[0];
    if(p == 1)
        b[0] = 1;
    else if(p == 3)
        a[0] = 1;
    
    c[0] = 1 - a[0] - b[0];
    d[0] = 3 * a[0] + b[0];
    for(int i = 1; i < nOrb; i++)
    {
        p = Orbital[i];
        if(p == 1)
            b[i] = b[i - 1] + 1;
        else if(p == 2)
        {
            a[i] = a[i - 1] + 1;
            b[i] = b[i - 1] - 1;
        }
        else if(p == 3)
            a[i] = a[i - 1] + 1;
        
        c[i] = i + 1 - a[i] - b[i];
        d[i] = 3 * (a[i] - a[i - 1]) + (b[i] - b[i - 1]);
    }

    
    for(int i = 0; i < n_Slater_CI; i++)
    {
        double f = 1;
        for(int k = 0; k < nOrb; k++)
        {
            int delta = Slater_CI[i].Orb(k, 0);
            
            if(d[k] == 0)
                f *= 1;
            else if(d[k] == 1)
            {
                int beta = Slater_CI[i].Orb(k, 1);
                int gamma = Slater_CI[i].gamma(k, beta);
                f *= sqrt(((double)(a[k] + b[k] - gamma))/(double)(b[k]));
            }
            else if(d[k] == 2)
            {
                int beta = Slater_CI[i].Orb(k, 1);
                int gamma = Slater_CI[i].gamma(k, beta);
                f *= sqrt(((double)(gamma - a[k] + 1))/(double)(b[k] + 2));
                f *= sgn(b[k] + delta);
            }
            else if(d[k] == 3)
            {
                f *= sgn(b[k]);
            }
        }
        coefficient[i] = f;
    }
    return 1;
}

int CSF_new(int nelec_ex, int nOrb_ex, double S_ex, std::vector<int> &Orbital_ex, std::vector<CSF> &CSF_Array, int start, double MS_ex)
{
    //std::cout << nelec_ex <<" "<< nOrb_ex<< endl;
    if (nelec_ex == 0)
    {   
        //for(int i = 0; i < 2; i++)
            //std::cout<<Orbital_ex[i]<<"  ";
        //std::cout<<std::endl;
        int n1 = count(Orbital_ex.begin(), Orbital_ex.end(), 1);//单占轨道数量（自旋下降）
        int n2 = count(Orbital_ex.begin(), Orbital_ex.end(), 2);//单占轨道数量（自旋上升）
        int n3 = count(Orbital_ex.begin(), Orbital_ex.end(), 3);//双占轨道数量
        double S = ((double)(n1 - n2))/2;
        if(fabs(S - S_ex) > 1e-6)
            return 0;


        //for(int i = 0; i < Orbital_ex.size(); i++)
            //std::cout<<Orbital_ex[i]<<"  ";
        //std::cout<<std::endl;

        int count = 0;
        if(fabs(MS_ex + 100) < 1e-6)
        {
            CSF new_CSF(Orbital_ex, S_ex, MS_ex);
            CSF_Array.push_back(new_CSF);
            count = 1;
        }
        else
        {        
            double MS = ((double)(n1 + n2))/2;
            while(count < n1 + n2 + 1)
            {
                //std::cout<<MS<<std::endl;
                CSF new_CSF(Orbital_ex, S_ex, MS);
                CSF_Array.push_back(new_CSF);
                MS--;
                count++;
            }
        }

        return count;

    }
    else if (start >= nOrb_ex)
    {
        return 0;
    }
    else if(nelec_ex < 0)
    {
        return 0;
    }


    //本次使用的i用于表示电子占据数
    int ncount = 0;
    int n = 0;
    int n1 = count(Orbital_ex.begin(), Orbital_ex.begin() + start, 1);
    int n2 = count(Orbital_ex.begin(), Orbital_ex.begin() + start, 2);
    for(int i = 0; i <= 3; i++)
    {
        if(i == 2)
        {
            if(n1 == n2)
                continue;
            
            n = 1;
        }
        else if(i == 3)
            n = 2;
        else
            n = i;
        
        if(n <= nelec_ex)
        {
            //std::cout<<start<<"  "<<i<<"  "<<nelec_ex - n<<"  "<<start + 1<<std::endl;
            Orbital_ex[start] = i;
            ncount += CSF_new(nelec_ex - n, nOrb_ex, S_ex, Orbital_ex, CSF_Array, start + 1);
            Orbital_ex[start] = 0;
        }
        else
            break;
    }

    return ncount;
    
}

std::ostream &operator<<(std::ostream &os, const CSF &k)
{
    os << endl;
    os << "nelec: " << k.nelec << endl;
    os << "nOrb:  " << k.nOrb << endl;
    os << "S:  " << k.S << endl;
    for (int i = 0; i < k.nOrb; i++)
    {
        os << k.Orbital[i] << "  ";
    }
    os << endl;
    
    return os;
}
