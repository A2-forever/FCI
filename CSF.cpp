#include "CSF.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using std::vector;
using std::string;
using std::endl;


CSF::CSF::CSF()                                   //默认构造函数
{
}

CSF::CSF(const std::vector<int> &Orbital_ex, const double &S_ex, const double &MS_ex)                                   //默认构造函数
{
    nelec = count(Orbital_ex.begin(), Orbital_ex.end(), 1) + 2 * count(Orbital_ex.begin(), Orbital_ex.end(), 2);
    nOrb = Orbital_ex.size();
    S = S_ex;
    MS = MS_ex;
    MS2 = 2 * MS;
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
    vector<int> index;
    for(int i = 0; i < nOrb; i++)
    {
        if(Orbital[i] == 2)
        {
            Orbital_Slater[i] = 1;
            Orbital_Slater[nOrb + i] = 1;
            Orbital_S[i] = 3;
        }
        else if(Orbital[i] == 1)
        {
            index.push_back(i);
        }
    }

    vector<int> nsigma(2);
    int count = index.size();
    nsigma[0] = (count + MS2)/2;
    nsigma[1] = (count - MS2)/2;

    this->vector2Slater(Orbital_Slater, index, count, nsigma[0]);

    

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
        vector2Slater(Orbital_Slater, index, i - 1, nalpha - 1);
        Orbital_Slater[p] = 0;
    }
}

int CSF_new(int nelec_ex, int nOrb_ex, const double &S_ex, std::vector<int> &Orbital_ex, std::vector<CSF> &CSF_Array)
{
    //cout << nelec_ex <<" "<< nOrb_ex<< endl;
    if (nelec_ex < 1e-6)
    {   
        int n1 = count(Orbital_ex.begin(), Orbital_ex.end(), 1);//单占轨道数量（自旋下降）
        int n2 = count(Orbital_ex.begin(), Orbital_ex.end(), 2);//单占轨道数量（自旋上升）
        int n3 = count(Orbital_ex.begin(), Orbital_ex.end(), 3);//双占轨道数量
        int nelec = n1 + 2 * n2;
        int n_beta = (nelec - MS)/2;
        if( n2 > n_beta)
            return 0;

        int n_alpha = n1;//i表示使自旋增大的电子数
        int n_beta = 0;//j表示使自旋减小的电子数
        while (n_alpha > n_beta)
        {
            double S= (n_up - n_down)/2;
            CSF new_CSF(Orbital_ex, S, MS);
            CSF_Array.push_back(new_CSF);
            n_up--;
            n_down++;
        }
    }
    else if (nOrb_ex < 1e-6)
        return 0;

    //本次使用的i用于表示电子占据数
    int count = 0;
    for(int i = 0; i<= 2; i++)
    {
        if(i < nelec_ex)
        {
            Orbital_ex[nOrb_ex - 1] = i;
            count += CSF_new(nelec_ex - i, nOrb_ex - 1, S_ex, Orbital_ex, CSF_Array);
            Orbital_ex[nOrb_ex - 1] = 0;
        }
        else//到达剩余电子数之后，就可以不用做了
        {
            Orbital_ex[nOrb_ex - 1] = nelec_ex;
            count += CSF_new(nelec_ex - i, nOrb_ex - 1, S_ex, Orbital_ex, CSF_Array);
            Orbital_ex[nOrb_ex - 1] = 0;
            break;
        }
    }

    return count;
    
}

std::ostream &operator<<(std::ostream &os, const CSF &k)
{
    string E_sigma[2];
    os << endl;
    os << "nelec: " << k.nelec << endl;
    os << "nOrb:  " << k.nOrb << endl;
    os << "MS:  " << k.MS << endl;
    for (int i = 0; i < k.nOrb; i++)
    {
        os << k.Orbital[i] << "  ";
    }
    os << endl;
    
    return os;
}
