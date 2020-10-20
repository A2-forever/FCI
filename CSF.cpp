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
    nelec = count(Orbital_ex.begin(), Orbital_ex.end(), 1) + count(Orbital_ex.begin(), Orbital_ex.end(), 2) + 2 * count(Orbital_ex.begin(), Orbital_ex.end(), 3);
    nOrb = Orbital_ex.size();
    S = S_ex;
    MS = MS_ex;
    MS2 = 2 * MS;
    for (int i = 0; i < nOrb; i++)
        Orbital[i] = Orbital_ex[i];
    
    this->CSF2Slater();
    this->cofcal();
    

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

    vector<int> index(nelec - n3);


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

    this->vector2Slater(Orbital_Slater, index, nelec -n3, nsigma[0]);
    n_Slater_CI = Slater_CI.size();

    

}

bool CSF::cofcal()
{
    vector<int> a(nOrb);
    vector<int> b(nOrb);
    vector<int> c(nOrb);
    vector<int> d(nOrb);

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
        for(int j = 0; j < nOrb; j++)
        {
            if(Slater_CI[i].Orb(j, 0))
                int delta = 1;
            else
                int delta = 0;
            
            if(d[j] == 0)
                f *= 1;
            else if(d[j] == 1)
            {
                if()
                f *= sqrt((a[j] + b[j] - gamma[j])/b[j]);
            }
            else if(d[j] == 2)
            {
                f *= 1;
            }
            else if(d[j] == 3)
            {
                f *= 1;
            }
        }
        coefficient[i] = f;
    }
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

int CSF_new(int nelec_ex, int nOrb_ex, const double &S_ex, std::vector<int> &Orbital_ex, std::vector<CSF> &CSF_Array, int start)
{
    //cout << nelec_ex <<" "<< nOrb_ex<< endl;
    if (nelec_ex < 1e-6)
    {   
        int n1 = count(Orbital_ex.begin(), Orbital_ex.end(), 1);//单占轨道数量（自旋下降）
        int n2 = count(Orbital_ex.begin(), Orbital_ex.end(), 2);//单占轨道数量（自旋上升）
        int n3 = count(Orbital_ex.begin(), Orbital_ex.end(), 3);//双占轨道数量
        double S = ((double)(n1 - n2))/2;
        if(abs(S - S_ex) > 1e-6)
            return 0;

        for(double MS = ((double)(n1 + n2))/2; MS > 0; MS--)
        {
            CSF new_CSF(Orbital_ex, S_ex, MS);
            CSF_Array.push_back(new_CSF);
        }
    }
    else if (start >= nOrb_ex)
    {
        return 0;
    }
    

    //本次使用的i用于表示电子占据数
    int ncount = 0;
    int n = 0;
    int n1 = count(Orbital_ex.begin(), Orbital_ex.begin() + start, 1);
    int n2 = count(Orbital_ex.begin(), Orbital_ex.begin() + start, 2);
    for(int i = 0; i<= 3; i++)
    {
        if(i == 2)
        {
            if(n1 < n2)
                continue;
            n = 1;
        }
        else if(i == 3)
            n = 2;
        else
            n = i;
            
        if(n <= nelec_ex)
        {
            Orbital_ex[start - 1] = i;
            ncount += CSF_new(nelec_ex - n, nOrb_ex, S_ex, Orbital_ex, CSF_Array, start + 1);
            Orbital_ex[start - 1] = 0;
        }
    }

    return ncount;
    
}

std::ostream &operator<<(std::ostream &os, const CSF &k)
{
    string E_sigma[2];
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
