#ifndef CSF_H_
#define CSF_H_

#include "FCI.h"

#include <iostream>
#include <string>
#include <vector>

class CSF{
    //构建所有符合条件的CI组态，输入指定电子数，轨道数与自旋z分量，CI组态存于AI_Array中
    int CSF_new(int nelec_ex, int nOrb_ex, const double &S_ex, std::vector<int> &Orbital_ex, std::vector<CSF> &CSF_Array, int start);
    friend std::ostream &operator<<(std::ostream &os, const CSF &k);

    private:
        int nOrb = 0;
        int nelec = 0;
        double S = 0;
        double MS = 0;
        int MS2 = 0;
        int n_Slater_CI = 0;
        std::vector<int> Orbital;//一个四进制数，表示 0/1/2/3分别表示空轨道/自旋升高/自旋降低/双占
        std::vector<Slater_det> Slater_CI;
        std::vector<double> coefficient;

        bool Slater_create(const double &MS);

    public:
        CSF();                                   //默认构造函数
        CSF(const std::vector<int> &Orbital_ex, const double &S_ex, const double &MS_ex);                                   //默认构造函数,Orbital_ex存储三进制数表示电子占据情况
        ~CSF();                                  //默认析构函数
        bool CSF2Slater();
        bool cofcal();
        bool vector2Slater(std::vector<int> &Orbital_Slater, std::vector<int> &index, int nOrb_ex, int nalpha);
};


#endif