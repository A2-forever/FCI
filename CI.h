#ifndef CI_H_
#define CI_H_

#define sgn(n) (n)%2 == 0 ? 1:-1//奇数返回-1，偶数返回1

#include <iostream>
#include <vector>

class Slater_det{
    //判断分子轨道i与分子轨道j之间alpha或beta电子的差别,index存放位置，N表示alpha和beta电子的数目
    friend bool find(Slater_det &k1, Slater_det &k2, std::vector<int> &Num, std::vector<int> &index);
    //构建所有符合条件的CI组态，输入指定电子数，轨道数与自旋z分量，CI组态存于AI_Array中
    friend int CI_new(int nelec_ex, int nOrb_ex, const double &MS, std::vector<int> &Orbital_ex, std::vector<Slater_det> &CI_Array);
    friend std::ostream &operator<<(std::ostream &os, const Slater_det &k);
    friend class CSF;

    private:
        int nOrb = 0;               //分子轨道数
        int nelec = 0;              //占据电子数
        double MS = 0;              //自旋角动量的z分量
        std::vector<int> Orbital;   //1表示alpha轨道，2表示beta轨道,此数组表示分子轨道的电子占据，0代表空，1代表占据电子
        std::vector<int> nelec_occ; //1表示alpha轨道，2表示beta轨道，此数组表示分子轨道i前的轨道的电子占据数

    public:
        Slater_det();                                   //默认构造函数
        Slater_det(const std::vector<int> &Orbital_ex); //将轨道信息写入k轨道中表示占据电子数
        ~Slater_det();                                  //默认析构函数
        int Orb(const int &I, const int &sigma);        //用于调用分子轨道的占据情况
        bool occ_cal();                                 //计算轨道i前的轨道的电子占据数
        int gamma(const int &I, const int &sigma);
        int Gamma(const int &I, const int &sigma);      //组态，I_sigma位置的Gamma

};

class CSF{
    //构建所有符合条件的CI组态，输入指定电子数，轨道数与自旋z分量，CI组态存于AI_Array中
    friend int CSF_new(int nelec_ex, int nOrb_ex, double S_ex, std::vector<int> &Orbital_ex, std::vector<CSF> &CSF_Array, int start);
    friend std::ostream &operator<<(std::ostream &os, const CSF &k);
    friend class FCI;

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

    public:
        CSF();                                   //默认构造函数
        CSF(const std::vector<int> &Orbital_ex, double S_ex, double MS_ex);                                   //默认构造函数,Orbital_ex存储三进制数表示电子占据情况
        ~CSF();                                  //默认析构函数
        bool CSF2Slater();
        bool vector2Slater(std::vector<int> &Orbital_Slater, std::vector<int> &index, int nOrb_ex, int nalpha);
        bool coff_cal();
};


#endif