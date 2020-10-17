#ifndef FCI_H_
#define FCI_H_

#define sgn(n) n%2==0 ? 1:-1//奇数返回-1，偶数返回1

#include <iostream>
#include <string>
#include <vector>

class Slater_det{
    //判断分子轨道i与分子轨道j之间alpha或beta电子的差别,index存放位置，N表示alpha和beta电子的数目
    friend bool find(Slater_det &k1, Slater_det &k2, std::vector<int> &Num, std::vector<int> &index);
    //构建所有符合条件的CI组态，输入指定电子数，轨道数与自旋z分量，CI组态存于AI_Array中
    friend int CI_new(int nelec_ex, int nOrb_ex, const double &MS, std::vector<int> &Orbital_ex, std::vector<Slater_det> &CI_Array);
    friend std::ostream &operator<<(std::ostream &os, const Slater_det &k);

    private:
        int nOrb = 0;               //分子轨道数
        int nelec = 0;              //占据电子数
        double MS = 0;              //自旋角动量的z分量
        std::vector<int> Orbital;   //1表示alpha轨道，2表示beta轨道,此数组表示分子轨道的电子占据，0代表空，1代表占据电子
        std::vector<int> nelec_occ; //1表示alpha轨道，2表示beta轨道，此数组表示分子轨道i前的轨道的电子占据数
                                    //eg:nelec[i][1]代表分子轨道i前的占据的beta电子数

    public:
        Slater_det();                                   //默认构造函数
        Slater_det(const std::vector<int> &Orbital_ex); //将轨道信息写入k轨道中表示占据电子数
        ~Slater_det();                                  //默认析构函数
        int Orb(const int &I, const int &sigma);        //用于调用分子轨道的占据情况
        bool occ_cal();                                 //计算轨道i前的轨道的电子占据数
        int gamma(const int &I, const int &sigma);      //组态，I_sigma位置的Gamma

};


class CI{
    private:
        int nOrb = 0;          //分子轨道数
        double h_nuc = 0;      //核积分项
        std::vector<double> h; //单电子积分
        std::vector<double> g; //双电子积分

        int dim = 0;
        int dim2 = 0;
        int dim3 = 0;

    public:
        CI(); //默认构造函数
        CI(const double &h_nuc_ex, const std::vector<double> &h_ex, const std::vector<double> &g_ex, const int nOrb_ex);
        ~CI(); //默认析构函数
        double get(const int &i, const int &j, const int &k, const int &l);

        //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
        double H_ij(Slater_det &k1, Slater_det &k2);
        double H(Slater_det &k1, Slater_det &k2);

        //计算分子轨道i与分子轨道j关于hamilton算符的耦合项单电子部分
        //sigma表示电子自旋，0为alpha，1为beta
        double F_ij(Slater_det &k1, Slater_det &k2, const int sigma, const std::vector<int> &Num, const std::vector<int> &index);

        //计算分子轨道i与分子轨道j关于hamilton算符的耦合项双电子部分
        //sigma表示电子自旋，0为alpha，1为beta
        double G_ij(Slater_det &k1, Slater_det &k2, const std::vector<int> &Num, const std::vector<int> &index);
        double G1_ij(Slater_det &k1, Slater_det &k2, const int sigma, const std::vector<int> &Num, const std::vector<int> &index); //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
        double G2_ij(Slater_det &k1, Slater_det &k2, const std::vector<int> &Num, const std::vector<int> &index);                  //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
};



#endif