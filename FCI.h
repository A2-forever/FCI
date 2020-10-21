#ifndef FCI_H_
#define FCI_H_

#include "CI.h"
#include <vector>


class FCI{
    private:
        int nOrb = 0;          //分子轨道数
        double h_nuc = 0;      //核积分项
        std::vector<double> h; //单电子积分
        std::vector<double> g; //双电子积分

        int dim = 0;
        int dim2 = 0;
        int dim3 = 0;

    public:
        FCI(); //默认构造函数
        FCI(const double &h_nuc_ex, const std::vector<double> &h_ex, const std::vector<double> &g_ex, const int nOrb_ex);
        ~FCI(); //默认析构函数
        double get(const int &i, const int &j, const int &k, const int &l);

        //计算分子轨道i与分子轨道j关于hamilton算符的耦合项
        double H_ij(CSF &k1, CSF &k2);
        double H_ij(Slater_det &k1, Slater_det &k2);

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