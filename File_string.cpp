#include "File_string.h"

#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::string;
using std::vector;


int COUNT_DEBUG = 0;
bool read_int(const string &file_name, int &nelec, int &nOrb, double &MS, double &h_nuc, vector<double> &h, vector<double> &g)//读取文件中的积分，电子数与分子轨道数
{
    std::ifstream file;
    file.open(file_name.c_str());

    if(!file.is_open()){
        cout << "未成功打开积分文件" << endl;
        return 0;
    }

    string line;
    string delim;
    
    //第一行格式如右，用于读取轨道个数与电子个数:"&FCI NORB=114,NELEC=42,MS2=0,"
    getline(file, line);
    delim = "=";
    vector<string> v = split(line, delim);

    delim = ",";
    vector<string> v1 = split(v[1], delim);
    nOrb = atof(v1[0].c_str());
    v1 = split(v[2], delim);
    nelec = atof(v1[0].c_str());
    v1 = split(v[3], delim);
    MS = atof(v1[0].c_str()) ;
    /*
    cout << nelec << endl
         << nOrb << endl;
    */

    //右为读取开始的标志:" &END"
    delim = "&END";

    h.resize(nOrb * nOrb);
    g.resize(nOrb * nOrb * nOrb * nOrb);

    bool flag=0;
    while(getline(file,line))
    {
        if(!flag)//用于判断是否可以开始读取积分,flag=1表示可以开始读取积分
            flag = judge_start(line, delim);
        else
        {
            flag = read_line_int(line, h_nuc, h, g);
            if(!flag)//读到空行，读取结束，即flag又返回0
                break;
        }
    }
    cout << "COUNT_DEBUG: " << COUNT_DEBUG << endl;
    file.close();

    return 1;
}

bool read_line_int(const string &line, double &h_nuc, vector<double> &h, vector<double> &g)//读取轨道积分
{
    string delim=" ";
    vector<string> v = split(line, delim);

    double INT = string2int(v[0]);

    int i = atof(v[1].c_str()) - 1;
    int j = atof(v[2].c_str()) - 1;
    int k = atof(v[3].c_str()) - 1;
    int l = atof(v[4].c_str()) - 1;

    int dim2 = h.size();
    int dim = sqrt(dim2); //数组的维度
    int dim3 = dim * dim2;

    //cout << i << " " << j << " " << k << " " << l << "\t" << INT << std::endl;
    if(k!=-1)//双电子积分
    {
    //用来判断是否读到了已修改过的位置，以防万一
    
        if(fabs(g[i * dim3 + j * dim2 + k * dim + l]) > 0)              //ijkl
            cout << "bad" << endl;
        if(fabs(g[i * dim3 + j * dim2 + l * dim + k]) > 0)              //jikl
            cout << "bad" << endl;
        if(fabs(g[j * dim3 + i * dim2 + k * dim + l]) > 0)              //ijlk
            cout << "bad" << endl;
        if(fabs(g[j * dim3 + i * dim2 + l * dim + k]) > 0)              //jilk
            cout << "bad" << endl;
    
        if(fabs(g[k * dim3 + l * dim2 + i * dim + j]) > 0)              //klij
            cout << "bad" << endl;
        if(fabs(g[k * dim3 + l * dim2 + j * dim + i]) > 0)              //klji
            cout << "bad" << endl;
        if(fabs(g[l * dim3 + k * dim2 + i * dim + j]) > 0)              //lkij
            cout << "bad" << endl;
        if(fabs(g[l * dim3 + k * dim2 + j * dim + i]) > 0)              //lkji
            cout << "bad" << endl;
        
        //分割线
            
        if(g[i * dim3 + j * dim2 + k * dim + l] == 0)               //ijkl
            COUNT_DEBUG++;
        g[i * dim3 + j * dim2 + k * dim + l] = INT;                 //ijkl

        if(g[i * dim3 + j * dim2 + l * dim + k] == 0)               //jikl
            COUNT_DEBUG++;
        g[i * dim3 + j * dim2 + l * dim + k] = INT;                 //jikl

        if(g[j * dim3 + i * dim2 + k * dim + l] == 0)               //ijlk
            COUNT_DEBUG++;
        g[j * dim3 + i * dim2 + k * dim + l] = INT;                 //ijlk

        if(g[j * dim3 + i * dim2 + l * dim + k] == 0)               //jilk
            COUNT_DEBUG++;
        g[j * dim3 + i * dim2 + l * dim + k] = INT;                 //jilk
    

    
        if(g[k * dim3 + l * dim2 + i * dim + j] == 0)               //klij
            COUNT_DEBUG++;
        g[k * dim3 + l * dim2 + i * dim + j] = INT;                 //klij

        if(g[k * dim3 + l * dim2 + j * dim + i] == 0)               //klji
            COUNT_DEBUG++;
        g[k * dim3 + l * dim2 + j * dim + i] = INT;                 //klji

        if(g[l * dim3 + k * dim2 + i * dim + j] == 0)               //lkij
            COUNT_DEBUG++;
        g[l * dim3 + k * dim2 + i * dim + j] = INT;                 //lkij

        if(g[l * dim3 + k * dim2 + j * dim + i] == 0)               //lkji
            COUNT_DEBUG++;
        g[l * dim3 + k * dim2 + j * dim + i] = INT;                 //lkji
        //cout << i + 1 << " " << j + 1 << " " << k + 1 << " " << l + 1 << "\t" << g[i * dim3 + j * dim2 + k * dim + l] << std::endl;
        
    }
    else if (i != -1)//单电子积分
    {
        h[i * dim + j] = INT;
        h[j * dim + i] = INT;
        //cout << i + 1 << " " << j + 1 << "\t" << h[i * dim + j] << std::endl;
    }
    else//核积分项
    {
        h_nuc = INT;
        //cout << h_nuc << endl;
        return 0;//表示读取结束，读到文件末尾
    }

    return 1;//表示读取正常，读取仍未结束
}

double string2int(const string &s)//将字符串读成数字，形式x.xxxxxExxx
{
    string delim="E";
    vector<string> v=split(s,delim);
    int n=v.size();
    double tmp=0;

    if(n==1)
    {
        tmp=atof(v[0].c_str());
    }
    else
    {
        double t1=atof(v[0].c_str());
        int t2=atof(v[1].c_str());
        tmp=t1*pow(10,t2);
    }
    return tmp;

}

bool judge_start(const string &line, const string &delim)//判断是否可以开始读取积分
{
    bool flag = 0;
    vector<string> v1 = split(line, delim);
    if(v1.size()==0)
    {
        flag = 1;
    }
    else{
        string delim2 = " ";
        vector<string> v2 = split(v1[0], delim2);
        if(v2.size()==0)
            flag = 1;
    }

    return flag;
}

vector<string> split(const string& s, const string& delim)
{
    vector<string> v;
    string::size_type pos1, pos2;
    pos2 = s.find(delim);
    pos1 = 0;
    while (string::npos != pos2)
    {
        if(pos1!=pos2)//防止出现连续的分隔符delim
            v.push_back(s.substr(pos1, pos2 - pos1)); //开始位置与获取字符串的长度，取出来放到vector中

        pos1 = pos2 + delim.size();
        pos2 = s.find(delim, pos1);

    }
        
    if(pos1 != s.length())
        v.push_back(s.substr(pos1));
    
    return v;
}

std::ostream &operator<<(std::ostream &os, const vector<double> &Matrix)
{
    int ndim = sqrt(Matrix.size());

    os << std::endl;
    for (int i = 0; i < ndim; i++)
    {
        for (int j = 0; j < ndim;j++)
        {
            os << Matrix[i * ndim + j];
            if(Matrix[i * ndim + j]<0)
                os << "       ";
            else
                os << "        ";
        }
        os << std::endl;
    }
    return os;
};
