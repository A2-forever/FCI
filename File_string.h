#ifndef FILE_STRING_H_
#define FILE_STRING_H_

#include <iostream>
#include <string>
#include <vector>

bool read_int(const std::string &file_name, int &nelec, int &n_Orb, double &MS, double &h_nuc, std::vector<double> &h, std::vector<double> &g); //读取文件中的积分，电子数与分子轨道数
bool read_line_int(const std::string &line, double &h_nuc, std::vector<double> &h, std::vector<double> &g);//读取轨道积分

double string2int(const std::string &s);//将字符串读成数字，形式x.xxxxxExxx
bool judge_start(const std::string &line, const std::string &delim);//判断是否可以开始读取积分
std::vector<std::string> split(const std::string &s, const std::string &delim);
std::ostream &operator<<(std::ostream &os, const std::vector<double> &Matrix);

#endif