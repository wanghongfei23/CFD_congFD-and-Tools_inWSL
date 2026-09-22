#pragma once
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <limits>     // 添加这一行以支持 std::numeric_limits
#include <algorithm>  // 添加这一行以支持 std::signbit
#include <cmath>      // 添加这一行以支持 std::isnan
#include "macro.hpp"  // 添加这一行以支持 real 类型定义

// 在这里声明全局变量
// 示例：
inline real global_beta_0 = 0.0;
inline real global_beta_1 = 0.0;
inline real global_beta_2 = 0.0;

inline real chi_0 = 0.0;
inline real chi_1 = 0.0;
inline real chi_2 = 0.0;
inline real max_chi = 0.0;
inline real global_y1 = 0.0;
inline real global_y2 = 0.0;

inline real global_theta = 0.0;

inline real tau = 0.0;

// 小量保护
inline real epsilon = 1e-40;

// 注意：使用 inline 关键字确保在多个编译单元中只有一个实例

// NaN检测和提示函数
inline void check_and_warn_nan(real value, const std::string& name) {
    if (std::isnan(value)) {
        std::cerr << "Warning: " << name << " is NaN!" << std::endl;
    }
}

inline void global_write_y1y2(){
    std::string statFileName_y1 = "globalTest/y1.txt";
    std::string statFileName_y2 = "globalTest/y2.txt";

    tau = std::abs(global_beta_0-global_beta_2);
    chi_0 = tau/(global_beta_0 + epsilon);
    chi_1 = tau/(global_beta_1 + epsilon);
    chi_2 = tau/(global_beta_2 + epsilon);
    max_chi = std::max(chi_0, std::max(chi_1, chi_2));
    global_y1 = max_chi/10.0;
    global_y2 = global_y1 + 1.0;
    
    // 检查NaN并提示
    check_and_warn_nan(global_y1, "global_y1");
    check_and_warn_nan(global_y2, "global_y2");

    std::ofstream statFile_y1(statFileName_y1, std::ios::app);
    statFile_y1 << global_y1 << "\n";
    statFile_y1.close();
    std::ofstream statFile_y2(statFileName_y2, std::ios::app);
    statFile_y2 << global_y2 << "\n";
    statFile_y2.close();
}

inline void global_write_theta(){
    std::string statFileName = "globalTest/theta.txt";

    tau = std::abs(global_beta_0 - global_beta_2);
    
    chi_0 = tau/(global_beta_0 + epsilon);
    chi_1 = tau/(global_beta_1 + epsilon);
    chi_2 = tau/(global_beta_2 + epsilon);
    
    max_chi = std::max(chi_0, std::max(chi_1, chi_2));
    global_theta = 1.0 / (1.0 + max_chi / 10.0);
    
    // 检查NaN并提示
    check_and_warn_nan(global_theta, "global_theta");
    
    std::ofstream statFile(statFileName, std::ios::app);
    statFile << global_theta << "\n";
    statFile.close();
}