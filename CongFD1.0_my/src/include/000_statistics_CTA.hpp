#pragma once
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>

// 在这里声明全局变量
// 示例：
inline int global_counter_5 = 0;
inline int global_counter_6 = 0;
inline int global_counter_7 = 0;
inline int global_counter_8 = 0;
inline int global_counter_9 = 0;
inline int global_counter_10 = 0;

inline bool pandaun_001 = false;

inline bool CTA_counter_on_off = false;
// inline bool CTA_counter_on_off = true;

// 注意：使用 inline 关键字确保在多个编译单元中只有一个实例

inline void output_CTA(){
    if (global_counter_5 || global_counter_6 || global_counter_7 || global_counter_8 || global_counter_9 || global_counter_10) {
        std::string statFileName = "CT-A_counter.txt";
        
        // 检查文件是否存在，如果存在则创建带编号的文件
        int counter = 1;
        std::string newStatFileName = statFileName;
        while (std::filesystem::exists(newStatFileName)) {
            size_t dotPos = statFileName.find_last_of(".");
            if (dotPos != std::string::npos) {
                // 在扩展名前插入编号
                newStatFileName = statFileName.substr(0, dotPos) + "_" + std::to_string(counter) + statFileName.substr(dotPos);
            } else {
                // 没有扩展名的情况
                newStatFileName = statFileName + "_" + std::to_string(counter);
            }
            counter++;
        }
        
        std::ofstream statFile(newStatFileName);
        statFile << "CT Value Statistics:\n";
        statFile << "CT=1e-05 count: " << global_counter_5 << "\n";
        statFile << "CT=1e-06 count: " << global_counter_6 << "\n";
        statFile << "CT=1e-07 count: " << global_counter_7 << "\n";
        statFile << "CT=1e-08 count: " << global_counter_8 << "\n";
        statFile << "CT=1e-09 count: " << global_counter_9 << "\n";
        statFile << "CT=1e-10 count: " << global_counter_10 << "\n";
        statFile.close();
    }
    if (pandaun_001)
    {
        std::cout << "有 非-5到-10的预期值" << "\n";
    }
    else
    {
        std::cout << "无 非-5到-10的预期值" << "\n";
    }
}