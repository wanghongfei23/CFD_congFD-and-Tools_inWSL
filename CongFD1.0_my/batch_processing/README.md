# 批量处理配置文件说明

## 重要提示

配置文件必须严格按照16行格式编写，缺少任何一行都可能导致程序读取错误或运行异常。

## 配置文件命名规范

配置文件应按照以下格式命名：
```
info_"算例名称"_"插值方法".txt
```

例如：
- `info_Lax_WCNS5CONGZ.txt`
- `info_Lax_WHFTCNSA.txt`

## 配置文件内容格式

配置文件**必须**包含16行，每行一个参数，顺序如下：

1. 插值方法编号 (interMethod)
2. 时间步数 (endStep)
3. 输出时间间隔 (outputDt)
4. CFL数 (CFL)
5. 算例编号 (nCase)
6. 计算区域左边界 (calZone[0])
7. 计算区域右边界 (calZone[1])
8. 计算区域下边界 (calZone[2])
9. 计算区域上边界 (calZone[3])
10. 计算区域z方向起始 (calZone[4])
11. 计算区域z方向结束 (calZone[5])
12. X方向网格点数 (iMax[0])
13. Y方向网格点数 (iMax[1])
14. Z方向网格点数 (iMax[2])
15. 维度 (dim)
16. 线程数 (omp_threads)

## 参数说明

### 1. 插值方法编号 (interMethod)

| 编号 | 插值方法       | 说明              |
|------|----------------|-------------------|
| 0    | FIRSTORDER     | 一阶迎风格式      |
| 1    | MUSCL          | MUSCL格式         |
| 2    | WCNS5          | WENO-JS           |
| 3    | WCNSZ5         | WENO-Z            |
| 4    | WCNS5Char      | WCNS5特征重构     |
| 5    | WCNSZ5Char     | WCNSZ5特征重构    |
| 6    | WCNS5CONG      | WCNS5改进格式     |
| 7    | TCNSCongA      | TENO-A格式        |
| 8    | WCNS5CONGZ     | TENO-S          |
| 9    | WCNS5CONGZCT4  | TENO-S CT4版本  |
| 10   | WCNS5CONGZCT7  | TENO-S CT7版本  |
| 11   | TCNS5          | TENO            |
| 12   | TCNS5CT4       | TENO CT4版本    |
| 13   | TCNS5CT7       | TENO CT7版本    |
| 14   | LINEAR5        | 五阶线性格式      |
| 15   | NICEST5        | NICEST格式        |
| 16   | INTERMAX       | 最大值插值格式    |
| 17   | WHFTCNSA       | TENO-A        |
| 18   | WHFTCNSASF002   | TENO-AS-myF002   |
| 19   | WHFTCNSAH002   | TENO-AS-myH002   |
| 20   | WHFTCNSASF102   | TENO-AS-myF102   |
| 21   | WHFTCNSASF103   | TENO-AS-myF103   |
| 22   | WHFTCNSASF102_reciprocal | TENO-AS-myF102_reciprocal |
| 23   | WHFTCNSASF103_reciprocal | TENO-AS-myF103_reciprocal |
| 24   | WHFTCNSAS_fx    | TENO-AS-fx       |
| 25   | WHFTCNSAS_initial | TENO-AS-initial  |
| 26   | WHFTCNSAS_approx_1 | TENO-AS-approx_1 |
| 27   | WHFTCNSAS_fx_real | TENO-AS-fx-real  |
| 28   | WHFTCNSAS_approx_2 | TENO-AS-approx_2 |
| 29   | WHFTCNSASF202_2S | TENO-AS-myF202_2S |
| 30   | WHFTCNSASF202 | TENO-AS-myF202 |
| 31   | WHFTCNSASFf_5_10 | TENO-AS-Ff_5_10 |

### 2. 算例编号 (nCase)

#### 一维算例

| 编号 | 算例名称           | 说明             |
|------|--------------------|------------------|
| 0    | Sod                | Sod激波管        |
| 1    | ShuOsher           | Shu-Osher问题    |
| 2    | Lax                | Lax激波管        |
| 3    | sedov              | Sedov爆炸波      |
| 4    | Woodward_Colella   | Woodward-Colella问题 |
| 5    | Double_sparse_wave | 双稀疏波问题     |

#### 二维算例

**参考文献**: Schulz-Rinne, C. W., Collins, J. P., & Glaz, H. M. (2004). 
"Numerical Solution of the Two-Dimensional Riemann Problems for Gas Dynamics". 
SIAM Journal on Scientific Computing, 25(6), 1895-1917.

| nCase | 算例名称              | 文献对应                    | 说明                      |
|-------|-----------------------|----------------------------|---------------------------|
| 0     | 2D_Riemann_Config3    | Configuration 3            | 双激波+接触间断           |
| 1     | 2D_Riemann_Config4    | Configuration 4 (Vortex)   | 涡旋结构                  |
| 2     | implosion             | -                          | 内爆问题                  |
| 3     | RTI                   | -                          | 瑞利-泰勒不稳定性         |
| 4     | Double_Mach           | -                          | 双马赫反射                |
| 5     | 2D_Riemann_Config12   | Configuration 12           | 混合波系                  |
| 6     | KHI                   | -                          | 开尔文-亥姆霍兹不稳定性   |
| 7     | 2D_Riemann_Config6    | Configuration 6            | 反向涡旋结构              |

**注意**: 
- Configuration 编号对应 Schulz-Rinne et al. (2004) 论文中的标准配置
- Config 4 和 Config 6 的主要区别在于水平速度场方向相反
- Config 3 和 Config 12 都是混合波系但初始条件不同

## 一维与二维问题参数设置差异

虽然一维和二维问题使用相同的16行配置文件格式，但在参数设置上有以下区别：

### 一维问题 (dim=1)
- calZone[2] 到 calZone[5] 应设置为 0
- iMax[1] 和 iMax[2] 通常设置为 2
- 示例：
  ```
  8        # interMethod
  100      # endStep
  0.01     # outputDt
  0.5      # CFL
  2        # nCase (Lax)
  -0.5     # calZone[0]
  0.5      # calZone[1]
  0        # calZone[2] (unused)
  0        # calZone[3] (unused)
  0        # calZone[4] (unused)
  0        # calZone[5] (unused)
  201      # iMax[0]
  2        # iMax[1] (unused)
  2        # iMax[2] (unused)
  1        # dim
  4        # omp_threads
  ```

### 二维问题 (dim=2)
- calZone[0] 到 calZone[3] 分别表示 x 和 y 方向的边界
- iMax[0] 和 iMax[1] 分别表示 x 和 y 方向的网格点数
- 示例：
  ```
  2        # interMethod
  100      # endStep
  0.01     # outputDt
  0.5      # CFL
  0        # nCase (2D_Riemann_1)
  -0.5     # calZone[0] (x_min)
  0.5      # calZone[1] (x_max)
  -0.5     # calZone[2] (y_min)
  0.5      # calZone[3] (y_max)
  0        # calZone[4] (unused)
  0        # calZone[5] (unused)
  201      # iMax[0] (x方向网格点数)
  201      # iMax[1] (y方向网格点数)
  2        # iMax[2] (unused)
  2        # dim
  4        # omp_threads
  ```

## 批量处理脚本使用方法

### run_all_cases.sh 脚本

该脚本会自动运行当前目录下所有以 `info_` 开头、以 `.txt` 结尾的配置文件。

使用步骤：
1. 确保所有配置文件按照命名规范放置在 `batch_processing` 目录中
2. 在终端中运行脚本：
   ```bash
   cd batch_processing
   ./run_all_cases.sh
   ```
3. 脚本会依次执行每个算例，并在控制台显示执行状态

注意事项：
- 脚本会将当前配置文件复制为 `info.txt` 供程序读取
- 程序执行结果将保存在 `batch_processing_data` 目录中
- 确保在运行脚本前已正确编译主程序 (`../build/congFD`)