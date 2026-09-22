#pragma once
#include "macro.hpp"
#include "000_statistics_CTA.hpp"
#include "100_LAD_for_beta0_1_2.hpp"
#include <array>

constexpr real u1(real q1, real q2, real q3)
{
    return 3.0 / 8.0 * q1 - 5.0 / 4.0 * q2 + 15.0 / 8.0 * q3;
}
constexpr real u2(real q1, real q2, real q3)
{
    return -1.0 / 8.0 * q1 + 3.0 / 4.0 * q2 + 3.0 / 8.0 * q3;
}
constexpr real u3(real q1, real q2, real q3)
{
    return 3.0 / 8.0 * q1 + 3.0 / 4.0 * q2 - 1.0 / 8.0 * q3;
}

inline real weno5_JSchen(std::array<real, 5> q)
{
    real eps = 1e-6;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

    real sumbeta = 0, result = 0;
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] / pow(eps + beta[i], 2);
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}

inline real Nicest5(std::array<real, 5> q)
{
    constexpr real eps = 1e-40;

    // 计算 beta 值
    const real beta_31 = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2)
        + eps;

    const real beta_32 = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2)
        + eps;

    const real beta_33 = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
        + eps;

    const real beta_51 = 0.06 * (std::pow(1. * q[0] - 4. * q[1] + 6. * q[2] - 4. * q[3] + 1. * q[4], 2) + std::pow(-0.5 * q[0] + 1. * q[1] + 0. * q[2] - 1. * q[3] + 0.5 * q[4], 2) + eps);

    // Stage 1: 找到最小的 beta 值
    const real min_beta_s1 = std::min({ beta_31, beta_32, beta_33, beta_51 });

    // 如果最小 beta 值对应 m=5，直接计算并返回 u_51
    if (min_beta_s1 == beta_51) {
        return 3.0 / 128.0 * q[0]
            - 5.0 / 32.0 * q[1]
            + 45.0 / 64.0 * q[2]
            + 15.0 / 32.0 * q[3]
            - 5.0 / 128.0 * q[4];
    }

    // Stage 2: 计算 beta_41 和 beta_42
    const real beta_41 = 0.06 * (std::pow(0 * q[0] + 1. * q[1] - 2. * q[2] + 1. * q[3], 2) + std::pow(-1. * q[0] + 3. * q[1] - 3. * q[2] + 1. * q[3], 2) + eps);

    const real beta_42 = 0.06 * (std::pow(0 * q[1] + 1. * q[2] - 2. * q[3] + 1. * q[4], 2) + std::pow(-1. * q[1] + 3. * q[2] - 3. * q[3] + 1. * q[4], 2) + eps);

    // 找到 Stage 2 中的最小 beta 值
    const real min_beta_s2 = std::min({ beta_41, beta_42, min_beta_s1 });

    // 根据最小 beta 值返回对应的 u 值
    if (min_beta_s2 == beta_41) {
        return 1.0 / 16.0 * q[0]
            - 5.0 / 16.0 * q[1]
            + 15.0 / 16.0 * q[2]
            + 5.0 / 16.0 * q[3];
    }

    if (min_beta_s2 == beta_42) {
        return -1.0 / 16.0 * q[1]
            + 9.0 / 16.0 * q[2]
            + 9.0 / 16.0 * q[3]
            - 1.0 / 16.0 * q[4];
    }

    if (min_beta_s2 == beta_31) {
        return 3.0 / 8.0 * q[0]
            - 5.0 / 4.0 * q[1]
            + 15.0 / 8.0 * q[2];
    }

    if (min_beta_s2 == beta_32) {
        return -1.0 / 8.0 * q[1]
            + 3.0 / 4.0 * q[2]
            + 3.0 / 8.0 * q[3];
    }

    if (min_beta_s2 == beta_33) {
        return 3.0 / 8.0 * q[2]
            + 3.0 / 4.0 * q[3]
            - 1.0 / 8.0 * q[4];
    }

    // 默认返回（理论上不会执行）
    return 0.0;
}

// inline real weno5_JSchen(std::array<real, 5> q)
// {
//     real eps = 1e-6;
//     std::array<real, 6> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0, 0.06, 0.06, 0.06 };
//     std::array<real, 6> beta;
//     std::array<real, 6> u;

//     // Third-order candidates
//     u[0] = (3.0 * q[0] - 10.0 * q[1] + 15.0 * q[2]) / 8.0;
//     u[1] = (-1.0 * q[1] + 6.0 * q[2] + 3.0 * q[3]) / 8.0;
//     u[2] = (3.0 * q[2] + 6.0 * q[3] - 1.0 * q[4]) / 8.0;

//     // Fourth-order candidates
//     u[3] = (u[0] + 5.0 * u[1]) / 6.0;
//     u[4] = (u[1] + u[2]) / 2.0;

//     // Fifth-order candidate
//     u[5] = (3.0 * u[3] + 5.0 * u[4]) / 8.0;

//     // Calculate beta values
//     beta[0] = pow(q[0] - 2.0 * q[1] + q[2], 2) + 0.25 * pow(q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
//     beta[1] = pow(q[1] - 2.0 * q[2] + q[3], 2) + 0.25 * pow(q[1] - q[3], 2);
//     beta[2] = pow(q[2] - 2.0 * q[3] + q[4], 2) + 0.25 * pow(3.0 * q[2] - 4.0 * q[3] + q[4], 2);

//     beta[3] = pow(u[0] - 2.0 * u[1] + u[2], 2) + 0.25 * pow(u[0] - 4.0 * u[1] + 3.0 * u[2], 2);
//     beta[4] = pow(u[1] - 2.0 * u[2] + u[3], 2) + 0.25 * pow(u[1] - u[3], 2);
//     beta[5] = pow(u[2] - 2.0 * u[3] + u[4], 2) + 0.25 * pow(3.0 * u[2] - 4.0 * u[3] + u[4], 2);

//     // Scale beta values
//     for (int i = 0; i < 6; i++) {
//         beta[i] = (beta[i] + eps) * gamma[i];
//     }

//     // Two-stage process to find the minimum beta
//     int minIndex1 = 0;
//     real minValue1 = beta[0];
//     for (int i = 1; i < 4; i++) {
//         if (beta[i] < minValue1) {
//             minValue1 = beta[i];
//             minIndex1 = i;
//         }
//     }

//     if (minIndex1 == 3) {
//         return u[5];
//     }

//     int minIndex2 = minIndex1;
//     real minValue2 = minValue1;
//     for (int i = 4; i < 6; i++) {
//         if (beta[i] < minValue2) {
//             minValue2 = beta[i];
//             minIndex2 = i;
//         }
//     }

//     return u[minIndex2];
// }

constexpr real weno5_Cong(std::array<real, 5> q)
{
    real eps = 1e-14;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta, u;
    // beta[0]= 0.0/1.0 *pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2)
    //          + 1.0/1.0 *pow(1.0*q[0]-4.0*q[1]+3.0*q[2],2);

    //  beta[1]= 0.0/1.0  *pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2)
    //          + 1.0/1.0 *pow(1.0*q[1]+0.0*q[2]-1.0*q[3],2);

    //  beta[2]= 0.0/1.0 *pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)
    //          + 1.0/1.0*pow(3.0*q[2]-4.0*q[3]+1.0*q[4],2);

    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) + 10 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2) / (q[2] * q[2]);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) + 10 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2) / (q[2] * q[2]);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) / (q[2] * q[2]) + 10 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2) / (q[2] * q[2]);

    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    real sumbeta = 0, result = 0;
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] / pow(eps + beta[i], 2.0);
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}
constexpr real weno5_Z(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

    real sumbeta = 0, result = 0;
    real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] * (C + pow(tau / (beta[i] + eps), qq));
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}

constexpr real Teno5_ZCT4(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    real sumbeta = 0;
    real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        real tempp = C + tau / (beta[i] + eps);
        tempp *= tempp;
        beta[i] = tempp * tempp * tempp;
        sumbeta += beta[i];
    }
    real CT = 1e-4 * sumbeta;
    // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
    unsigned short flag = 0;
    if (beta[0] < CT)
        flag += 1;
    if (beta[1] < CT)
        flag += 2;
    if (beta[2] < CT)
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

constexpr real Teno5_ZCT7(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    real sumbeta = 0;
    real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        real tempp = C + tau / (beta[i] + eps);
        tempp *= tempp;
        beta[i] = tempp * tempp * tempp;
        sumbeta += beta[i];
    }
    real CT = 1e-7 * sumbeta;
    // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
    unsigned short flag = 0;
    if (beta[0] < CT)
        flag += 1;
    if (beta[1] < CT)
        flag += 2;
    if (beta[2] < CT)
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

// x 的平方（替代 std::pow(x, 2)：后者在 -O3 下仍是运行时库调用）
constexpr real pow2(real x) { return x * x; }

constexpr real Teno5_Z(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta;
    beta[0] = 1.0/1.0 * pow2(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2]) +
              1.0/4.0 * pow2(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2]);
    beta[1] = 1.0/1.0 * pow2(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3]) + 
              1.0/4.0 * pow2(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3]);
    beta[2] = 1.0/1.0 * pow2(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4]) + 
              1.0/4.0 * pow2(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4]);

    real sumbeta = 0;
    real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        real tempp = C + tau / (beta[i] + eps);
        tempp *= tempp;
        beta[i] = tempp * tempp * tempp;
        sumbeta += beta[i];
    }
    real CT = 1e-5 * sumbeta;
    // volatile unsigned flag=(beta[0]<CT)+((beta[1]<CT)<<1)+((beta[2]<CT)<<2);
    unsigned short flag = 0;
    if (beta[0] < CT)
        flag += 1;
    if (beta[1] < CT)
        flag += 2;
    if (beta[2] < CT)
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0/128.0*q[0] - 5.0/32.0*q[1] + 45.0/64.0*q[2] + 15.0/32.0*q[3] - 5.0/128.0*q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0/16.0*q[1] + 9.0/16.0*q[2] + 9.0/16.0*q[3] - 1.0/16.0*q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0/8.0*q[2] + 3.0/4.0*q[3] - 1.0/8.0*q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0/8.0*q[2] + 3.0/4.0*q[3] - 1.0/8.0*q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0/16.0*q[0] - 5.0/16.0*q[1] + 15.0/16.0*q[2] + 5.0/16.0*q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0/8.0*q[1] + 3.0/4.0*q[2] + 3.0/8.0*q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0/8.0*q[0] - 5.0/4.0*q[1] + 15.0/8.0*q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }

}



constexpr real Teno5_ZConvex(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};
    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    real sumbeta = 0, result = 0;
    real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        real tempp = C + tau / (beta[i] + eps);
        tempp *= tempp;
        beta[i] = tempp * tempp * tempp;
        sumbeta += beta[i];
    }

    real CT = 1e-5 * sumbeta, sumGamma = 0;
    for (int i = 0; i < 3; i++) {
        if (beta[i] > CT) {
            // result+=gamma[i]*(*u[i])(q[i],q[i+1],q[i+2]);
            result += gamma[i] * u[i];
            sumGamma += gamma[i];
        }
    }
    result /= sumGamma;
    return result;
}

const static real CTi = pow(1.5 * 1e-10, 1.0 / 6.0);

constexpr real Teno5_CongZ(std::array<real, 5> q)
{
    real eps = 1e-40; // 1e-10;
    std::array<real, 3> beta ;
    beta[0] = 1.0/1.0 * pow2(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2]) +
              1.0/4.0 * pow2(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2]);
    beta[1] = 1.0/1.0 * pow2(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3]) + 
              1.0/4.0 * pow2(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3]);
    beta[2] = 1.0/1.0 * pow2(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4]) + 
              1.0/4.0 * pow2(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4]);

    // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1: 2):((beta[2]>beta[0])? 0 : 2);
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // constexpr real CT=0.23050581003334941;//4
    // constexpr real CT = 0.15704178024750198; // 5
    constexpr real CT = 0.15704178024750198; // 5  std::pow(1.5 * 1e-5, 1.0 / 6.0) 预计算值
    // constexpr real CT=0.08;//5
    // constexpr real CT=0.10699131939336631;//6
    // constexpr real CT=0.072892337360747711; //7
    // constexpr real CT=0.033833625914958219;//9
    // constexpr real CT=0.023050581003334944;//10
    constexpr real CT_1 = 1 - CT;
    real tau = std::abs(beta[2] - beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];
    // unsigned flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

constexpr real Teno5_CongZCT4(std::array<real, 5> q)
{
    real eps = 1e-40; // 1e-10;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
            + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
            + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
            + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };

    // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1: 2):((beta[2]>beta[0])? 0 : 2);
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    constexpr real CT = 0.23050581003334941; // 4
    // constexpr real CT=0.15704178024750198;//5
    // constexpr real CT=0.157;//5
    // constexpr real CT=0.10699131939336631;//6
    // constexpr real CT=0.072892337360747711; //7
    // constexpr real CT=0.033833625914958219;//9
    // constexpr real CT=0.023050581003334944;//10
    constexpr real CT_1 = 1 - CT;
    real tau = std::abs(beta[2] - beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];
    // unsigned flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}
constexpr real Teno5_CongZCT7(std::array<real, 5> q)
{
    real eps = 1e-40; // 1e-10;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
            + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
            + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
            + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };

    // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1: 2):((beta[2]>beta[0])? 0 : 2);
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // constexpr real CT=0.23050581003334941;//4
    // constexpr real CT=0.15704178024750198;//5
    // constexpr real CT=0.157;//5
    // constexpr real CT=0.10699131939336631;//6
    constexpr real CT = 0.072892337360747711; // 7
    // constexpr real CT=0.033833625914958219;//9
    //  constexpr real CT=0.023050581003334944;//10
    constexpr real CT_1 = 1 - CT;
    real tau = std::abs(beta[2] - beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];
    // unsigned flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}

constexpr real Teno5_CongA(std::array<real, 5> q)
{
    real eps = 1e-40; // 1e-10;
    std::array<real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
            + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
            + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
            + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };

    // std::array<real,3> beta={
    //         pow(1.0*q[0]-2.0*q[1]+1.0*q[2],2),
    //         pow(1.0*q[1]-2.0*q[2]+1.0*q[3],2),
    //         pow(1.0*q[2]-2.0*q[3]+1.0*q[4],2)};

    // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1: 2):((beta[2]>beta[0])? 0 : 2);
    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // real tau=std::abs(beta[2]+beta[0]-2*beta[1]);
    real tau = std::abs(beta[2] - beta[0]);

    // constexpr real CT=0.23050581003334941;//4
    // constexpr real CT=0.15704178024750198;//5
    constexpr real CT = 0.1; // 6
    // constexpr real CT=0.072892337360747711;//7
    // constexpr real CT=0.033833625914958219;//9
    //  constexpr real CT=0.023050581003334944;//10
    //  real CT=0.15704178024750198;//5
    real CT_1 = 1 - CT;

    // real rr=CT/1/(beta[minBeta]+eps)-CT_1/(tau+eps);//CT*tau-CT_1*beta[minBeta];
    real mulbeta = beta[0] * beta[1] * beta[2];
    real rr = CT * tau * (beta[0] * beta[1] + beta[1] * beta[2] + beta[0] * beta[2]) - CT_1 * mulbeta; // CT*tau-CT_1*beta[minBeta];
    real ll = tau * mulbeta; // tau*beta[minBeta];
    // unsigned flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
    unsigned short flag = 0;
    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;
    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    default:
        /* 0,0,0 */
        return q[2];
        break;
    }
}
constexpr real Teno5_CongC(std::array<real, 5> q)
{
    real eps = 1e-40; // 1e-10;
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // real tau=std::abs(beta[2]+beta[0]-2*beta[1]);

    real tau = std::abs(beta[2] - beta[0]);

    // constexpr real CT=0.23050581003334941;//4
    // constexpr real CT=0.15704178024750198;//5
    // constexpr real CT=0.10699131939336631;//6
    // constexpr real CT=0.072892337360747711; //7
    // constexpr real CT=0.033833625914958219;//9
    //  constexpr real CT=0.023050581003334944;//10
    real CT = 0.15704178024750198; // 5
    real CT_1 = 1 - CT;

    real rr = CT * tau - CT_1 * beta[minBeta];
    real ll = tau * beta[minBeta];

    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    // std::array<real(*)(real,real,real),3> u={&u1,&u2,&u3};

    real sumbeta = 0, result = 0, sumGamma = 0;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    for (int i = 0; i < 3; i++) {
        if (ll >= rr * beta[i]) {
            sumGamma += gamma[i];
            result += gamma[i] * u[i];
        }
    }
    return result / sumGamma;
}

constexpr real Teno5_Cong(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    real sumbeta = 0, result = 0;

    real minBeta = *std::min_element(beta.begin(), beta.end());

    std::array<real (*)(real, real, real), 3> u = { &u1, &u2, &u3 };
    real CT = 3.5, sumGamma = 0;
    for (int i = 0; i < 3; i++) {
        if (beta[i] < CT * (minBeta + eps)) {
            sumGamma += gamma[i];
            result += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
        }
    }
    result /= sumGamma;
    return result;
}

constexpr real Teno5_Cong2(std::array<real, 5> q)
{
    real eps = 1e-40;
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    real sumbeta = 0, result = 0;

    real minBeta = *std::min_element(beta.begin(), beta.end());
    if (minBeta < eps)
        return q[2];
    real CT = 2, sumGamma = 0;
    int flag = 0, flags[3] = { 1, 2, 4 };

    for (int i = 0; i < 3; i++) {
        if (beta[i] < CT * (minBeta + eps)) {
            flag += flags[i];
        }
    }
    switch (flag) {
    case 1:
        /* 1,0,0 */
        result = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    case 2:
        /* 0,1,0 */
        result = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 3:
        /* 1,1,0 */
        result = 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 4:
    case 5:
        /* 0,0,1 */
        /* 1,0,1 */
        result = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 6:
        /* 0,1,1 */
        result = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 7:
        /* 1,1,1 */
        result = 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    default:
        /* 0,0,0 */
        result = q[2];
        break;
    }

    return result;
}

constexpr std::array<real, 3> Teno5_BVDCong(std::array<real, 5> q, bool& flag)
{
    real eps = 1e-40;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    std::array<real, 3> result = { 0, 0 };
    real sumbeta = 0;
    real C = 1, qq = 6, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        beta[i] = (C + pow(tau / (beta[i] + eps), qq));
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        beta[i] /= sumbeta;

    std::array<real (*)(real, real, real), 3> u = { &u1, &u2, &u3 };
    real CT1 = 1e-7, CT2 = 1e-3, sumGamma = 0;
    for (int i = 0; i < 3; i++) {
        if (beta[i] > CT1) {
            sumGamma += gamma[i];
            result[0] += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
        }
    }
    result[0] /= sumGamma;

    real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
    if (!(extremPoint >= 1 || extremPoint <= 0)) {
        result[1] = q[2];
        result[2] = result[0];
        flag = false;
    } else {
        flag = true;
        sumGamma = 0;
        for (int i = 0; i < 3; i++) {
            if (beta[i] > CT2) {
                sumGamma += gamma[i];
                result[1] += gamma[i] * (*u[i])(q[i], q[i + 1], q[i + 2]);
            }
        }
        result[1] /= sumGamma;
        // auto minBetai=std::max_element(beta.begin(),beta.end())-beta.begin();
        // result[1]=(*u[minBetai])(q[minBetai],q[minBetai+1],q[minBetai+2]);
        // {result[2]=q[2];}
    }

    return result;
}

constexpr std::array<real, 3> Teno5_BVDMR(std::array<real, 5> q, bool& flag)
{
    real eps = 1e-6;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    std::array<real, 3> result = { 0, 0 };
    real sumbeta = 0;

    auto minBetap = std::min_element(beta.begin(), beta.end());
    auto minBetai = minBetap - beta.begin();
    auto minBeta = *minBetap;
    if (minBeta < eps)
        return { q[2], q[2], q[2] };

    real CT1 = 20, CT2 = 1.5, sumGamma = 0;
    int flag1 = 0, flag2 = 0, flags[] = { 1, 2, 4 };
    for (int i = 0; i < 3; i++) {
        if (beta[i] < CT1 * (minBeta + eps)) {
            flag1 += flags[i];
        }
        if (beta[i] < CT2 * (minBeta + eps)) //(extremPoint>=2-minBetai || extremPoint<=1-minBetai))
        {
            flag2 += flags[i];
        }
    }
    switch (flag1) {
    case 1:
        /* 1,0,0 */
        result[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
        break;
    case 2:
        /* 0,1,0 */
        result[0] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        break;
    case 3:
        /* 1,1,0 */
        result[0] = 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
        break;
    case 4:
    case 5:
        /* 0,0,1 */
        /* 1,0,1 */
        result[0] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
        break;
    case 6:
        /* 0,1,1 */
        result[0] = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 7:
        /* 1,1,1 */
        result[0] = 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
        break;
    default:
        /* 0,0,0 */
        result[0] = q[2];
        break;
    }
    switch (flag2) {
    case 1:
        /* 1,0,0 */
        result[1] = -5.0 / 12.0 * q[0] + 1.0 / 3.0 * q[1] + 13.0 / 12.0 * q[2];
        break;
    case 2:
        /* 0,1,0 */
        result[1] = 1.0 / 12.0 * q[1] + 1.0 / 3.0 * q[2] + 7.0 / 12.0 * q[3];
        break;
    case 3:
        /* 1,1,0 */
        result[1] = -9.0 / 80.0 * q[0] - 17.0 / 80.0 * q[1] + 33.0 / 80.0 * q[2] + 39.0 / 80.0 * q[3];
        break;
    case 4:
    case 5:
        /* 0,0,1 */
        /* 1,0,1 */
        result[1] = 7.0 / 12.0 * q[2] + 1.0 / 3.0 * q[3] + 1.0 / 12.0 * q[4];
        break;
    case 6:
        /* 0,1,1 */
        result[1] = -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
        break;
    case 7:
        /* 1,1,1 */
        result[1] = -3.0 / 160.0 * q[0] + 1.0 / 80.0 * q[1] + 9.0 / 20.0 * q[2] + 51.0 / 80.0 * q[3] - 13.0 / 160.0 * q[4];
        break;
    default:
        /* 0,0,0 */
        result[1] = q[2];
        break;
    }
    result[1] = result[0];

    return result;
}

constexpr real Teno5_CongSort(std::array<real, 5> q)
{
    real eps = 1e-12;
    std::array<real (*)(real, real, real), 3> u = { &u1, &u2, &u3 };

    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    real CC = (beta[1] * 2 - beta[2] - beta[0]) / 2;
    // beta[0]+=2.0/3.0*CC;
    // beta[1]-=1.0/3.0*CC;
    // beta[2]+=2.0/3.0*CC;

    // 排序
    std::array<real, 3> index = { 0, 1, 2 };
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        });
    // if(beta[index[2]]<eps) return q[2];
    int ii = index[0];
    real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];

    real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
    real extremPoint2 = 0.5 + (q[2] - q[3]) / (q[2] - 2 * q[3] + q[4] + eps);

    bool flags[3] = { true //(u1<*std::max_element(q.begin(),q.end())&&(u1>*std::min_element(q.begin(),q.end())))
        ,
        (extremPoint >= 1 || extremPoint <= 0), (extremPoint2 >= 1 || extremPoint2 <= 0) }; //;|| (factor>6);
    // bool flag=flags[1]||flags[2];//;|| (factor>6);
    real CT2, CT = 3, sumGamma = gamma[index[0]];
    // if(!flag)
    // {
    //     CT2=3;
    //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
    // }
    // else
    // {
    //     CT2=3;
    //     CT=CT2;//-beta[index[1]]/(beta[index[0]]+eps)/1.5;
    // }
    // return result/sumGamma;

    ii = index[1];
    if (beta[ii] < CT * (beta[index[0]] + eps)) {
        sumGamma += gamma[ii];
        result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
    } else {
        // real extremPoint=0.5*(q[1]-q[3])/(q[1]-2*q[2]+q[3]+eps);
        //  int i0=index[0];
        //  if(!flags[i0]) return q[2];
        return result / sumGamma;
    }
    // {if(flag) return result/sumGamma;
    //  else
    //  return q[2];}

    ii = index[2];
    if (beta[ii] < (CT) * (beta[index[0]] + eps)) {
        sumGamma += gamma[ii];
        result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
    }
    result /= sumGamma;

    return result;
}

constexpr real Teno5_CongSortPositive(std::array<real, 5> q)
{
    real eps = 1e-12;
    std::array<real (*)(real, real, real), 3> u = { &u1, &u2, &u3 };

    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    // 排序
    std::array<real, 3> index = { 0, 1, 2 };
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        });
    // if(beta[index[2]]<eps) return q[2];
    int ii = index[0];
    real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];
    if (result < 0)
        return q[2];

    real extremPoint = 0.5 * (q[1] - q[3]) / (q[1] - 2 * q[2] + q[3] + eps);
    real extremPoint2 = 0.5 + (q[2] - q[3]) / (q[2] - 2 * q[3] + q[4] + eps);
    real factor = std::abs((q[2] - q[1] + eps) / (q[2] - q[3] + eps));
    bool flag = (extremPoint >= 1 + eps && extremPoint >= 0 - eps) || (extremPoint2 <= 1 + eps && extremPoint2 >= 0 - eps); //;|| (factor>6);
    real CT2 = factor, CT, sumGamma = gamma[index[0]];
    if (flag) {
        CT2 = 4;
        CT = CT2 - beta[index[1]] / (beta[index[0]] + eps);
    } else {
        CT2 = 20;
        CT = CT2 - beta[index[1]] / (beta[index[0]] + eps) / 1.5;
    }

    // return result/sumGamma;

    ii = index[1];
    real u1 = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
    if (beta[ii] < CT * (beta[index[0]]) && u1 > 0) {
        sumGamma += gamma[ii];
        result += gamma[ii] * u1;
    } else {
        if (!flag)
            return result / sumGamma;
        else
            return q[2];
    }

    ii = index[2];
    real u2 = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
    if (beta[ii] < (CT) * (beta[index[0]] + eps) && u1 > 0) {
        sumGamma += gamma[ii];
        result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
    }
    result /= sumGamma;

    return result;
}

constexpr real Teno5_CongSortabs(std::array<real, 5> q)
{
    real eps = 1e-12;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    // 排序
    std::array<real, 3> index = { 0, 1, 2 };
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        });

    std::array<real (*)(real, real, real), 3> u = { &u1, &u2, &u3 };

    real factor = std::abs((q[2] - q[1] + eps) / (q[2] - q[3] + eps));
    real CT2 = factor, CT, sumGamma = gamma[index[0]];
    if (factor > 6) {
        CT2 = 6;
        CT = CT2 - beta[index[1]] / (beta[index[0]] + eps);
    } else {
        CT2 = 20;
        CT = CT2 - beta[index[1]] / (beta[index[0]] + eps) / 1.5;
    }
    int ii = index[0];
    real result = (*u[ii])(q[ii], q[ii + 1], q[ii + 2]) * gamma[index[0]];
    // return result/sumGamma;

    if (beta[index[2]] < eps)
        return q[2];
    // real critical=pow(q[2]-q[1],2);
    // CT2=std::max(20.0*pow(critical,2.0),4.0);

    // CT=4;
    ii = index[1];
    if (beta[ii] < CT * (beta[index[0]])) {
        sumGamma += gamma[ii];
        result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
    } else {
        if (factor < 10)
            return result / sumGamma;
        else
            return q[2];
    }

    ii = index[2];
    if (beta[ii] < (CT) * (beta[index[0]] + eps)) {
        sumGamma += gamma[ii];
        result += gamma[ii] * (*u[ii])(q[ii], q[ii + 1], q[ii + 2]);
    }
    result /= sumGamma;

    return result;
}

constexpr real Teno5_CongIncrease(std::array<real, 5> q)
{
    real eps = 1e-20;
    std::array<real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    std::array<real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2)
        + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2)
        + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);

    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2)
        + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    // 排序
    std::array<real, 3> index = { 0, 1, 2 };
    std::sort(index.begin(), index.end(),
        [&](const int& a, const int& b) {
            return (beta[a] < beta[b]);
        });

    // 过于光滑
    if (beta[index[2]] <= 1e-10)
        return q[2];

    std::array<real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    real u0 = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    real sumbeta = beta[index[0]], sumGamma = gamma[index[0]];
    real result = u[index[0]] * gamma[index[0]];

    int ii = index[1];
    real CT = 10, a1 = (beta[index[0]] + eps);
    if (beta[ii] < CT * a1) {
        sumGamma += gamma[ii];
        result += gamma[ii] * u[ii];
        sumbeta += beta[ii];
    } else {
        return result / sumGamma;
    }

    real CT2 = beta[ii];
    ii = index[2];
    real ddd = 3;
    if (pow(beta[ii] * a1, ddd) - 5 * pow(CT2 * CT2, ddd) < 0) {
        sumGamma += gamma[ii];
        result += gamma[ii] * u[ii];
        sumbeta += beta[ii];
    } else {
        return result / sumGamma;
    }

    return result / sumGamma;
}

constexpr real musclInterpolation(real q1, real q2, real q3)
{

    real delta;
    real deltam, deltap;
    deltam = q2 - q1;
    deltap = q3 - q2;

    // minmod
    real beta = 1.0;
    if (deltap > 0) {
        delta = std::max(0.0, std::max(std::min(beta * deltam, deltap), std::min(deltam, beta * deltap)));
    } else {
        delta = std::min(0.0, std::min(std::max(beta * deltam, deltap), std::max(deltam, beta * deltap)));
    }
    return q2 + delta * 0.5;
}

constexpr std::array<real, 2> THINC(real q1, real q2, real q3)
{
    if ((q1 - q2) * (q2 - q3) < 0)
        return { q2, q2 };

    real qmax, qmin;
    if (q1 > q3) {
        qmax = q1;
        qmin = q3;
    } else {
        qmax = q3;
        qmin = q1;
    }
    real beta = 2;
    real T1 = exp(2.0 * beta), T3 = exp(-2.0 * beta);
    real f = (qmax - q2) / (q2 - qmin);
    T1 *= f;
    T3 *= f;
    real qBar = (q1 + q3) / 2, dq = (q1 - q3) / 2;
    return {
        qBar + dq * ((1 - T1) / (1 + T1)),
        qBar + dq * ((1 - T3) / (1 + T3))
    };
}

constexpr real THINC1(real q1, real q2, real q3)
{
    if ((q1 - q2) * (q2 - q3) < 1e-20)
        return q2;

    real qmax, qmin;
    if (q1 > q3) {
        qmax = q1;
        qmin = q3;
    } else {
        qmax = q3;
        qmin = q1;
    }
    real beta = 1.2;
    real T1 = exp(2.0 * beta), T3 = exp(-2.0 * beta);
    real f = (qmax - q2) / (q2 - qmin);
    T1 *= f;
    T3 *= f;
    real qBar = (q1 + q3) / 2, dq = (q1 - q3) / 2;
    return qBar + dq * ((1 - T3) / (1 + T3));
}

constexpr real Linear5(std::array<real, 5> q)
{
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
}


// 【王鸿飞】begin-1
constexpr real whf_TCNS_A(std::array<real, 5> q)
{
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0/1.0 * pow2(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2]) +
              1.0/4.0 * pow2(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2]);
    beta[1] = 1.0/1.0 * pow2(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3]) + 
              1.0/4.0 * pow2(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3]);
    beta[2] = 1.0/1.0 * pow2(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4]) + 
              1.0/4.0 * pow2(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4]);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);

    // 计算自适应阈值CT_A
    // real alpha1 = 10.0, alpha2 = 5.0;  //参数
    real alpha1 = 9.0, alpha2 = 4.0;  //参数
    // real xi = 1e-3, C_r = 0.24;
    // real epsilon_A = (0.9*C_r)/(1.0 - 0.9*C_r)*xi*xi;
    real epsilon_A = 2.7551*1e-7;

    std::array<real, 4> delta_q;
    delta_q[0] = q[0] - q[1];
    delta_q[1] = q[1] - q[2];
    delta_q[2] = q[2] - q[3];
    delta_q[3] = q[3] - q[4];

    std::array<real, 3> eta;
        eta[0] = (std::abs(2.0*delta_q[0]*delta_q[1]) + epsilon_A)
              / (delta_q[0]*delta_q[0] + delta_q[1]*delta_q[1] + epsilon_A);
        eta[1] = (std::abs(2.0*delta_q[1]*delta_q[2]) + epsilon_A)
              / (delta_q[1]*delta_q[1] + delta_q[2]*delta_q[2] + epsilon_A);
        eta[2] = (std::abs(2.0*delta_q[2]*delta_q[3]) + epsilon_A)
              / (delta_q[2]*delta_q[2] + delta_q[3]*delta_q[3] + epsilon_A);
    
    // 计算η_min值
    real eta_min = std::min({eta[0],eta[1],eta[2]});

    // real min = std::min(1.0, eta_min/0.24);
    real min = std::min(1.0, 4.166667*eta_min);
    real gm = min*min*min*min*(5.0 - 4*min);
    int beta_A = std::floor(alpha1 - alpha2*(1.0 - gm));

    // real m = 1.0 - std::min(1.0, eta_min/0.24);
    // real gm = (1.0-m)*(1.0-m)*(1.0-m)*(1.0-m)*(1.0 + 4*m);
    // int beta_A = std::floor(alpha1 - alpha2*(1.0 - gm));
    
    real CT_A;

    
    // 求光滑度量gamma
    const real C = 1.0;
    const real eps = 1e-40;

    std::array<real, 3> gamma;

    gamma[0] = C + tau/(beta[0] + eps);
    gamma[1] = C + tau/(beta[1] + eps);
    gamma[2] = C + tau/(beta[2] + eps);

    gamma[0] = gamma[0]*gamma[0];
    gamma[1] = gamma[1]*gamma[1];
    gamma[2] = gamma[2]*gamma[2];

    gamma[0] = gamma[0]*gamma[0]*gamma[0];
    gamma[1] = gamma[1]*gamma[1]*gamma[1];
    gamma[2] = gamma[2]*gamma[2]*gamma[2];

    real gamma_sum = gamma[0] + gamma[1] + gamma[2];
    
    switch(beta_A) {
    case 5: CT_A = 1e-5; 
    break;
    case 6: CT_A = 1e-6;
    break;
    case 7: CT_A = 1e-7;
    break;
    case 8: CT_A = 1e-8;
    break;
    case 9: CT_A = 1e-9;
    break;
    case 10: CT_A = 1e-10;
    break;
    default: pandaun_001 = true ;
    break;
    }

    // if (CTA_counter_on_off)
    // {
    //   switch(beta_A) {
    //   case 5: global_counter_5++; break;
    //   case 6: global_counter_6++; break;
    //   case 7: global_counter_7++; break;
    //   case 8: global_counter_8++; break;
    //   case 9: global_counter_9++; break;
    //   case 10: global_counter_10++; break;
    //   }
    // }
    
    real rr= CT_A*gamma_sum;
    // 构造截止函数
    unsigned short flag = 0;
    if (gamma[0] < rr) flag += 1;
    if (gamma[1] < rr) flag += 2;
    if (gamma[2] < rr) flag += 4;

    // 插值
    switch (flag) {
    case 0:  // 111
        return 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] + 15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    case 1:  // 011
        return -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] - 1.0/16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    case 4:  // 110
        return 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] + 5.0/16.0 * q[3];
    case 5:  // 010
        return -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    case 6:  // 100
        return 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}





















// 【王鸿飞】F2.0时代













constexpr real whf_TCNS_AS_Ff3_5_9_time_improve(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow2(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2]) +
          1.0 / 4.0 * pow2(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2]),

      1.0 / 1.0 * pow2(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3]) +
          1.0 / 4.0 * pow2(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3]),

      1.0 / 1.0 * pow2(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4]) +
          1.0 / 4.0 * pow2(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4])};

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);
  unsigned short minBeta =
      std::min_element(beta.begin(), beta.end()) - beta.begin();
  //CT adapt begin
  // real xi = 1e-3;
  // real Cr = 0.24;
  // 此时 epsilon_A = 2.7551*1e-7 ， 1/Cr = 4.1666667

  real epsilon_A = 2.7551*1e-7;

    std::array<real, 4> delta_q;
    delta_q[0] = q[0] - q[1];
    delta_q[1] = q[1] - q[2];
    delta_q[2] = q[2] - q[3];
    delta_q[3] = q[3] - q[4];
    
  // 计算η值
  real eta_im1 = (std::abs(2.0*delta_q[1]*delta_q[0]) + epsilon_A) / 
                  (pow2(delta_q[1]) + pow2(delta_q[0]) + epsilon_A);
  
  real eta_i = (std::abs(2.0*delta_q[2]*delta_q[1]) + epsilon_A) / 
                (pow2(delta_q[2]) + pow2(delta_q[1]) + epsilon_A);
  
  real eta_ip1 = (std::abs(2.0*delta_q[3]*delta_q[2]) + epsilon_A) / 
                  (pow2(delta_q[3]) + pow2(delta_q[2]) + epsilon_A);
  
  real eta_min = std::min({eta_im1, eta_i, eta_ip1});

  // 计算min
  // real min = std::min(1.0, 4.166667*eta_min);
  // real m = 1 - min;
  
  real x = std::min(0.24, eta_min);

  // 计算C_T_prime
  real xx=x*x;
  real xxx=x*xx;
  // // 三次拟合
  // // Ff3_5_10

  // // Ff3_5_9
  real C_T = 26.778*xxx - 11.3915*xx + 0.700688*x + 0.148606;
  // // Ff3_4_10
  // // Ff3_4_9
  // //CT adapt end
  
  real CT_1 = 1 - C_T;
  real tau = std::abs(beta[2] -
                      beta[0]); //,KK=0.15704178024750198*(beta[minBeta]+tau);
  real rr = C_T * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}




constexpr real whf_TCNS_LADS_Ff2_5_10(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  real tau = std::abs(beta[2] - beta[0]); 

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);

  unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();

  real HminBeta = 10*minBeta;
  
  // real theta = HminBeta/(HminBeta + tau + 1e-40);
  real theta = HminBeta/(HminBeta + tau );

  // 计算C_T_prime
  // ff2-5-10
  real C_T = 0.117931*theta*theta - 0.243913*theta + 0.152466;
  // ff2-5-9
  // real C_T = 0.089408*theta*theta - 0.207764*theta + 0.154337;
  
  //CT adapt end
  
  real CT_1 = 1 - C_T;
  real rr = C_T * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}


constexpr real whf_TCNS_LAD(std::array<real, 5> q) {
    // 局部光滑因子β
    std::array<real, 3> beta;
    beta[0] = 1.0/1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
              1.0/4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0/1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 
              1.0/4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0/1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 
              1.0/4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);
    // 全局光滑因子τ
    real tau = std::abs(beta[2] - beta[0]);

    // 计算自适应阈值CT_A
    const real eps = 1e-40;
    real b_u = 10.0, b_l = 5.0;  //参数
    

    // 计算 chi_max 值(方案1)
    // std::array<real, 3> chi;
    // chi[0] = tau/(beta[0] + eps);
    // chi[1] = tau/(beta[1] + eps);
    // chi[2] = tau/(beta[2] + eps);
    // real chi_max = std::max({chi[0],chi[1],chi[2]});
    // real theta = 1/(1 + chi_max/10);
    // int m = b_l + std::floor((b_u - b_l)*theta);
    
    // 计算 chi_max 值(方案1的优化)
    real minbeta = std::min({beta[0],beta[1],beta[2]});
    real temp = 10*(minbeta + eps);
    real theta = temp/(temp + tau);
    int m = b_l + std::floor((b_u - b_l)*theta);


    real CT_A;
    
    // 求光滑度量gamma
    const real C = 1.0;

    std::array<real, 3> gamma;

    gamma[0] = C + tau/(beta[0] + eps);
    gamma[1] = C + tau/(beta[1] + eps);
    gamma[2] = C + tau/(beta[2] + eps);

    gamma[0] = gamma[0]*gamma[0];
    gamma[1] = gamma[1]*gamma[1];
    gamma[2] = gamma[2]*gamma[2];

    gamma[0] = gamma[0]*gamma[0]*gamma[0];
    gamma[1] = gamma[1]*gamma[1]*gamma[1];
    gamma[2] = gamma[2]*gamma[2]*gamma[2];

    real gamma_sum = gamma[0] + gamma[1] + gamma[2];
    
    switch(m) {
    case 5: CT_A = 1e-5; 
    break;
    case 6: CT_A = 1e-6;
    break;
    case 7: CT_A = 1e-7;
    break;
    case 8: CT_A = 1e-8;
    break;
    case 9: CT_A = 1e-9;
    break;
    case 10: CT_A = 1e-10;
    break;
    default: pandaun_001 = true ;
    break;
    }
    
    real rr= CT_A*gamma_sum;
    // 构造截止函数
    unsigned short flag = 0;
    if (gamma[0] < rr) flag += 1;
    if (gamma[1] < rr) flag += 2;
    if (gamma[2] < rr) flag += 4;

    // 插值
    switch (flag) {
    case 0:  // 111
        return 3.0/128.0 * q[0] - 5.0/32.0 * q[1] + 45.0/64.0 * q[2] + 15.0/32.0 * q[3] - 5.0/128.0 * q[4];
    case 1:  // 011
        return -1.0/16.0 * q[1] + 9.0/16.0 * q[2] + 9.0/16.0 * q[3] - 1.0/16.0 * q[4];
    case 2:  // 101
    case 3:  // 001
        return 3.0/8.0 * q[2] + 3.0/4.0 * q[3] - 1.0/8.0 * q[4];
    case 4:  // 110
        return 1.0/16.0 * q[0] - 5.0/16.0 * q[1] + 15.0/16.0 * q[2] + 5.0/16.0 * q[3];
    case 5:  // 010
        return -1.0/8.0 * q[1] + 3.0/4.0 * q[2] + 3.0/8.0 * q[3];
    case 6:  // 100
        return 3.0/8.0 * q[0] - 5.0/4.0 * q[1] + 15.0/8.0 * q[2];
    default: // 000: 强间断区(直接取值)
        return q[2];
    }
}


constexpr real whf_TCNS_LADS_g1(std::array<real, 5> q) {
  real eps = 1e-40; // 1e-10;
  std::array<real, 3> beta = {
      1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) +
          1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),

      1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) +
          1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),

      1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) +
          1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)};

  real tau = std::abs(beta[2] - beta[0]); 

  // int minBeta=(beta[0]>beta[1])? ((beta[2]>beta[1])? 1:
  // 2):((beta[2]>beta[0])? 0 : 2);

  unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();

  real y = 1.0 + tau/(10*minBeta + eps);


  // 计算C_T_prime
  real C_T = (0.157042*y - 0.096862)/(y + 1.6115);
  
  //CT adapt end
  
  real CT_1 = 1 - C_T;
  real rr = C_T * tau - CT_1 * beta[minBeta];
  real ll = tau * beta[minBeta];
  // unsigned
  // flag=(minBeta!=0&&ll<rr*beta[0])+((minBeta!=1&&ll<rr*beta[1])<<1)+((minBeta!=2&&ll<rr*beta[2])<<2);
  unsigned short flag = 0;
  if (ll < rr * beta[0])
    flag += 1;
  if (ll < rr * beta[1])
    flag += 2;
  if (ll < rr * beta[2])
    flag += 4;
  switch (flag) {
  case 0:
    /* 1,1,1 */
    return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] +
           15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    break;
  case 1:
    /* 0,1,1 */
    return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] -
           1.0 / 16.0 * q[4];
    break;
  case 2:
    /* 1,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 3:
    /* 0,0,1 */
    return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    break;
  case 4:
    /* 1,1,0 */
    return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] +
           5.0 / 16.0 * q[3];
    break;
  case 5:
    /* 0,1,0 */
    return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    break;
  case 6:
    /* 1,0,0 */
    return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    break;
  default:
    /* 0,0,0 */
    return q[2];
    break;
  }
}












































// 【王鸿飞】end
