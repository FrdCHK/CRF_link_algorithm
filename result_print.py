"""
description
@Author: Jingdong Zhang
@DATE  : 2023/12/14
"""
import numpy as np

from units import rad2mas


def result_print(pref, res_x, res_cov, res_q, res_niu):
    res_cov_mas = res_cov * ((180 / np.pi * 3.6e6) ** 2)
    std = np.diag(res_cov_mas) ** 0.5
    chi2 = np.sum(res_q) / (np.sum(res_niu))
    std_chi2 = std*np.sqrt(chi2)
    x_mas = rad2mas(res_x)
    print(f"{pref}\n"
          f"{x_mas[0, 0]:+6.3f}±{std[0]:.3f}({std_chi2[0]:.3f}) {x_mas[1, 0]:+6.3f}±{std[1]:.3f}({std_chi2[1]:.3f}) "
          f"{x_mas[2, 0]:+6.3f}±{std[2]:.3f}({std_chi2[2]:.3f}) {x_mas[3, 0]:+6.3f}±{std[3]:.3f}({std_chi2[3]:.3f}) "
          f"{x_mas[4, 0]:+6.3f}±{std[4]:.3f}({std_chi2[4]:.3f}) {x_mas[5, 0]:+6.3f}±{std[5]:.3f}({std_chi2[5]:.3f})\n"
          f"chi2:{chi2:7.3f}")


def ori_print(pref, res_x, res_cov, res_q, res_niu):
    res_cov_mas = res_cov * ((180 / np.pi * 3.6e6) ** 2)
    std = np.diag(res_cov_mas) ** 0.5
    chi2 = np.sum(res_q) / (np.sum(res_niu))
    std_chi2 = std*np.sqrt(chi2)
    x_mas = rad2mas(res_x)
    print(f"{pref}\n"
          f"{x_mas[0, 0]:+6.3f}±{std[0]:.3f}({std_chi2[0]:.3f}) {x_mas[1, 0]:+6.3f}±{std[1]:.3f}({std_chi2[1]:.3f}) "
          f"{x_mas[2, 0]:+6.3f}±{std[2]:.3f}({std_chi2[2]:.3f})\n"
          f"chi2:{chi2:7.3f}")


if __name__ == "__main__":
    pass
