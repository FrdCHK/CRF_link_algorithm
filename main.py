"""
description
@Author: Jingdong Zhang
@DATE  : 2025/3/3
"""
import numpy as np
import pandas as pd
import sqlite3

from getData import get_by_name_list
from result_print import result_print
from connection import solveAll


if __name__ == '__main__':
    fname = "sample"
    connect = sqlite3.connect(f'data/{fname}.db')
    sel = pd.read_csv(f"list/solution_0G10.5.txt", comment='#')

    dt = get_by_name_list(connect, sel)
    connect.close()
    res_x, res_cov, res_q, res_niu, res_e, res_o = solveAll(dt)
    q = np.array([item[0] for item in res_q])
    result_print("res: value±std(normalized)", res_x, res_cov, q, res_niu)

    print("---------------------------------------------")

    # 计算相关系数矩阵
    unit_mul = (3.6e6*180/np.pi)**2
    res_cov_mas = res_cov * unit_mul
    inv_d = np.linalg.inv(np.diag(np.diag(res_cov_mas) ** 0.5))
    r = inv_d @ res_cov_mas @ inv_d
    blank = "   ... "
    print(f"cov matrix:\n{r[0, 0]:+7.3f}{r[0, 1]:+7.3f}{r[0, 2]:+7.3f}{r[0, 3]:+7.3f}{r[0, 4]:+7.3f}{r[0, 5]:+7.3f}\n"
          f"{blank}{r[1, 1]:+7.3f}{r[1, 2]:+7.3f}{r[1, 3]:+7.3f}{r[1, 4]:+7.3f}{r[1, 5]:+7.3f}\n"
          f"{blank}{blank}{r[2, 2]:+7.3f}{r[2, 3]:+7.3f}{r[2, 4]:+7.3f}{r[2, 5]:+7.3f}\n"
          f"{blank}{blank}{blank}{r[3, 3]:+7.3f}{r[3, 4]:+7.3f}{r[3, 5]:+7.3f}\n"
          f"{blank}{blank}{blank}{blank}{r[4, 4]:+7.3f}{r[4, 5]:+7.3f}\n"
          f"{blank}{blank}{blank}{blank}{blank}{r[5, 5]:+7.3f}")

    print("---------------------------------------------")

    # 计算orientation和spin权中位数
    res_e_mas2 = res_e/unit_mul
    res_o_mas2 = res_o/unit_mul
    print(f"median E: {np.median(res_e_mas2):6.1f} median O: {np.median(res_o_mas2):6.1f}")

    print("---------------------------------------------")

    # normalized Q for each source
    Q_n = q/res_niu
    print("SOURCE_NAME         nu    Q/nu       E   Omega")
    for j in range(Q_n.size):
        print(f"{dt[j]['SOURCE_NAME']:<18s} {res_niu[j]:3d} {Q_n[j]:7.2f} {res_e_mas2[j]:7.1f} {res_o_mas2[j]:7.1f}")
