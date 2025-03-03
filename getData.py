"""
从数据库获取VLBI多历元数据, 并生成对应的Gaia外推数据
@Author: Jingdong Zhang
@DATE  : 2022/2/21
"""
import sqlite3
import numpy as np
import pandas as pd

from units import mas2rad
from crud import get_gaia, get_vlbi_obs, get_vlbi_source


def get_by_name_list(conn, names, gmag=None):
    data_list = []
    for i, item in names.iterrows():
        g = get_gaia(conn, item[0], gmag)
        v5 = get_vlbi_source(conn, item[0])
        v2 = get_vlbi_obs(conn, item[0])

        # 如果name list里的源找不到，避免append空item
        if g.empty or (v5.empty and v2.empty):
            continue

        # 转为弧度制
        g[['RA', 'DEC']] = np.deg2rad(g[['RA', 'DEC']])
        g[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR', 'PMDEC',
           'PMDEC_ERROR']] = mas2rad(g[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR',
                                        'PMDEC', 'PMDEC_ERROR']])
        v5[['RA', 'DEC']] = np.deg2rad(v5[['RA', 'DEC']])
        if not v5.empty:
            if (v5.at[0, 'RA_ERROR'] is None) | (v5.at[0, 'DEC_ERROR'] is None):
                v5[['PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR', 'PMDEC',
                    'PMDEC_ERROR']] = mas2rad(v5[['PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR',
                                                  'PMDEC', 'PMDEC_ERROR']])
            else:
                v5[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR', 'PMDEC',
                    'PMDEC_ERROR']] = mas2rad(v5[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR',
                                                  'PMDEC', 'PMDEC_ERROR']])
        v2[['RA', 'DEC']] = np.deg2rad(v2[['RA', 'DEC']])
        v2[['RA_ERROR', 'DEC_ERROR']] = mas2rad(v2[['RA_ERROR', 'DEC_ERROR']])
        data_list.append({"SOURCE_NAME": item[0], "REF_EPOCH": g.at[0, "EPOCH"], "GAIA": g, "VLBI_5P": v5, "VLBI_2P": v2})
    return data_list


def get_simu(conn, data_type):
    """
    get simulated data from db
    :param conn: db connection
    :param data_type: 'five'/'double'/'single'
    :return: simulated data
    """
    data_list = []
    g = get_gaia(conn)  # 必然有全的Gaia表，因此先读gaia表，并按本表中SOURCE_NAME顺序获取VLBI数据
    g[['RA', 'DEC']] = np.deg2rad(g[['RA', 'DEC']])
    g[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR',
       'PMRA', 'PMRA_ERROR', 'PMDEC', 'PMDEC_ERROR']] = mas2rad(g[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR',
                                                                   'PMRA', 'PMRA_ERROR', 'PMDEC', 'PMDEC_ERROR']])
    if data_type == "five":
        for i, item in g["SOURCE_NAME"].items():
            v5 = get_vlbi_source(conn, item)
            v5[['RA', 'DEC']] = np.deg2rad(v5[['RA', 'DEC']])
            v5[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR', 'PMDEC',
                'PMDEC_ERROR']] = mas2rad(v5[['RA_ERROR', 'DEC_ERROR', 'PI', 'PI_ERROR', 'PMRA', 'PMRA_ERROR',
                                              'PMDEC', 'PMDEC_ERROR']])
            v2 = pd.DataFrame(columns=['SOURCE_NAME', 'EPOCH', 'COMPONENT', 'CALIBRATOR',
                                       'RA', 'RA_ERROR', 'DEC', 'DEC_ERROR', 'RA_DEC_CORR'])
            g_row = g.iloc[i:i+1]
            g_row.reset_index(inplace=True, drop=True)
            data_list.append({"SOURCE_NAME": item, "REF_EPOCH": g.at[i, "EPOCH"],
                              "GAIA": g_row, "VLBI_5P": v5, "VLBI_2P": v2})
    elif data_type == "double":
        for i, item in g["SOURCE_NAME"].items():
            v2 = get_vlbi_obs(conn, item)
            v2 = v2.iloc[[0, -1]]
            v2.reset_index(inplace=True, drop=True)
            v2[['RA', 'DEC']] = np.deg2rad(v2[['RA', 'DEC']])
            v2[['RA_ERROR', 'DEC_ERROR']] = mas2rad(v2[['RA_ERROR', 'DEC_ERROR']])
            v5 = pd.DataFrame(columns=['SOURCE_NAME', 'EPOCH', 'RA', 'RA_ERROR', 'DEC', 'DEC_ERROR', 'PI', 'PI_ERROR',
                                       'PMRA', 'PMRA_ERROR', 'PMDEC', 'PMDEC_ERROR', 'RA_DEC_CORR', 'RA_PI_CORR',
                                       'RA_PMRA_CORR', 'RA_PMDEC_CORR', 'DEC_PI_CORR', 'DEC_PMRA_CORR',
                                       'DEC_PMDEC_CORR', 'PI_PMRA_CORR', 'PI_PMDEC_CORR', 'PMRA_PMDEC_CORR'])
            g_row = g.iloc[i:i+1]
            g_row.reset_index(inplace=True, drop=True)
            data_list.append({"SOURCE_NAME": item, "REF_EPOCH": g.at[i, "EPOCH"],
                              "GAIA": g_row, "VLBI_5P": v5, "VLBI_2P": v2})
    elif data_type == "single":
        for i, item in g["SOURCE_NAME"].items():
            v2 = get_vlbi_obs(conn, item)
            rand_index = np.random.choice([v2.index[0], v2.index[-1]])
            v2 = v2.iloc[rand_index:rand_index+1]
            v2.reset_index(inplace=True, drop=True)
            v2[['RA', 'DEC']] = np.deg2rad(v2[['RA', 'DEC']])
            v2[['RA_ERROR', 'DEC_ERROR']] = mas2rad(v2[['RA_ERROR', 'DEC_ERROR']])
            v5 = pd.DataFrame(columns=['SOURCE_NAME', 'EPOCH', 'RA', 'RA_ERROR', 'DEC', 'DEC_ERROR', 'PI', 'PI_ERROR',
                                       'PMRA', 'PMRA_ERROR', 'PMDEC', 'PMDEC_ERROR', 'RA_DEC_CORR', 'RA_PI_CORR',
                                       'RA_PMRA_CORR', 'RA_PMDEC_CORR', 'DEC_PI_CORR', 'DEC_PMRA_CORR',
                                       'DEC_PMDEC_CORR', 'PI_PMRA_CORR', 'PI_PMDEC_CORR', 'PMRA_PMDEC_CORR'])
            g_row = g.iloc[i:i+1]
            g_row.reset_index(inplace=True, drop=True)
            data_list.append({"SOURCE_NAME": item, "REF_EPOCH": g.at[i, "EPOCH"],
                              "GAIA": g_row, "VLBI_5P": v5, "VLBI_2P": v2})
    return data_list


if __name__ == '__main__':
    # connect = sqlite3.connect('data/Lind2020+Lunz2023.db')
    # sel = pd.read_csv("input/selected_sources_Lind2020+Lunz2023.txt", comment='#')
    # dt = get_by_name_list(connect, sel)

    connect = sqlite3.connect('data/simu_1x_50_2024_plx.db')
    dt = get_simu(connect, data_type="double")

    connect.close()

    # from connection import AllIn1
    # test = AllIn1(dt[21])
    # test = AllIn1(dt[28])
