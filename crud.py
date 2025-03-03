# 对数据库CRUD操作的包装函数族
# 张景栋，创建于2022.2.14，最后修改于2022.2.17
import sqlite3

import numpy as np
import pandas as pd


def add_vlbi_source(conn, paras):
    """向vlbi_source表添加新的射电星

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        paras (array-like): 参数
    """
    # 可能误差为nan的情况
    ra_e = 'NULL'
    if not np.isnan(paras[3]):
        ra_e = f"\'{paras[3]}\'"
    dec_e = 'NULL'
    if not np.isnan(paras[5]):
        dec_e = f"\'{paras[5]}\'"
    # 如果输入paras不包含相关系数，则全部置零
    if len(paras) < 13:
        sql = f"INSERT INTO vlbi_source(SOURCE_NAME, EPOCH, RA, RA_ERROR, DEC, DEC_ERROR, PI, PI_ERROR, PMRA, " \
              f"PMRA_ERROR, PMDEC, PMDEC_ERROR, RA_DEC_CORR, RA_PI_CORR, RA_PMRA_CORR, RA_PMDEC_CORR, DEC_PI_CORR, " \
              f"DEC_PMRA_CORR, DEC_PMDEC_CORR, PI_PMRA_CORR, PI_PMDEC_CORR, PMRA_PMDEC_CORR) " \
              f"VALUES (\'{paras[0]}\', {paras[1]}, {paras[2]}, {ra_e}, {paras[4]}, {dec_e}, {paras[6]}, " \
              f"{paras[7]}, {paras[8]}, {paras[9]}, {paras[10]}, {paras[11]}, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)"
    else:
        sql = f"INSERT INTO vlbi_source(SOURCE_NAME, EPOCH, RA, RA_ERROR, DEC, DEC_ERROR, PI, PI_ERROR, PMRA, " \
              f"PMRA_ERROR, PMDEC, PMDEC_ERROR, RA_DEC_CORR, RA_PI_CORR, RA_PMRA_CORR, RA_PMDEC_CORR, DEC_PI_CORR, " \
              f"DEC_PMRA_CORR, DEC_PMDEC_CORR, PI_PMRA_CORR, PI_PMDEC_CORR, PMRA_PMDEC_CORR) " \
              f"VALUES (\'{paras[0]}\', {paras[1]}, {paras[2]}, {ra_e}, {paras[4]}, {dec_e}, {paras[6]}, " \
              f"{paras[7]}, {paras[8]}, {paras[9]}, {paras[10]}, {paras[11]}, {paras[12]}, {paras[13]}, {paras[14]}, " \
              f"{paras[15]}, {paras[16]}, {paras[17]}, {paras[18]}, {paras[19]}, {paras[20]}, {paras[21]})"
    conn.execute(sql)
    conn.commit()


def delete_vlbi_source(conn, name=None):
    """删除VLBI source

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        name (str, optional): 指定的星的名字. Defaults to None.
    """
    if name is None:
        sql = f"DELETE FROM vlbi_source"
    else:
        sql = f"DELETE FROM vlbi_source WHERE SOURCE_NAME=\'{name}\'"
    conn.execute(sql)
    conn.commit()
        

def add_vlbi_obs(conn, paras):
    """向vlbi_observation表添加新的观测

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        paras (array-like): 参数
    """
    component = 'NULL'
    if isinstance(paras[2], str):
        component = f"\'{paras[2]}\'"
    # component = f"\'{paras[2]}\'"
    sql = f"INSERT INTO vlbi_observation(SOURCE_NAME, EPOCH, COMPONENT, CALIBRATOR, RA, RA_ERROR, DEC, DEC_ERROR, RA_DEC_CORR) " \
          f"VALUES (\'{paras[0]}\', {paras[1]}, {component}, \'{paras[3]}\', {paras[4]}, {paras[5]}, {paras[6]}, {paras[7]}, {paras[8]})"
    conn.execute(sql)
    conn.commit()


def get_vlbi_obs_source_list(conn, min_num=None):
    """获取全部至少有min_num次观测数据的射电星名字列表

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        min_num (int): 最低观测历元数. Defaults to None.

    Returns:
        List: 射电星名字
    """
    if min_num is None:
        sql = f"SELECT DISTINCT SOURCE_NAME FROM vlbi_observation WHERE RA is not null"
    else:
        sql = f"SELECT SOURCE_NAME FROM (SELECT SOURCE_NAME, COUNT(*) c FROM vlbi_observation GROUP BY SOURCE_NAME) WHERE c>={min_num:d};"
    cur = conn.execute(sql)
    s_name_list = []
    for row in cur:
        s_name_list.append(row[0])
    return s_name_list


def get_vlbi_source(conn, name=None):
    """获取VLBI五参数观测数据

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        name (str): 指定的星的名字. Defaults to None.

    Returns:
        Dataframe: 五参数数据表
    """
    if name is None:
        sql = f"SELECT * FROM vlbi_source WHERE RA is not null"
    else:
        sql = f"SELECT * FROM vlbi_source WHERE (RA is not null) and (SOURCE_NAME = \'{name}\')"
    cur = conn.execute(sql)
    df = pd.DataFrame(columns=['SOURCE_NAME', 'EPOCH', 'RA', 'RA_ERROR', 'DEC', 'DEC_ERROR', 'PI', 'PI_ERROR',
                               'PMRA', 'PMRA_ERROR', 'PMDEC', 'PMDEC_ERROR', 'RA_DEC_CORR', 'RA_PI_CORR',
                               'RA_PMRA_CORR', 'RA_PMDEC_CORR', 'DEC_PI_CORR', 'DEC_PMRA_CORR', 'DEC_PMDEC_CORR',
                               'PI_PMRA_CORR', 'PI_PMDEC_CORR', 'PMRA_PMDEC_CORR'])
    for row in cur:
        df.loc[df.index.size] = row
    return df


def get_vlbi_obs(conn, name=None):
    """获取VLBI单历元观测数据

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        name (str): 指定的星的名字. Defaults to None.

    Returns:
        Dataframe: 单历元数据表
    """
    if name is None:
        sql = f"SELECT * FROM vlbi_observation WHERE RA is not null"
    else:
        sql = f"SELECT * FROM vlbi_observation WHERE (RA is not null) and (SOURCE_NAME = \'{name}\')"
    cur = conn.execute(sql)
    df = pd.DataFrame(columns=['SOURCE_NAME', 'EPOCH', 'COMPONENT', 'CALIBRATOR',
                               'RA', 'RA_ERROR', 'DEC', 'DEC_ERROR', 'RA_DEC_CORR'])
    for row in cur:
        df.loc[df.index.size] = row
    return df


def add_gaia(conn, g_data):
    """向数据库添加Gaia数据行

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        g_data (Series): Gaia数据行
    """
    sql = f"INSERT INTO gaia(SOURCE_NAME, DESIGNATION, EPOCH, GMAG, RA, RA_ERROR, DEC, DEC_ERROR, PI, PI_ERROR, PMRA, " \
          f"PMRA_ERROR, PMDEC, PMDEC_ERROR, RA_DEC_CORR, RA_PI_CORR, RA_PMRA_CORR, RA_PMDEC_CORR, DEC_PI_CORR, " \
          f"DEC_PMRA_CORR, DEC_PMDEC_CORR, PI_PMRA_CORR, PI_PMDEC_CORR, PMRA_PMDEC_CORR) " \
          f"VALUES (\'{g_data['SOURCE_NAME']}\', \'{g_data['DESIGNATION']}\', \'{g_data['EPOCH']}\', \'{g_data['GMAG']}\', {g_data['RA']}, " \
          f"{g_data['RA_ERROR']}, {g_data['DEC']}, {g_data['DEC_ERROR']}, {g_data['PI']}, {g_data['PI_ERROR']}, " \
          f"{g_data['PMRA']}, {g_data['PMRA_ERROR']}, {g_data['PMDEC']}, {g_data['PMDEC_ERROR']}, " \
          f"{g_data['RA_DEC_CORR']}, {g_data['RA_PI_CORR']}, {g_data['RA_PMRA_CORR']}, {g_data['RA_PMDEC_CORR']}, " \
          f"{g_data['DEC_PI_CORR']}, {g_data['DEC_PMRA_CORR']}, {g_data['DEC_PMDEC_CORR']}, " \
          f"{g_data['PI_PMRA_CORR']}, {g_data['PI_PMDEC_CORR']}, {g_data['PMRA_PMDEC_CORR']})"
    conn.execute(sql)
    conn.commit()


def get_gaia(conn, name=None, gmag=None):
    """获取Gaia数据

    Args:
        conn (sqlite3.Connection): 数据库连接实例
        name (str, optional): 指定的星的名字. Defaults to None.
        gmag (tuple, optional): 星等范围. Defaults to None.

    Returns:
        Dataframe: Gaia数据表
    """
    if (name is None) and (gmag is None):
        sql = f"SELECT * FROM gaia"
    elif (name is not None) and (gmag is None):
        sql = f"SELECT * FROM gaia WHERE SOURCE_NAME = \'{name}\'"
    elif (name is None) and (gmag is not None):
        if gmag[0] is None:
            sql = f"SELECT * FROM gaia WHERE GMAG <= \'{gmag[1]}\'"
        elif gmag[1] is None:
            sql = f"SELECT * FROM gaia WHERE GMAG >= \'{gmag[0]}\'"
        else:
            sql = f"SELECT * FROM gaia WHERE (GMAG >= \'{gmag[0]}\') AND (GMAG <= \'{gmag[1]}\')"
    else:
        if gmag[0] is None:
            sql = f"SELECT * FROM gaia WHERE (SOURCE_NAME = \'{name}\') AND (GMAG <= \'{gmag[1]}\')"
        elif gmag[1] is None:
            sql = f"SELECT * FROM gaia WHERE (SOURCE_NAME = \'{name}\') AND (GMAG >= \'{gmag[0]}\')"
        else:
            sql = f"SELECT * FROM gaia WHERE (SOURCE_NAME = \'{name}\') AND (GMAG >= \'{gmag[0]}\') AND (GMAG <= \'{gmag[1]}\')"
    cur = conn.execute(sql)
    df = pd.DataFrame(columns=['SOURCE_NAME', 'DESIGNATION', 'EPOCH', 'GMAG', 'RA', 'RA_ERROR', 'DEC', 'DEC_ERROR', 'PI',
                               'PI_ERROR', 'PMRA', 'PMRA_ERROR', 'PMDEC', 'PMDEC_ERROR', 'RA_DEC_CORR', 'RA_PI_CORR',
                               'RA_PMRA_CORR', 'RA_PMDEC_CORR', 'DEC_PI_CORR', 'DEC_PMRA_CORR', 'DEC_PMDEC_CORR',
                               'PI_PMRA_CORR', 'PI_PMDEC_CORR', 'PMRA_PMDEC_CORR'])
    for row in cur:
        df.loc[df.index.size] = row
    return df


def dbClear(conn):
    """清空数据库

    Args:
        conn (sqlite3.Connection): 数据库连接实例
    """ 
    sql = "DELETE FROM gaia"
    conn.execute(sql)
    sql = "DELETE FROM vlbi_observation"
    conn.execute(sql)
    sql = "DELETE FROM vlbi_source"
    conn.execute(sql)
    conn.commit()


if __name__ == "__main__":
    connection = sqlite3.connect('db/com.db')
    # add_source(connection, 'test2')
    # add_obs(connection, 'test2', '2022-02-14 14:09')
    # add_obs(connection, 'test2', 2459602.3)
    # obs_para_update(connection, 'test1', 2459599.20763889)
    # rua = gregory2julian('2022-02-14 14:09')
    # s_list = get_vlbi_source_list(connection)
    # rua = get_vlbi_obs(connection)
    # r_name, ra, ra_e, dec, dec_e, d_ra, d_ra_e, d_dec, d_dec_e
    # ll = [['ref_test1', 1.1, 0.2, -0.9, 0.3], ['ref_test2', -1.05, 0.8, 0.98, 0.6]]
    # add_ref_obs(connection, 'test1', '2022-02-15 11:09', ll)
    # g = get_gaia(connection)
    # v = get_vlbi_obs(connection)
    v = get_vlbi_source(connection)
    connection.close()
    pass
