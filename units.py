"""
description
@Author: Jingdong Zhang
@DATE  : 2022/8/26
"""
import numpy as np


def mas2rad(mas_data):
    """unit convertion: mas to rad

    Args:
        mas_data (float or array-like): data in mas

    Returns:
        float or ndarray: data in rad
    """
    return np.deg2rad(mas_data/3.6e6)


def mas2deg(mas_data):
    """unit convertion: mas to deg

    Args:
        deg_data (float or array-like): data in mas

    Returns:
        float or ndarray: data in deg
    """
    return mas_data/3.6e6


def rad2mas(rad_data):
    """unit convertion: rad to mas

    Args:
        rad_data (float or array-like): data in rad

    Returns:
        float or ndarray: data in mas
    """
    return np.rad2deg(rad_data)*3.6e6


def deg2mas(deg_data):
    """unit convertion: deg to mas

    Args:
        deg_data (float or array-like): data in deg

    Returns:
        float or ndarray: data in mas
    """
    return deg_data * 3.6e6


if __name__ == '__main__':
    import pandas as pd
    rua = pd.DataFrame([[None]])
    rua.fillna(np.nan, inplace=True)
    a = mas2rad(rua)
    pass
