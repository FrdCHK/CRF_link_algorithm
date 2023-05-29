"""
由相关系数和标准差，计算协方差阵
@Author: Jingdong Zhang
@DATE  : 2022/3/28
"""
import numpy as np


def gaia2ParaCorr2Cov(para):
    """由Gaia单历元的标准差及相关系数计算协方差阵

    Args:
        para (Series): getData()从数据库CRUD导出的数据行

    Returns:
        ndarray: 协方差阵
    """
    R = np.array([[1, para['g.RA_DEC_CORR']],
                  [para['g.RA_DEC_CORR'], 1]],
                 dtype=np.float64)
    D = np.diag(para[['g.RA_ERROR', 'g.DEC_ERROR']].astype(np.float64))
    return np.dot(np.dot(D, R), D)


def gaiaPmCorr2Cov(para):
    """由Gaia自行的标准差及相关系数计算协方差阵

    Args:
        para (Series): getData()从数据库CRUD导出的数据行

    Returns:
        ndarray: 协方差阵
    """
    R = np.array([[1, para['g.PMRA_PMDEC_CORR']],
                  [para['g.PMRA_PMDEC_CORR'], 1]],
                 dtype=np.float64)
    D = np.diag(para[['g.PMRA_ERROR', 'g.PMDEC_ERROR']].astype(np.float64))
    return np.dot(np.dot(D, R), D)


def gaia5ParaCorr2Cov(para):
    """由Gaia五参数的标准差及相关系数计算协方差阵

    Args:
        para (Series): getData()从数据库CRUD导出的数据行

    Returns:
        ndarray: 协方差阵
    """
    R = np.array([[1, para['g.RA_DEC_CORR'], para['g.RA_PI_CORR'], para['g.RA_PMRA_CORR'], para['g.RA_PMDEC_CORR']],
                  [para['g.RA_DEC_CORR'], 1, para['g.DEC_PI_CORR'], para['g.DEC_PMRA_CORR'], para['g.DEC_PMDEC_CORR']],
                  [para['g.RA_PI_CORR'], para['g.DEC_PI_CORR'], 1, para['g.PI_PMRA_CORR'], para['g.PI_PMDEC_CORR']],
                  [para['g.RA_PMRA_CORR'], para['g.DEC_PMRA_CORR'], para['g.PI_PMRA_CORR'], 1,
                   para['g.PMRA_PMDEC_CORR']],
                  [para['g.RA_PMDEC_CORR'], para['g.DEC_PMDEC_CORR'], para['g.PI_PMDEC_CORR'],
                   para['g.PMRA_PMDEC_CORR'], 1]], dtype=np.float64)
    D = np.diag(para[['g.RA_ERROR', 'g.DEC_ERROR', 'g.PI_ERROR', 'g.PMRA_ERROR', 'g.PMDEC_ERROR']].astype(np.float64))
    return np.dot(np.dot(D, R), D)


def gaia3ParaCorr2Cov(para):
    """由VLBI五参数的标准差及相关系数计算协方差阵

    Args:
        para (Series): getData()从数据库CRUD导出的数据行

    Returns:
        ndarray: 协方差阵
    """
    R = np.array([[1, para['g.PI_PMRA_CORR'], para['g.PI_PMDEC_CORR']],
                  [para['g.PI_PMRA_CORR'], 1, para['g.PMRA_PMDEC_CORR']],
                  [para['g.PI_PMDEC_CORR'], para['g.PMRA_PMDEC_CORR'], 1]], dtype=np.float64)
    D = np.diag(para[['g.PI_ERROR', 'g.PMRA_ERROR', 'g.PMDEC_ERROR']].astype(np.float64))
    return np.dot(np.dot(D, R), D)


def vlbi2ParaCorr2Cov(para):
    """由VLBI单历元的标准差构造协方差阵

    Args:
        para (Series): get_vlbi_obs()从数据库CRUD导出的数据行

    Returns:
        ndarray: 协方差阵
    """
    R = np.array([[1, 0],
                  [0, 1]])
    D = np.diag(para[['RA_ERROR', 'DEC_ERROR']].astype(np.float64))
    return np.dot(np.dot(D, R), D)


def vlbi5ParaCorr2Cov(para):
    """由VLBI五参数的标准差及相关系数计算协方差阵

    Args:
        para (Series): getData()从数据库CRUD导出的数据行

    Returns:
        ndarray: 协方差阵
    """
    R = np.array([[1, para['v.RA_DEC_CORR'], para['v.RA_PI_CORR'], para['v.RA_PMRA_CORR'], para['v.RA_PMDEC_CORR']],
                  [para['v.RA_DEC_CORR'], 1, para['v.DEC_PI_CORR'], para['v.DEC_PMRA_CORR'], para['v.DEC_PMDEC_CORR']],
                  [para['v.RA_PI_CORR'], para['v.DEC_PI_CORR'], 1, para['v.PI_PMRA_CORR'], para['v.PI_PMDEC_CORR']],
                  [para['v.RA_PMRA_CORR'], para['v.DEC_PMRA_CORR'], para['v.PI_PMRA_CORR'], 1,
                   para['v.PMRA_PMDEC_CORR']],
                  [para['v.RA_PMDEC_CORR'], para['v.DEC_PMDEC_CORR'], para['v.PI_PMDEC_CORR'],
                   para['v.PMRA_PMDEC_CORR'], 1]], dtype=np.float64)
    D = np.diag(para[['v.RA_ERROR', 'v.DEC_ERROR', 'v.PI_ERROR', 'v.PMRA_ERROR', 'v.PMDEC_ERROR']].astype(np.float64))
    return np.dot(np.dot(D, R), D)


def vlbi3ParaCorr2Cov(para):
    """由VLBI五参数的标准差及相关系数计算协方差阵

    Args:
        para (Series): getData()从数据库CRUD导出的数据行

    Returns:
        ndarray: 协方差阵
    """
    R = np.array([[1, para['v.PI_PMRA_CORR'], para['v.PI_PMDEC_CORR']],
                  [para['v.PI_PMRA_CORR'], 1, para['v.PMRA_PMDEC_CORR']],
                  [para['v.PI_PMDEC_CORR'], para['v.PMRA_PMDEC_CORR'], 1]], dtype=np.float64)
    D = np.diag(para[['v.PI_ERROR', 'v.PMRA_ERROR', 'v.PMDEC_ERROR']].astype(np.float64))
    return np.dot(np.dot(D, R), D)


if __name__ == '__main__':
    pass
