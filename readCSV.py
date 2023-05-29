"""
read data from csv file directly
@Author: Jingdong Zhang
@DATE  : 2022/12/13
"""
import numpy as np
import pandas as pd

from units import mas2rad


def readCSV(fname):
    data_in = pd.read_csv(fname, comment="#")
    df = pd.DataFrame(columns=['SOURCE_NAME', 'v.EPOCH', 'v.RA', 'v.RA_ERROR', 'v.DEC', 'v.DEC_ERROR', 'v.PI',
                               'v.PI_ERROR', 'v.PMRA', 'v.PMRA_ERROR', 'v.PMDEC', 'v.PMDEC_ERROR', 'v.RA_DEC_CORR',
                               'v.RA_PI_CORR',
                               'v.RA_PMRA_CORR', 'v.RA_PMDEC_CORR', 'v.DEC_PI_CORR', 'v.DEC_PMRA_CORR',
                               'v.DEC_PMDEC_CORR', 'v.PI_PMRA_CORR', 'v.PI_PMDEC_CORR', 'v.PMRA_PMDEC_CORR',
                               'DESIGNATION', 'g.EPOCH', 'g.RA', 'g.RA_ERROR', 'g.DEC', 'g.DEC_ERROR', 'g.PI',
                               'g.PI_ERROR',
                               'g.PMRA', 'g.PMRA_ERROR', 'g.PMDEC', 'g.PMDEC_ERROR', 'g.RA_DEC_CORR', 'g.RA_PI_CORR',
                               'g.RA_PMRA_CORR', 'g.RA_PMDEC_CORR', 'g.DEC_PI_CORR', 'g.DEC_PMRA_CORR',
                               'g.DEC_PMDEC_CORR', 'g.PI_PMRA_CORR', 'g.PI_PMDEC_CORR', 'g.PMRA_PMDEC_CORR'])
    df[['SOURCE_NAME', 'v.EPOCH', 'v.RA', 'v.RA_ERROR', 'v.DEC', 'v.DEC_ERROR', 'v.PI',
        'v.PI_ERROR', 'v.PMRA', 'v.PMRA_ERROR', 'v.PMDEC', 'v.PMDEC_ERROR',
        'DESIGNATION', 'g.EPOCH', 'g.RA', 'g.RA_ERROR', 'g.DEC', 'g.DEC_ERROR', 'g.PI',
        'g.PI_ERROR',
        'g.PMRA', 'g.PMRA_ERROR', 'g.PMDEC', 'g.PMDEC_ERROR', 'g.RA_DEC_CORR', 'g.RA_PI_CORR',
        'g.RA_PMRA_CORR', 'g.RA_PMDEC_CORR', 'g.DEC_PI_CORR', 'g.DEC_PMRA_CORR',
        'g.DEC_PMDEC_CORR', 'g.PI_PMRA_CORR', 'g.PI_PMDEC_CORR', 'g.PMRA_PMDEC_CORR']] = data_in[
        ['SOURCE_NAME', 'EPOCH_VLBI', 'RA_VLBI', 'RA_ERR_VLBI', 'DEC_VLBI', 'DEC_ERR_VLBI', 'PI_VLBI',
         'PI_ERR_VLBI', 'PMRA_VLBI', 'PMRA_ERR_VLBI', 'PMDEC_VLBI', 'PMDEC_ERR_VLBI',
         'DESIGNATION', 'EPOCH', 'RA', 'RA_ERROR', 'DEC', 'DEC_ERROR', 'PI', 'PI_ERROR',
         'PMRA', 'PMRA_ERROR', 'PMDEC', 'PMDEC_ERROR', 'RA_DEC_CORR', 'RA_PI_CORR',
         'RA_PMRA_CORR', 'RA_PMDEC_CORR', 'DEC_PI_CORR', 'DEC_PMRA_CORR',
         'DEC_PMDEC_CORR', 'PI_PMRA_CORR', 'PI_PMDEC_CORR', 'PMRA_PMDEC_CORR']]
    # to rad
    df[['v.RA', 'v.DEC', 'g.RA', 'g.DEC']] = np.deg2rad(df[['v.RA', 'v.DEC', 'g.RA', 'g.DEC']])
    df[['v.RA_ERROR', 'v.DEC_ERROR', 'v.PI', 'v.PI_ERROR', 'v.PMRA', 'v.PMRA_ERROR', 'v.PMDEC',
        'v.PMDEC_ERROR', 'g.RA_ERROR', 'g.DEC_ERROR', 'g.PI', 'g.PI_ERROR', 'g.PMRA', 'g.PMRA_ERROR',
        'g.PMDEC', 'g.PMDEC_ERROR']] = mas2rad(df[['v.RA_ERROR', 'v.DEC_ERROR', 'v.PI', 'v.PI_ERROR', 'v.PMRA',
                                                   'v.PMRA_ERROR', 'v.PMDEC', 'v.PMDEC_ERROR', 'g.RA_ERROR',
                                                   'g.DEC_ERROR', 'g.PI', 'g.PI_ERROR', 'g.PMRA', 'g.PMRA_ERROR',
                                                   'g.PMDEC', 'g.PMDEC_ERROR']])
    # VLBI corr init
    df[['v.RA_DEC_CORR', 'v.RA_PI_CORR', 'v.RA_PMRA_CORR', 'v.RA_PMDEC_CORR', 'v.DEC_PI_CORR', 'v.DEC_PMRA_CORR',
        'v.DEC_PMDEC_CORR', 'v.PI_PMRA_CORR', 'v.PI_PMDEC_CORR', 'v.PMRA_PMDEC_CORR']] = 0
    return df


if __name__ == '__main__':
    pass
