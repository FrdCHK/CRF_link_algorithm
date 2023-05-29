"""
CRF link algorithm
@Author: Jingdong Zhang
@DATE  : 2022/3/24
"""
import numpy as np

from corr2Cov import gaia2ParaCorr2Cov, gaia5ParaCorr2Cov, vlbi5ParaCorr2Cov, gaia3ParaCorr2Cov, vlbi3ParaCorr2Cov, gaiaPmCorr2Cov


"""
Individual-position solution
"""
class Indi:
    def __init__(self, para, ref_epoch):
        """init

        Args:
            para (Dataframe): input parameters
            ref_epoch (float): reference epoch of CRF link
        """
        self.source_name = para['SOURCE_NAME']
        self.vlbi = para[['v.RA', 'v.DEC']].astype(np.float64)
        self.vlbi.index = ['RA', 'DEC']
        self.gaia = para[['g.RA', 'g.DEC']].astype(np.float64)
        self.gaia.index = ['RA', 'DEC']
        self.delta_epoch = para['EPOCH']-ref_epoch
        self.V_matrix = np.diag(para[['v.RA_ERROR', 'v.DEC_ERROR']].astype(np.float64))**2
        self.C_matrix = gaia2ParaCorr2Cov(para)
        self.delta_f = None
        self.deltaF()
        self.K_matrix = None
        self.kMatrix()
        self.D_matrix = None
        self.dMatrix()
        self.inv_D = np.linalg.inv(self.D_matrix)
        self.N_matrix = None
        self.nMatrix()
        self.b_vector = None
        self.bVector()
        self.E = np.trace(self.N_matrix[:3, :3])
        self.omega = np.trace(self.N_matrix[3:, 3:])

    def deltaF(self):
        self.delta_f = np.array([[(self.vlbi['RA'] - self.gaia['RA']) * np.cos(self.gaia['DEC'])],
                                 [self.vlbi['DEC'] - self.gaia['DEC']]])

    def kMatrix(self):
        tmp = np.zeros([2, 6])
        tmp[0, 0] = np.cos(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[0, 1] = np.sin(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[0, 2] = -np.cos(self.gaia['DEC'])
        tmp[0, 3] = np.cos(self.gaia['RA']) * np.sin(self.gaia['DEC'])*self.delta_epoch
        tmp[0, 4] = np.sin(self.gaia['RA']) * np.sin(self.gaia['DEC'])*self.delta_epoch
        tmp[0, 5] = -np.cos(self.gaia['DEC'])*self.delta_epoch
        tmp[1, 0] = -np.sin(self.gaia['RA'])
        tmp[1, 1] = np.cos(self.gaia['RA'])
        tmp[1, 3] = -np.sin(self.gaia['RA'])*self.delta_epoch
        tmp[1, 4] = np.cos(self.gaia['RA'])*self.delta_epoch
        self.K_matrix = tmp

    def dMatrix(self):
        self.D_matrix = self.V_matrix + self.C_matrix

    def nMatrix(self):
        self.N_matrix = np.dot(np.dot(self.K_matrix.T, self.inv_D), self.K_matrix)

    def bVector(self):
        self.b_vector = np.dot(np.dot(self.K_matrix.T, self.inv_D), self.delta_f)

    def qScalar(self, x_vector):
        """calculate the chi2

        Args:
            x_vector (ndarray): parameter estimatating result

        Returns:
            float: chi2
        """
        R = self.delta_f - np.dot(self.K_matrix, x_vector)
        return np.dot(np.dot(R.T, self.inv_D), R)[0, 0]


"""
Proper-motion solution
"""
class Pm:
    def __init__(self, para):
        """init

        Args:
            para (Dataframe): input parameters
        """
        self.source_name = para['SOURCE_NAME']
        self.vlbi = para[['v.PMRA', 'v.PMDEC']].astype(np.float64)
        self.vlbi.index = ['PMRA', 'PMDEC']
        self.gaia = para[['g.RA', 'g.DEC', 'g.PMRA', 'g.PMDEC']].astype(np.float64)
        self.gaia.index = ['RA', 'DEC', 'PMRA', 'PMDEC']
        self.V_matrix = np.diag(para[['v.PMRA_ERROR', 'v.PMDEC_ERROR']].astype(np.float64))**2
        self.C_matrix = gaiaPmCorr2Cov(para)
        self.delta_f = None
        self.deltaF()
        self.K_matrix = None
        self.kMatrix()
        self.D_matrix = None
        self.dMatrix()
        self.inv_D = np.linalg.inv(self.D_matrix)
        self.N_matrix = None
        self.nMatrix()
        self.b_vector = None
        self.bVector()
        self.E = 0
        self.omega = np.trace(self.N_matrix[3:, 3:])

    def deltaF(self):
        self.delta_f = np.array([[(self.vlbi['PMRA'] - self.gaia['PMRA']) * np.cos(self.gaia['DEC'])],
                                 [self.vlbi['PMDEC'] - self.gaia['PMDEC']]])

    def kMatrix(self):
        tmp = np.zeros([2, 6])
        tmp[0, 3] = np.cos(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[0, 4] = np.sin(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[0, 5] = -np.cos(self.gaia['DEC'])
        tmp[1, 3] = -np.sin(self.gaia['RA'])
        tmp[1, 4] = np.cos(self.gaia['RA'])
        self.K_matrix = tmp

    def dMatrix(self):
        self.D_matrix = self.V_matrix + self.C_matrix

    def nMatrix(self):
        self.N_matrix = np.dot(np.dot(self.K_matrix.T, self.inv_D), self.K_matrix)

    def bVector(self):
        self.b_vector = np.dot(np.dot(self.K_matrix.T, self.inv_D), self.delta_f)

    def qScalar(self, x_vector):
        """calculate the chi2

        Args:
            x_vector (ndarray): parameter estimatating result

        Returns:
            float: chi2
        """
        R = self.delta_f - np.dot(self.K_matrix, x_vector)
        return np.dot(np.dot(R.T, self.inv_D), R)[0, 0]


class Lin5Para:
    def __init__(self, para):
        """init

        Args:
            para (Dataframe): input parameters
        """
        self.source_name = para['SOURCE_NAME']
        self.vlbi = para[['v.EPOCH', 'v.RA', 'v.DEC', 'v.PI', 'v.PMRA', 'v.PMDEC']].astype(np.float64)
        self.vlbi.index = ['EPOCH', 'RA', 'DEC', 'PI', 'PMRA', 'PMDEC']
        self.gaia = para[['g.EPOCH', 'g.RA', 'g.DEC', 'g.PI', 'g.PMRA', 'g.PMDEC']].astype(np.float64)
        self.gaia.index = ['EPOCH', 'RA', 'DEC', 'PI', 'PMRA', 'PMDEC']
        self.V_matrix = vlbi5ParaCorr2Cov(para)
        self.C_matrix = gaia5ParaCorr2Cov(para)
        self.M_matrix = None
        self.mMatrix()
        self.delta_f = None
        self.deltaF()
        self.K_matrix = None
        self.kMatrix()
        self.D_matrix = None
        self.dMatrix()
        self.inv_D = np.linalg.inv(self.D_matrix)
        self.N_matrix = None
        self.nMatrix()
        self.b_vector = None
        self.bVector()
        self.E = np.trace(self.N_matrix[:3, :3])
        self.omega = np.trace(self.N_matrix[3:, 3:])

    def mMatrix(self):
        tmp = np.eye(5)
        tmp[0, 3] = self.vlbi['EPOCH'] - self.gaia['EPOCH']
        tmp[1, 4] = self.vlbi['EPOCH'] - self.gaia['EPOCH']
        self.M_matrix = tmp
    
    def deltaF(self):
        delta_epoch = self.vlbi['EPOCH'] - self.gaia['EPOCH']
        g_propagated = self.gaia.copy()
        g_propagated.loc['EPOCH'] = self.vlbi['EPOCH']
        g_propagated.loc['RA'] = self.gaia['RA'] + (self.gaia['PMRA'] / np.cos(self.gaia['DEC'])) * delta_epoch
        g_propagated.loc['DEC'] = self.gaia['DEC'] + self.gaia['PMDEC'] * delta_epoch
        self.delta_f = np.array([[(self.vlbi['RA'] - g_propagated['RA']) * np.cos(self.gaia['DEC'])],
                                 [self.vlbi['DEC'] - g_propagated['DEC']],
                                 [self.vlbi['PI'] - g_propagated['PI']],
                                 [(self.vlbi['PMRA'] - g_propagated['PMRA'])],
                                 [self.vlbi['PMDEC'] - g_propagated['PMDEC']]])
    
    def kMatrix(self):
        tmp = np.zeros([5, 6])
        tmp[0, 0] = np.cos(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[0, 1] = np.sin(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[0, 2] = -np.cos(self.gaia['DEC'])
        tmp[1, 0] = -np.sin(self.gaia['RA'])
        tmp[1, 1] = np.cos(self.gaia['RA'])
        tmp[3, 3] = np.cos(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[3, 4] = np.sin(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[3, 5] = -np.cos(self.gaia['DEC'])
        tmp[4, 3] = -np.sin(self.gaia['RA'])
        tmp[4, 4] = np.cos(self.gaia['RA'])
        self.K_matrix = tmp

    def dMatrix(self):
        self.D_matrix = self.V_matrix + np.dot(np.dot(self.M_matrix, self.C_matrix), self.M_matrix.T)

    def nMatrix(self):
        self.N_matrix = np.dot(np.dot(np.dot(np.dot(self.K_matrix.T, self.M_matrix.T), self.inv_D), self.M_matrix), self.K_matrix)

    def bVector(self):
        self.b_vector = np.dot(np.dot(np.dot(self.K_matrix.T, self.M_matrix.T), self.inv_D), self.delta_f)

    def qScalar(self, x_vector):
        """calculate the chi2

        Args:
            x_vector (ndarray): parameter estimatating result

        Returns:
            float: chi2
        """
        R = self.delta_f - np.dot(np.dot(self.M_matrix, self.K_matrix), x_vector)
        return np.dot(np.dot(R.T, self.inv_D), R)[0, 0]


class Lin3Para:
    def __init__(self, para):
        """init

        Args:
            para (Dataframe): input parameters
        """
        self.source_name = para['SOURCE_NAME']
        self.vlbi = para[['v.EPOCH', 'v.RA', 'v.DEC', 'v.PI', 'v.PMRA', 'v.PMDEC']].astype(np.float64)
        self.vlbi.index = ['EPOCH', 'RA', 'DEC', 'PI', 'PMRA', 'PMDEC']
        self.gaia = para[['g.EPOCH', 'g.RA', 'g.DEC', 'g.PI', 'g.PMRA', 'g.PMDEC']].astype(np.float64)
        self.gaia.index = ['EPOCH', 'RA', 'DEC', 'PI', 'PMRA', 'PMDEC']
        self.V_matrix = vlbi3ParaCorr2Cov(para)
        self.C_matrix = gaia3ParaCorr2Cov(para)
        self.delta_f = None
        self.deltaF()
        self.K_matrix = None
        self.kMatrix()
        self.D_matrix = None
        self.dMatrix()
        self.inv_D = np.linalg.inv(self.D_matrix)
        self.N_matrix = None
        self.nMatrix()
        self.b_vector = None
        self.bVector()
        self.E = np.trace(self.N_matrix[:3, :3])
        self.omega = np.trace(self.N_matrix[3:, 3:])

    def deltaF(self):
        self.delta_f = np.array([[self.vlbi['PI'] - self.gaia['PI']],
                                 [(self.vlbi['PMRA'] - self.gaia['PMRA'])],
                                 [self.vlbi['PMDEC'] - self.gaia['PMDEC']]])

    def kMatrix(self):
        tmp = np.zeros([3, 6])
        tmp[1, 3] = np.cos(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[1, 4] = np.sin(self.gaia['RA']) * np.sin(self.gaia['DEC'])
        tmp[1, 5] = -np.cos(self.gaia['DEC'])
        tmp[2, 3] = -np.sin(self.gaia['RA'])
        tmp[2, 4] = np.cos(self.gaia['RA'])
        self.K_matrix = tmp

    def dMatrix(self):
        self.D_matrix = self.V_matrix + self.C_matrix

    def nMatrix(self):
        self.N_matrix = np.dot(np.dot(self.K_matrix.T, self.inv_D), self.K_matrix)

    def bVector(self):
        self.b_vector = np.dot(np.dot(self.K_matrix.T, self.inv_D), self.delta_f)

    def qScalar(self, x_vector):
        """calculate the chi2

        Args:
            x_vector (ndarray): parameter estimatating result

        Returns:
            float: chi2
        """
        R = self.delta_f - np.dot(self.K_matrix, x_vector)
        return np.dot(np.dot(R.T, self.inv_D), R)[0, 0]


def solveMix(five_para=None, two_para=None, ref_epoch=None, pm=None):
    """Hybrid solution
    if two_para is not None, ref_epoch must be passed.

    Args:
        five_para (Dataframe): five-parameter data. Defaults to None.
        two_para (Dataframe): two-parameter data. Defaults to None.
        ref_epoch (float): reference epoch of CRF link. Defaults to None.
        pm (Dataframe): proper motion data. Defaults to None.

    Returns:
        tuple[ndarray, ndarray, ndarray, ndarray, ndarray, ndarray]: [parameter estimating result, covariance matrix, Q, niu, E, Omega]
    """
    niu_list = []
    e_list = []
    o_list = []
    n_sum_5p = np.zeros([6, 6])
    b_sum_5p = np.zeros([6, 1])
    data_list = []
    if five_para is not None:
        for _, row in five_para.iterrows():
            if np.isnan(row['v.RA_ERROR']) | np.isnan(row['v.DEC_ERROR']):
                tmp_data_i = Lin3Para(row)
                niu_list.append(3)
            else:
                tmp_data_i = Lin5Para(row)
                niu_list.append(5)
            data_list.append(tmp_data_i)
            e_list.append(tmp_data_i.E)
            o_list.append(tmp_data_i.omega)
            n_sum_5p = n_sum_5p + tmp_data_i.N_matrix
            b_sum_5p = b_sum_5p + tmp_data_i.b_vector

    n_sum_2p = np.zeros([6, 6])
    b_sum_2p = np.zeros([6, 1])
    if two_para is not None:
        for _, row in two_para.iterrows():
            tmp_data_i = Indi(row, ref_epoch)
            niu_list.append(2)
            e_list.append(tmp_data_i.E)
            o_list.append(tmp_data_i.omega)
            data_list.append(tmp_data_i)
            n_sum_2p = n_sum_2p + tmp_data_i.N_matrix
            b_sum_2p = b_sum_2p + tmp_data_i.b_vector

    n_sum_pm = np.zeros([6, 6])
    b_sum_pm = np.zeros([6, 1])
    if pm is not None:
        for _, row in pm.iterrows():
            tmp_data_i = Pm(row)
            niu_list.append(2)
            e_list.append(tmp_data_i.E)
            o_list.append(tmp_data_i.omega)
            data_list.append(tmp_data_i)
            n_sum_pm = n_sum_pm + tmp_data_i.N_matrix
            b_sum_pm = b_sum_pm + tmp_data_i.b_vector

    n_sum = n_sum_5p+n_sum_2p+n_sum_pm
    b_sum = b_sum_5p+b_sum_2p+b_sum_pm
    x_result = np.linalg.solve(n_sum, b_sum)
    cov = np.linalg.inv(n_sum)

    # Q
    q_list = []
    for item in data_list:
        q_list.append(item.qScalar(x_result))
    return x_result, cov, np.array(q_list), np.array(niu_list), np.array(e_list), np.array(o_list)


if __name__ == '__main__':
    pass
