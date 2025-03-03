"""
description
@Author: Jingdong Zhang
@DATE  : 2025/2/26
"""
import numpy as np

from ppm2pos import sunPosition
from units import mas2rad


class OriOnly:
    def __init__(self, input_data, spin_prior):
        self.source_name = input_data['SOURCE_NAME']
        self.ref_epo = input_data['REF_EPOCH']
        self.vlbi_5p = input_data['VLBI_5P']
        self.vlbi_2p = input_data['VLBI_2P']
        self.gaia = input_data['GAIA'].copy(deep=True)
        self.epo_2p = self.vlbi_2p.index.size

        self.spin = spin_prior.reshape(3, 1)

        self.niu = 0

        # 计算全过程中使用的三角函数必须一致
        self.cDEC = np.cos(input_data['GAIA']['DEC'])
        self.sDEC = np.sin(input_data['GAIA']['DEC'])
        self.cRA = np.cos(input_data['GAIA']['RA'])
        self.sRA = np.sin(input_data['GAIA']['RA'])
        self.gaia_spin()
        self.K_matrix = self.calcK()

        self.C_matrix = self.gaiaCov()
        self.inv_C = np.linalg.inv(self.C_matrix)
        self.gaia_value = np.array(self.gaia[['RA', 'DEC', 'PI', 'PMRA', 'PMDEC']]).T
        self.gaia_value[0] *= self.cDEC

        self.M = []
        self.V = []
        self.df = []
        for i, row in self.vlbi_5p.iterrows():
            self.vlbi5p(row)
            self.niu += 5

        for i, row in self.vlbi_2p.iterrows():
            self.vlbi2p(row)
            self.niu += 2

        tmp_1 = []
        tmp_2 = []
        for i in range(len(self.M)):
            tmp_1.append(self.M[i].T @ np.linalg.inv(self.V[i]) @ self.M[i])
            tmp_2.append(self.M[i].T @ np.linalg.inv(self.V[i]) @ self.df[i])
        self.sum_1 = tmp_1[0].copy()
        for i in range(1, len(tmp_1)):
            self.sum_1 += tmp_1[i]
        self.sum_2 = tmp_2[0].copy()
        for i in range(1, len(tmp_2)):
            self.sum_2 += tmp_2[i]
        self.H = self.inv_C + self.sum_1
        self.inv_H = np.linalg.inv(self.H)
        self.N_matrix = self.K_matrix.T @ self.inv_C @ (self.K_matrix - self.inv_H @ self.inv_C @ self.K_matrix)
        self.b_vector = self.K_matrix.T @ self.inv_C @ self.inv_H @ self.sum_2

        self.E = np.trace(self.N_matrix)
        # self.omega = np.trace(self.N_matrix[3:, 3:])

    def gaia_spin(self):
        tmp = np.zeros([2, 3])
        tmp[0, 0] = self.cRA * self.sDEC
        tmp[0, 1] = self.sRA * self.sDEC
        tmp[0, 2] = -self.cDEC
        tmp[1, 0] = -self.sRA
        tmp[1, 1] = self.cRA
        s = tmp @ self.spin
        self.gaia.loc[0, 'PMRA'] -= s[0, 0]
        self.gaia.loc[0, 'PMDEC'] -= s[1, 0]

    def gaiaCov(self):
        R = np.array([[1, self.gaia.at[0, 'RA_DEC_CORR'], self.gaia.at[0, 'RA_PI_CORR'],
                       self.gaia.at[0, 'RA_PMRA_CORR'],  self.gaia.at[0, 'RA_PMDEC_CORR']],
                      [self.gaia.at[0, 'RA_DEC_CORR'], 1, self.gaia.at[0, 'DEC_PI_CORR'],
                       self.gaia.at[0, 'DEC_PMRA_CORR'], self.gaia.at[0, 'DEC_PMDEC_CORR']],
                      [self.gaia.at[0, 'RA_PI_CORR'], self.gaia.at[0, 'DEC_PI_CORR'], 1,
                       self.gaia.at[0, 'PI_PMRA_CORR'], self.gaia.at[0, 'PI_PMDEC_CORR']],
                      [self.gaia.at[0, 'RA_PMRA_CORR'], self.gaia.at[0, 'DEC_PMRA_CORR'],
                       self.gaia.at[0, 'PI_PMRA_CORR'], 1, self.gaia.at[0, 'PMRA_PMDEC_CORR']],
                      [self.gaia.at[0, 'RA_PMDEC_CORR'], self.gaia.at[0, 'DEC_PMDEC_CORR'],
                       self.gaia.at[0, 'PI_PMDEC_CORR'], self.gaia.at[0, 'PMRA_PMDEC_CORR'], 1]], dtype=np.float64)
        D = np.diag(np.array(self.gaia[['RA_ERROR', 'DEC_ERROR', 'PI_ERROR', 'PMRA_ERROR', 'PMDEC_ERROR']].astype(np.float64)).flatten())
        return np.dot(np.dot(D, R), D)

    def calcK(self):
        tmp = np.zeros([5, 3])
        tmp[0, 0] = self.cRA * self.sDEC
        tmp[0, 1] = self.sRA * self.sDEC
        tmp[0, 2] = -self.cDEC
        tmp[1, 0] = -self.sRA
        tmp[1, 1] = self.cRA
        return tmp

    def vlbi5p(self, data):
        m = np.eye(5)
        m[0, 3] = data['EPOCH'] - self.ref_epo
        m[1, 4] = data['EPOCH'] - self.ref_epo
        self.M.append(m)

        R = np.array([[1, data['RA_DEC_CORR'], data['RA_PI_CORR'],
                       data['RA_PMRA_CORR'], data['RA_PMDEC_CORR']],
                      [data['RA_DEC_CORR'], 1, data['DEC_PI_CORR'],
                       data['DEC_PMRA_CORR'], data['DEC_PMDEC_CORR']],
                      [data['RA_PI_CORR'], data['DEC_PI_CORR'], 1,
                       data['PI_PMRA_CORR'], data['PI_PMDEC_CORR']],
                      [data['RA_PMRA_CORR'], data['DEC_PMRA_CORR'],
                       data['PI_PMRA_CORR'], 1, data['PMRA_PMDEC_CORR']],
                      [data['RA_PMDEC_CORR'], data['DEC_PMDEC_CORR'],
                       data['PI_PMDEC_CORR'], data['PMRA_PMDEC_CORR'], 1]], dtype=np.float64)
        D = np.diag(np.array(
            data[['RA_ERROR', 'DEC_ERROR', 'PI_ERROR', 'PMRA_ERROR', 'PMDEC_ERROR']].astype(np.float64)).flatten())
        v = np.dot(np.dot(D, R), D)
        self.V.append(v)

        vlbi_value = np.expand_dims(np.array(data[['RA', 'DEC', 'PI', 'PMRA', 'PMDEC']], dtype=np.float64), axis=1)
        vlbi_value[0, 0] *= self.cDEC
        df = vlbi_value - m @ self.gaia_value
        self.df.append(df)

    def vlbi2p(self, data):
        m = np.zeros((2, 5))
        sun_X, sun_Y, sun_Z = sunPosition(data['EPOCH'])
        m[0, 0] = 1
        m[1, 1] = 1
        m[0, 2] = sun_Y * self.cRA - sun_X * self.sRA
        m[0, 3] = data['EPOCH'] - self.ref_epo
        m[1, 2] = sun_Z * self.cDEC - sun_X * self.cRA * self.sDEC - sun_Y * self.sRA * self.sDEC
        m[1, 4] = data['EPOCH'] - self.ref_epo
        self.M.append(m)

        R = np.array([[1, data['RA_DEC_CORR']], [data['RA_DEC_CORR'], 1]], dtype=np.float64)
        D = np.diag(np.array(data[['RA_ERROR', 'DEC_ERROR']].astype(np.float64)).flatten())
        v = np.dot(np.dot(D, R), D)
        self.V.append(v)

        vlbi_value = np.expand_dims(np.array(data[['RA', 'DEC']], dtype=np.float64), axis=1)
        vlbi_value[0, 0] *= self.cDEC
        df = vlbi_value - m @ self.gaia_value
        self.df.append(df)

    def qScalar(self, x_vector):
        """计算卡方Q

        Args:
            x_vector (ndarray): 参数拟合结果

        Returns:
            该源所有数据的q之和, 该源每条数据的q
        """
        # y = np.linalg.inv(self.inv_C + self.sum_1) @ (self.inv_C @ self.K_matrix @ x_vector + self.sum_2)
        # tmp = np.zeros((1, 1))
        # for i in range(len(self.M)):
        #     tmp = tmp + (self.df[i] - self.M[i] @ y).T @ self.V[i] @ (self.df[i] - self.M[i] @ y)
        # q = (y - self.K_matrix @ x_vector).T @ self.inv_C @ (y - self.K_matrix @ x_vector) + tmp
        # return q[0, 0]

        q_list = []
        for i in range(len(self.M)):
            D = self.V[i] + self.M[i] @ self.C_matrix @ self.M[i].T
            tmp = self.df[i] - self.M[i] @ self.K_matrix @ x_vector
            q_list.append((tmp.T @ np.linalg.inv(D) @ tmp)[0, 0])
        return np.sum(q_list), np.array(q_list)


def find_prior(d, v):
    """
    查找DataFrame中GMIN<v<GMAX的行，并返回对应行的X, Y, Z列的值(微角秒转为弧度)。

    参数:
    d (pd.DataFrame): 包含GMAX, GMIN, X, Y, Z五列的DataFrame。
    v (float): 要检查的浮点数。

    返回:
    np.ndarray: 包含X, Y, Z列值的ndarray，如果没有找到符合条件的行，则返回空ndarray。
    """
    # 遍历DataFrame的每一行
    for _, row in d.iterrows():
        if row['GMIN'] < v < row['GMAX']:
            return mas2rad(np.array([row['X'], row['Y'], row['Z']])/1e3)


def solveAll(data, spin_prior):
    niu_list = []
    e_list = []
    n_sum = np.zeros([3, 3])
    b_sum = np.zeros([3, 1])
    data_list = []
    for row in data:
        tmp_data_i = OriOnly(row, find_prior(spin_prior, row['GAIA'].at[0, 'GMAG']))
        data_list.append(tmp_data_i)
        niu_list.append(tmp_data_i.niu)
        e_list.append(tmp_data_i.E)
        n_sum = n_sum + tmp_data_i.N_matrix
        b_sum = b_sum + tmp_data_i.b_vector

    x_result = np.linalg.solve(n_sum, b_sum)
    cov = np.linalg.inv(n_sum)

    # Q
    q_list = []
    for item in data_list:
        q_list.append(item.qScalar(x_result))
    return x_result, cov, q_list, np.array(niu_list), np.array(e_list)


def data_screen(data):
    data_return = []
    for item in data:
        item['VLBI_5P'] = item['VLBI_5P'].dropna(subset=['RA_ERROR', 'DEC_ERROR'])
        if (item['VLBI_5P'].shape[0] > 0) | (item['VLBI_2P'].shape[0] > 0):
            data_return.append(item)
    return data_return
