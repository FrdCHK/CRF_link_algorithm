"""
calculate the apparent motion caused by parallax & proper motion
@Author: Jingdong Zhang
@DATE  : 2022/8/26
"""
import numpy as np
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
from astropy import units as u

from units import mas2rad, rad2mas


def sunPosition(epoch, ephem='de440'):
    """calculate solar system barycenter position relative to Earth

    Args:
        epoch (float or array-like): J2000 epoch
        ephem (str): ephemeris version. Defaults to 'de440'.

    Returns:
        tuple: cartesian coordinate [x, y, z] in AU
    """
    t = Time(epoch, format='jyear')
    with solar_system_ephemeris.set(ephem):
        earth = get_body_barycentric('earth', t)
        x_au = earth.x.to(u.au)
        y_au = earth.y.to(u.au)
        z_au = earth.z.to(u.au)
        return -x_au.value, -y_au.value, -z_au.value


def ppm2pos(epoch0, ra0, dec0, pi, pmra, pmdec, epoch):
    """calculate the apparent motion caused by parallax & proper motion

    Args:
        epoch0 (float): J2000 reference epoch
        ra0 (float): RA in deg
        dec0 (float): DEC in deg
        pi (float): parallax in mas
        pmra (float): east direction proper motion in mas
        pmdec (float): north direction proper motion in mas
        epoch (float or array-like): J2000 target epoch

    Returns:
        tuple: east & north direction offset [RA*, DEC] in mas
    """
    
    # convert unit to rad
    ra0 = np.deg2rad(ra0)
    dec0 = np.deg2rad(dec0)
    pi = mas2rad(pi)
    pmra = mas2rad(pmra)
    pmdec = mas2rad(pmdec)

    # get sun position
    sun_X, sun_Y, sun_Z = sunPosition(epoch)

    # calculate position offset
    delta_epo = epoch - epoch0
    sRA = np.sin(ra0)
    cRA = np.cos(ra0)
    sDEC = np.sin(dec0)
    cDEC = np.cos(dec0)
    x_offset = pi * (sun_Y * cRA - sun_X * sRA) + delta_epo * pmra
    y_offset = pi * (sun_Z * cDEC - sun_X * cRA * sDEC - sun_Y * sRA * sDEC) + delta_epo * pmdec
    return rad2mas(x_offset), rad2mas(y_offset)


if __name__ == '__main__':
    # test code
    import matplotlib.pyplot as plt
    e = np.linspace(2020.3535, 2021.4346, 24)
    # a, b = ppm2pos(2016.0, 332.1696408540, 45.7425247583, 23.55882773476, -51.86451045170, 46.8518414858, e)
    a, b = ppm2pos(2016.0, 332.1696408540, 45.7425247583, 23.55882773476, 0, 0, e)
    plt.plot(a, b, '.-', c='#00BFFF')
    for i in range(len(a)):
        plt.text(a[i], b[i], f"{e[i]:.2f}")
    plt.xlabel('$\\Delta\\alpha$cos$\\delta$ (mas)')
    plt.ylabel('$\\Delta\\delta$ (mas)')
    plt.gca().set_aspect('equal')
    plt.show()
    pass
