#!/usr/bin/env python
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

from qharv.inspect import volumetric, axes_pos
from qharv.inspect.grsk import sk2gr, gr2sk
from qharv.plantation import kyrt, sugar

def init_pos(natom):
  from qharv.inspect import crystal
  from ase import Atoms
  from ase.build import make_supercell
  rho = 2.2e22*1e6/1e30  # atoms/A^3
  lbox = (natom/rho)**(1./3)
  nxf = (natom/4.)**(1./3)
  nx = int(round(nxf))
  if not np.isclose(nx, nxf):
    raise RuntimeError('natom=%d nxf=%3.2f!=%d' % (natom, nxf, nx))
  # create FCC crystal
  alat = lbox/nx
  axes0 = alat/2*(np.ones(3)-np.eye(3))
  tmat = nx*(np.ones(3)-2*np.eye(3))
  s0 = Atoms('H', cell=axes0, positions=[[0, 0, 0]], pbc=[1, 1, 1])
  s1 = make_supercell(s0, tmat)
  pos = s1.get_positions()
  axes = s1.get_cell()
  # check density
  rho1 = natom/axes_pos.volume(axes)
  if not np.isclose(rho, rho1):
    raise RuntimeError('supercell density is wrong')
  # save/view crystal
  fig, ax = volumetric.figax3d()
  crystal.draw_cell(ax, axes)
  crystal.draw_atoms(ax, pos)
  plt.show()
  return pos

def mean_and_error(t1, t2, nobs):
  """Calculate mean and error given accumulated first and second moments.

  Args:
    t1 (np.array): accumulated first moment \sum x_i
    t2 (np.array): accumulated second moment \sum x_i^2
    nobs (int): number of accumulated observations
  Return:
    (np.array, np.array): (mean, error)
  """
  tm = t1/nobs
  te = (t2-tm**2)**0.5/((nobs-1)*nobs)**0.5
  return tm, te

def get_sofk(posl, lbox, nsh):
  """Calculate static structure factor S(k) from sampled configurations.

  Args:
    posl (list): a list of sampled configurations
    lbox (float): side length of cubic box
    nsh (int): number of kshells to use
  Return:
    (np.array,)*3: (kvecs, skm, ske), (kvetors, S(k) mean, S(k) error)
  """
  from forlib.grsk import calc_sofk
  nvecs = axes_pos.cubic_pos(nsh)
  kvecs = 2*np.pi/lbox*nvecs[1:]  # skip k=0
  kmags = np.linalg.norm(kvecs, axis=-1)
  nk = len(kmags)
  skm = np.zeros(nk)
  sks1 = np.zeros(nk)
  sks2 = np.zeros(nk)
  nobs = 0
  for iblock, pos1 in enumerate(posl):
    if iblock > nequil:
      sk1 = calc_sofk(kvecs, pos1, lbox)
      sks1 += sk1
      sks2 += sk1**2
      nobs += 1
  myskm, myske = mean_and_error(sks1, sks2, nobs)
  return kvecs, myskm, myske

def get_gofr(posl, lbox, ngr):
  """Calculate pair correlation g(r) from sampled configurations.

  Args:
    posl (list): a list of sampled configurations
    lbox (float): side length of cubic box
    ngr (int): number linear grid points from [0, lbox/2)
  Return:
    (np.array,)*3: (myr, grm, gre), (r grid, g(r) mean, g(r) error)
  """
  from forlib.grsk import calc_gofr
  dr = lbox/2./ngr
  myr = np.arange(0, lbox/2-1e-10, dr)
  assert len(myr) == ngr
  myr += dr/2.
  grnorm = np.concatenate([
    [4*np.pi/3*(dr/2.)**3],
    np.diff(4*np.pi/3*myr**3)
  ], axis=0)/(2*lbox**3/natom/(natom-1))
  grs1 = np.zeros(ngr)
  grs2 = np.zeros(ngr)
  nobs = 0
  for iblock, pos1 in enumerate(posl):
    if iblock > nequil:
      gr1 = calc_gofr(pos1, lbox, ngr)
      grs1 += gr1
      grs2 += gr1**2
      nobs += 1
  mygrm, mygre = mean_and_error(grs1, grs2, nobs)
  mygrm /= grnorm
  mygre /= grnorm
  return myr, mygrm, mygre

if __name__ == '__main__':
  from pybind.mmh import McMillanHe
  mmh = McMillanHe()
  lbox = mmh.get_lbox()
  # simulation parameters
  #  initial configuration
  natom = 32
  #natom = 108
  nx = int(round((natom/4.)**(1./3)))
  lbox *= nx/2.
  mmh.set_lbox(lbox)
  pdat = 'pos%d.dat' % natom
  if os.path.isfile(pdat):
    pos = np.loadtxt(pdat)
  else:  # require ase
    pos = init_pos(natom)
    np.savetxt(pdat, pos)
  #  VMC parameters
  nsubstep = 16  # number of steps per block
  #   N108 solid
  a1 = 3.1
  nequil = 64    # number of equilibration blocks
  nblock = 4096  # number of blocks to run
  tau = 0.01     # time step
  #   N108 liquid
  a1 = 2.1
  #a1 = 3.5
  tau = 0.1      # time step
  nblock = 4096  # number of blocks to run
  #   N32 liquid
  a1 = 2.1
  tau = 0.5      # time step
  nblock = 256  # number of blocks to run

  mmh.set_a1(a1)
  # sample and cache configuartions
  view = False  # view sampled configurations in 3D (can be slow)
  prefix = 'n%d-a%3.2f-tau%3.2f-nb%d' % (natom, a1, tau, nblock)
  fpos = 'cache/%s_all-pos.h5' % prefix
  if not os.path.isfile(fpos):
    sugar.mkdir(os.path.dirname(fpos))
    from qharv.reel import config_h5
    fp = config_h5.open_write(fpos)
    for iblock in range(nblock):
      path = 'pos%06d' % iblock
      pos1 = mmh.diffuse(pos, nsubstep, tau)
      config_h5.save_vec(pos1, fp, fp.root, path)
      pos = pos1
    fp.close()
    print('acceptance: ', mmh.get_acc())
  # read sampled configurations
  posl = []
  fp = h5py.File(fpos, 'r')
  snaps = fp.keys()
  snaps.sort()
  for snap in snaps:
    pos1 = fp[snap][()]
    posl.append(pos1)
  fp.close()

  if view:
    # view configurations
    fig, ax = volumetric.figax3d()
    for iblock, pos1 in enumerate(posl):
      volumetric.color_scatter(ax, pos1)
    plt.show()

  # calculate observables
  # sofk
  nsh = 6
  kvecs, skm, ske = get_sofk(posl, lbox, nsh)
  import static_correlation as sc  # shell average
  uk, uskm, uske = sc.shavg(kvecs, skm, ske)
  #  gofr
  ngr = 28
  myr, grm, gre = get_gofr(posl, lbox, ngr)

  # plot together
  nskip = 8  # skip a few bad points
  gr1 = sk2gr(myr[nskip:], uk, uskm, natom/lbox**3)
  sk1 = gr2sk(uk[nskip:], myr, grm, natom/lbox**3)
  fig, axl = plt.subplots(1, 2)
  ax = axl[0]
  ax.set_xlabel(r'r ($\AA$)')
  ax.set_ylabel('g(r)')
  ax.errorbar(myr[:-1], grm[:-1], gre[:-1], ls='', marker='x')
  ax.plot(myr[nskip:], gr1, label='FT[sk]')
  ax.legend()
  ax = axl[1]
  ax.set_xlabel(r'k ($\AA^{-1}$)')
  ax.set_ylabel('S(k)')
  ax.errorbar(uk, uskm, uske, ls='', marker='x')
  ax.plot(uk[nskip:], sk1, label='FT[gr]')
  ax.legend()
  kyrt.yright(ax)
  fig.subplots_adjust(wspace=0.01)
  plt.show()
# end __main__
