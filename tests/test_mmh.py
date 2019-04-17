#!/usr/bin/env python
import numpy as np
import sys
sys.path.insert(0, '..')

def test_init():
  from pybind.mmh import McMillanHe
  mmh = McMillanHe()

def test_get_lbox():
  from pybind.mmh import McMillanHe
  mmh = McMillanHe()
  lbox = mmh.get_lbox()
  assert np.isclose(11.3303267, lbox)

def test_set_lbox():
  from pybind.mmh import McMillanHe
  mmh = McMillanHe()
  lbox = 22.66
  mmh.set_lbox(lbox)
  assert np.isclose(lbox, mmh.get_lbox())

def test_get_a1():
  from pybind.mmh import McMillanHe
  mmh = McMillanHe()
  a1 = mmh.get_a1()
  assert np.isclose(2.6, a1)

def test_set_a1():
  from pybind.mmh import McMillanHe
  mmh = McMillanHe()
  a1 = 5.0
  mmh.set_a1(a1)
  assert np.isclose(a1, mmh.get_a1())

def test_wfval():
  from pybind.mmh import McMillanHe
  mmh = McMillanHe()
  pos = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 2.0, 3.0],
  ])
  assert np.isclose(0.850431168292, mmh.wfval(pos))
