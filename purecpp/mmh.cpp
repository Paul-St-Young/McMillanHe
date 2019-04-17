#include <mmh.h>
#include <iostream>
using namespace std;

void McMillanHe::pos_in_box(Matrix& pos)
{
  int natom = pos.rows();
  int ndim = pos.cols();
  for (int iatom=0; iatom<natom; iatom++)
  {
    for (int idim=0; idim<ndim; idim++)
    {
      pos(iatom, idim) -= _lbox*std::round(pos(iatom, idim)/_lbox);
    }
  }
}

McMillanHe::Vector
McMillanHe::displacement(const Vector& r1, const Vector& r2)
{
  int ndim = r1.size();
  int nint;
  Vector dr = r1-r2;
  for (int idim=0; idim<ndim; idim++)
  {
    nint = std::round(dr(idim)/_lbox);
    dr(idim) -= nint*_lbox;
  }
  return dr;
}

double McMillanHe::lnf(const double r)
{
  return -std::pow(_a1/r, _a2);
}

double McMillanHe::fval(const double r)
{
  return std::exp(lnf(r));
}

double McMillanHe::lnwf(const Matrix& pos)
{
  int natom = pos.rows();
  double dist;  // temporary pair distance
  double val = 0.0;
  for (int i=0; i<natom; i++)
  {
    for (int j=i+1; j<natom; j++)
    {
      dist = displacement(pos.row(i), pos.row(j)).norm();
      val += lnf(dist);
    }
  }
  return val;
}

double McMillanHe::wfval(const Matrix& pos)
{
  return std::exp(lnwf(pos));
}

double McMillanHe::diff_lnwf(const Matrix& pos, const Vector& move, const int i)
{
  int natom = pos.rows();
  double diff = 0.0;
  double rij, rij1;
  Vector curpos = pos.row(i);
  Vector newpos = curpos + move;
  for (int j=0; j<natom; j++)
  {
    if (i==j) continue;
    rij = displacement(curpos, pos.row(j)).norm();
    rij1 = displacement(newpos, pos.row(j)).norm();
    diff += pow(_a1/rij, _a2) - pow(_a1/rij1, _a2);
  }
  return diff;
}

double McMillanHe::ratio(const Matrix& pos, const Vector& move, const int i)
{
  return std::exp(diff_lnwf(pos, move, i));
}

McMillanHe::Vector McMillanHe::grad_lnwf(const Matrix& pos, const int i)
{
  int natom = pos.rows();
  int ndim = pos.cols();
  Vector grad(ndim);
  Vector disp(ndim); // temporary pair displacement
  double dist;  // temporary pair distance
  for (int j=0; j<natom; j++)
  {
    if (j==i) continue;
    disp = displacement(pos.row(i), pos.row(j));
    dist = disp.norm();
    grad += lnf(dist)/std::pow(dist, 2) * disp;
  }
  return -_a2*grad;
}

double McMillanHe::lap_lnwf(const Matrix& pos, const int i)
{
  int natom = pos.rows();
  int ndim = pos.cols();
  double dist;  // temporary pair distance
  double lap = 0.0;
  for (int j=0; j<natom; j++)
  {
    if (j==i) continue;
    dist = displacement(pos.row(i), pos.row(j)).norm();
    lap += lnf(dist)/std::pow(dist, 2);
  }
  return _a2*(_a2+2-ndim)*lap;
}

McMillanHe::Matrix
McMillanHe::diffuse(const Matrix& pos, const int nstep, const double tau)
{
  Matrix pos1(pos); // make a copy of initial positions
  _natt = 0;
  _nacc = 0;
  int natom = pos1.rows();
  int ndim = pos1.cols();
  double sig = sqrt(tau);
  // temporary variables
  double lna, prob;
  Vector move(ndim);
  for (int istep=0; istep<nstep; istep++)
  {
    for (int iatom=0; iatom<natom; iatom++)
    {
      // make move vector
      for (int idim=0; idim<ndim; idim++)
      {
        move(idim) = sig*_rng.randn();
      }
      // calculate acceptance ratio
      lna = diff_lnwf(pos1, move, iatom);
      prob = exp(2*lna);
      // attempt move
      _natt += 1;
      if (prob>_rng.rand())
      { // accept move
        pos1.row(iatom) += move;
        _nacc += 1;
      } else { // reject move
      }
    }
  }
  pos_in_box(pos1);
  return pos1;
}
