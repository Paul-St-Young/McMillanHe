#include <mmh.h>

McMillanHe::Vector McMillanHe::displacement(Vector r1, Vector r2)
{
  int ndim = r1.size();
  int nint;
  Vector dr(ndim);
  for (int idim=0; idim<ndim; idim++)
  {
    dr(idim) = r1(idim)-r2(idim);
    nint = std::round(dr(idim)/_lbox);
    dr(idim) -= nint*_lbox;
  }
  return r1-r2;
}

double McMillanHe::lnf(double r)
{
  return -std::pow(_a1/r, _a2);
}

double McMillanHe::pairf(double r)
{
  return std::exp(lnf(r));
}

double McMillanHe::lnwf(Matrix pos)
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

double McMillanHe::wfval(Matrix pos)
{
  return std::exp(lnwf(pos));
}

double McMillanHe::diff_lnwf(Matrix pos, Vector move, int i)
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

double McMillanHe::ratio(Matrix pos, Vector move, int i)
{
  return std::exp(diff_lnwf(pos, move, i));
}

McMillanHe::Vector McMillanHe::grad_lnwf(Matrix pos, int i)
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

double McMillanHe::lap_lnwf(Matrix pos, int i)
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
