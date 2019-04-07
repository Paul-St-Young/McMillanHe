#include <mmh.h>

McMillanHe::Vector McMillanHe::displacement(Vector r1, Vector r2)
{
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
