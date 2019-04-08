#include <iostream>
#include <fstream>
#include "mmh.h"

using namespace std;

void test_gl()
{
  int natom=2, ndim=3;
  double wf0, wf1, wf2, dr=1e-6;
  McMillanHe mmh = McMillanHe();
  McMillanHe::Matrix pos(natom, ndim);
  pos(0, 0) = -1.0;
  pos(1, 0) = 1.0;
  wf0 = mmh.lnwf(pos);
  int iatom = 0; // atom to move
  cout << mmh.grad_lnwf(pos, iatom)(0) << endl;
  // calculate FD gradient
  pos(iatom, 0) += dr;
  wf1 = mmh.lnwf(pos);
  cout << (wf1-wf0)/dr << endl;
  pos(iatom, 0) -= dr;
  // calculate FD laplacian
  cout << mmh.lap_lnwf(pos, iatom) << endl;
  double lap1=0.0, lap=0.0;
  for (int idim=0; idim<ndim; idim++)
  {
    pos(iatom, idim) -= dr;
    wf1 = mmh.lnwf(pos);
    pos(iatom, idim) += 2*dr;
    wf2 = mmh.lnwf(pos);
    pos(iatom, idim) -= dr;
    lap1 = (wf1+wf2-2*wf0)/pow(dr, 2);
    lap += lap1;
  }
  cout << lap << endl;
}

void test_ratio()
{
  int natom=2, ndim=3;
  double wf0, wf1, wf2, dr=0.2;
  McMillanHe mmh = McMillanHe();
  McMillanHe::Matrix pos(natom, ndim);
  McMillanHe::Vector move(ndim);
  move(0) = dr;
  pos(0, 0) = -1.0;
  pos(1, 0) = 1.0;
  wf0 = mmh.lnwf(pos);
  int iatom = 0; // atom to move
  cout << mmh.ratio(pos, move, iatom) << endl;
  // calculate FD ratio
  pos.row(iatom) += move;
  wf1 = mmh.lnwf(pos);
  cout << exp(wf1-wf0) << endl;
}

int main(int argc, char** argv)
{
  //test_gl();
  test_ratio();
  return 0;
}
