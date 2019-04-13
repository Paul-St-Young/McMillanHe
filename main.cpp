#include <iostream>
#include <fstream>
#include "mmh.h"

#include "loadtxt.h"
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <time.h>
#include <stdio.h>

using namespace std;

typedef McMillanHe::Matrix Matrix;
typedef McMillanHe::Vector Vector;

void test_gl()
{
  int natom=2, ndim=3;
  double wf0, wf1, wf2, dr=1e-6;
  McMillanHe mmh = McMillanHe();
  Matrix pos(natom, ndim);
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
  Matrix pos(natom, ndim);
  Vector move(ndim);
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

Matrix get_fcc_pos(int natom, int ndim=3)
{
  Matrix pos(natom, ndim);
  char fname[100];
  sprintf(fname, "../pos%d.dat", natom);
  vector<vector<double>> pos0 = loadtxt(fname);
  for (int iatom=0; iatom<natom; iatom++)
    for (int idim=0; idim<ndim; idim++)
      pos(iatom, idim) = pos0[iatom][idim];
  return pos;
}

void vmc()
{
  int natom=32, ndim=3;
  McMillanHe mmh = McMillanHe();
  // read initial atomic positions
  Matrix pos = get_fcc_pos(natom, ndim);
  // VMC simulation
  //  initialize RNG rand() and randn()
  boost::random::mt19937 gen;
  boost::uniform_real<> udist(0.0, 1.0);  // uniform distribution
  boost::normal_distribution<> ndist(0.0, 1.0); // normal dist.
  boost::variate_generator<
    boost::mt19937&, boost::uniform_real<>
  > rand(gen, udist); // uniform [0, 1)
  boost::variate_generator<
    boost::mt19937&, boost::normal_distribution<>
  > randn(gen, ndist); // normal
  //  set simulation parameters
  mmh.set_a1(2.1);
  double lbox;
  if (natom == 32) lbox = 11.3303267;
  if (natom == 108) lbox = 16.99549;
  mmh.set_lbox(lbox);
  int nstep = 2560;
  int nacc = 0;
  double tau = 0.02;
  tau = 2.0;
  double lam = 1.0;
  double sig = sqrt(tau)*lam;
  double x2_forward, x2_backward;
  // temporary variables
  double lna, lnt, prob;
  bool use_drift=false;
  int iconf=10;
  Matrix pos1(natom, ndim), vel1(natom, ndim);
  Vector move(ndim), drift(ndim), drift1(ndim);
  Vector curpos(ndim), newpos(ndim);
  Vector x2_backvec(ndim);
  pos1 = pos; // make a copy of initial positions
  // save configurations
  ofstream fpos, fvel;
  fpos.open("all_pos.dat");
  fvel.open("all_vel.dat");
  fpos << "# natom=" << natom << " ndim=" << ndim << " lbox=" << lbox << endl;
  // time simulation
  clock_t begin, end;
  begin = clock();
  for (int istep=0; istep<nstep; istep++)
  {
    for (int iatom=0; iatom<natom; iatom++)
    {
      curpos = pos.row(iatom);
      // make move vector
      for (int idim=0; idim<ndim; idim++)
      {
        move(idim) = sig*randn();
      }
      if (use_drift)
      {
        x2_forward = move.squaredNorm();
        drift = tau*lam*2*mmh.grad_lnwf(pos1, iatom);
        vel1.row(iatom) = drift;
        move += drift;
      }
      // calculate acceptance ratio
      lna = mmh.diff_lnwf(pos1, move, iatom);
      newpos = curpos + move;
      pos1.row(iatom) = newpos;
      if (use_drift)
      {
        drift1 = tau*lam*2*mmh.grad_lnwf(pos1, iatom);
        x2_backvec = -drift1-move;
        for (int idim=0; idim<ndim; idim++)
        {
          x2_backvec(idim) -= round(x2_backvec(idim)/lbox)*lbox;
        }
        x2_backward = x2_backvec.squaredNorm();
        lnt = (x2_forward - x2_backward)/(2*pow(sig, 2));
        //lnt = (2*(newpos-curpos)+(drift1-drift)).squaredNorm();
      } else {
        lnt = 0.0;
      }
      prob = exp(2*lna+lnt);
      if (prob>rand())
      { // accept move
        nacc += 1;
      } else { // reject move
        pos1.row(iatom) = curpos;
      }
    }
    if (istep%iconf==0)
    {
      // bring configurations back into box
      for (int iatom=0; iatom<natom; iatom++)
        for (int idim=0; idim<ndim; idim++)
          pos1(iatom, idim) -= lbox*round(pos(iatom, idim)/lbox);
      fpos << pos1 << endl;
      if (use_drift) fvel << vel1 << endl;
    }
  }
  end = clock();
  cout << "elapsed time: " << (double)(end-begin)/CLOCKS_PER_SEC << " s." << endl;
  cout << ((double)nacc)/(nstep*natom) << endl;
  fpos.close();
  fvel.close();
}

int main(int argc, char** argv)
{
  //test_gl();
  //test_ratio();
  vmc();
  return 0;
}
