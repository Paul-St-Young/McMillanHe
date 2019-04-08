#include <iostream>
#include <fstream>
#include "mmh.h"

#include "loadtxt.h"
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

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

void vmc()
{
  int natom=32, ndim=3;
  McMillanHe mmh = McMillanHe();
  McMillanHe::Matrix pos(natom, ndim);
  // read initial atomic positions
  vector<vector<double>> pos0 = loadtxt("../pos.dat");
  for (int iatom=0; iatom<natom; iatom++)
    for (int idim=0; idim<ndim; idim++)
      pos(iatom, idim) = pos0[iatom][idim];
  //cout << pos << endl;
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
  double lbox = mmh.get_lbox();
  int nstep = 256;
  int nacc = 0;
  double tau = 0.01;
  double lam = 0.5;
  // temporary variables
  double ratio_wf2, lna;
  McMillanHe::Matrix pos1(natom, ndim);
  McMillanHe::Vector move(ndim), drift(ndim), drift1(ndim);
  pos1 = pos; // make a copy of initial positions
  // save configurations
  ofstream fpos, fvel;
  fpos.open("all_pos.dat");
  //fvel.open("all_vel.dat");
  fpos << "# natom=" << natom << " ndim=" << ndim << " lbox=" << lbox << endl;
  for (int istep=0; istep<nstep; istep++)
  {
    for (int iatom=0; iatom<natom; iatom++)
    {
      // make move vector
      //drift = tau*lam*2*mmh.grad_lnwf(pos1);
      for (int idim=0; idim<ndim; idim++)
      {
        move(idim) = sqrt(tau)*lam*randn();// + drift;
      }
      ratio_wf2 = pow(mmh.ratio(pos1, move, iatom), 2);
      if (ratio_wf2>rand())
      { // accept move
        pos1.row(iatom) += move;
        nacc += 1;
      } else { // reject move
      }
    }
    fpos << pos1 << endl;
  }
  cout << ((double)nacc)/(nstep*natom) << endl;
  fpos.close();
}

int main(int argc, char** argv)
{
  //test_gl();
  //test_ratio();
  vmc();
  return 0;
}
