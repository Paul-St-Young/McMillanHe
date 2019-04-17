#ifndef MCMILLAN_HE_H
#define MCMILLAN_HE_H

#include <Eigen/Dense>
#include "rng.h"

class McMillanHe
{
public:
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::VectorXd Vector;
  McMillanHe(int seed=-1)
    : _a1(2.6), _a2(5.0), _lbox(11.3303267),
    _natt(0), _nacc(0), _rng(seed)
  {}
  void pos_in_box(Matrix& pos);
  Vector displacement(const Vector& r1, const Vector& r2);
  // wf value
  double lnf(const double r);
  double fval(const double r);
  double lnwf(const Matrix& pos);
  double wfval(const Matrix& pos);
  // wf ratio
  double diff_lnwf(const Matrix& pos, const Vector& move, const int i);
  double ratio(const Matrix& pos, const Vector& move, const int i);
  // wf derivative
  Vector grad_lnwf(const Matrix& pos, const int i);
  double lap_lnwf(const Matrix& pos, const int i);
  // QMC driver
  Matrix diffuse(const Matrix& pos,
    const int nstep, const double tau);
  double get_acc(){return ((double)_nacc)/_natt;}
  // getters
  double get_a1(){return _a1;}
  double get_a2(){return _a2;}
  double get_lbox(){return _lbox;}
  // setters
  void set_a1(double a1){_a1=a1;}
  void set_a2(double a2){_a2=a2;}
  void set_lbox(double lbox){_lbox=lbox;}
private:
  double _a1, _a2, _eps, _sig, _lbox;
  int _natt, _nacc;  // number of attempts and accepts
  RandomNumberGenerator _rng;
};
#endif
