#include <Eigen/Dense>
#include <iostream>

class McMillanHe
{
public:
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::VectorXd Vector;
  McMillanHe()
    : _a1(2.6), _a2(5.0)
  {}
  double get_a1(){return _a1;}
  double get_a2(){return _a2;}
  Vector displacement(Vector r1, Vector r2);
  double lnf(double r);
  double pairf(double r);
  double lnwf(Matrix pos);
  double wfval(Matrix pos);
  Vector grad_lnwf(Matrix pos, int i);
  double lap_lnwf(Matrix pos, int i);
private:
  double _a1, _a2;
};
