#include <Eigen/Dense>
#include <iostream>

class McMillanHe
{
public:
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::VectorXd Vector;
  McMillanHe()
    : _a1(2.6), _a2(5.0), _lbox(11.3303267)
  {}
  Vector displacement(Vector r1, Vector r2);
  // wf value
  double lnf(double r);
  double fval(double r);
  double lnwf(Matrix pos);
  double wfval(Matrix pos);
  // wf ratio
  double diff_lnwf(Matrix pos, Vector move, int i);
  double ratio(Matrix pos, Vector move, int i);
  // wf derivative
  Vector grad_lnwf(Matrix pos, int i);
  double lap_lnwf(Matrix pos, int i);
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
};
