#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace boost;

class RandomNumberGenerator
{
public:
  RandomNumberGenerator(const int iseed=-1);
  ~RandomNumberGenerator();
  double rand();
  double randn();
private:
  int _iseed;
  variate_generator<mt19937&, uniform_real<>>* _grand;
  variate_generator<mt19937&, normal_distribution<>>* _grandn;
  random::mt19937 _gen;
  uniform_real<> _udist;
  normal_distribution<> _ndist;
};
