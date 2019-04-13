#include <ctime>
#include "rng.h"

RandomNumberGenerator::RandomNumberGenerator(const int iseed)
  : _udist(0.0, 1.0), _ndist(0.0, 1.0), _iseed(iseed)
{
  if (iseed==-1)
  {
    _iseed = static_cast<unsigned int>(std::time(0));
    _gen.seed(_iseed);
  }
  variate_generator<mt19937&, uniform_real<>>
    rand(_gen, _udist);
  _grand = new variate_generator<mt19937&, uniform_real<>>
    (_gen, _udist);
  _grandn = new variate_generator<mt19937&, normal_distribution<>>
    (_gen, _ndist);
}

RandomNumberGenerator::~RandomNumberGenerator()
{
  delete _grand;
  delete _grandn;
}

double RandomNumberGenerator::rand()
{
  return (*_grand)();
}

double RandomNumberGenerator::randn()
{
  return (*_grandn)();
}
