#ifndef RANDOM_NUMBER_GENERATOR_H
#define RANDOM_NUMBER_GENERATOR_H

#include <random>
#include <ctime>

class RandomNumberGenerator
{
public:
  RandomNumberGenerator(const int iseed)
  : _udist(0.0, 1.0), _ndist(0.0, 1.0), _iseed(iseed)
  {
    if (iseed==-1)
    {
      _iseed = static_cast<unsigned int>(std::time(0));
    }
    _gen().seed(_iseed);
  }
  ~RandomNumberGenerator() {};

  inline double rand()
  {
    return _udist(_gen());
  }

  inline double randn()
  {
    return _ndist(_gen());
  }
  int get_seed(){return _iseed;};
private:
  int _iseed;
  std::uniform_real_distribution<double> _udist;
  std::normal_distribution<double> _ndist;
  std::mt19937& _gen(){
    static thread_local std::mt19937 generator;
    return generator;
  }
};
#endif
