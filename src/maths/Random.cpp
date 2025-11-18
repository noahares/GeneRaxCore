#include "Random.hpp"

#include <cmath>

std::mt19937_64 Random::_rng;
std::uniform_int_distribution<int> Random::_uniint(0);
std::uniform_real_distribution<double> Random::_uniproba(0.0, 1.0);

void Random::setSeed(unsigned int seed) { _rng.seed(seed); }

int Random::getInt() { return _uniint(_rng); }

int Random::getInt(unsigned int min, unsigned int max) {
  std::uniform_int_distribution<int> distr(min, max);
  return distr(_rng);
}

bool Random::getBool() { return getInt() % 2; }

double Random::getProba() {
  // sometimes produces 1.0, though should not; see the link:
  // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  double proba = _uniproba(_rng);
  // convert 1.0 to the closest smaller double
  if (proba == 1.0)
    proba = std::nextafter(1.0, 0.0);
  return proba;
}
