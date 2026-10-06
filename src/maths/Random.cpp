#include "Random.hpp"

#include <cassert>
#include <cmath>

std::mt19937_64 Random::_rng;
std::uniform_int_distribution<unsigned int> Random::_uniuint(0);
std::uniform_real_distribution<double> Random::_uniproba(0.0, 1.0);

void Random::setSeed(unsigned int seed) { _rng.seed(seed); }

bool Random::getBool() { return getUInt() % 2; }

unsigned int Random::getUInt() { return _uniuint(_rng); }

unsigned int Random::getUInt(unsigned int min, unsigned int max) {
  assert(min <= max);
  std::uniform_int_distribution<unsigned int> distr(min, max);
  return distr(_rng);
}

double Random::getProba() {
  // sometimes produces 1.0, though should not; see the link:
  // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  double proba = _uniproba(_rng);
  // convert 1.0 to the closest smaller double
  if (proba == 1.0) {
    proba = std::nextafter(proba, 0.0);
  }
  return proba;
}
