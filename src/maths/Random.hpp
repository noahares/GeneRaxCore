#pragma once

#include <random>

class Random {
public:
  Random() = delete;
  static void setSeed(unsigned int seed);
  // return a random bool
  static bool getBool();
  // return a random unsigned int
  static unsigned int getUInt();
  // return a uniform random unsigned int from the [min,max] interval
  static unsigned int getUInt(unsigned int min, unsigned int max);
  // return a uniform random double from the [0,1) interval
  static double getProba();

private:
  static std::mt19937_64 _rng;
  static std::uniform_int_distribution<unsigned int> _uniuint;
  static std::uniform_real_distribution<double> _uniproba;
};
