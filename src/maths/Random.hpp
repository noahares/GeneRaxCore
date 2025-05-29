#pragma once

#include <random>

class Random {
public:
  Random() = delete;
  static void setSeed(unsigned int seed);
  // return a non-negative random int
  static int getInt();
  // return a uniform random int from the [min,max] interval
  static int getInt(unsigned int min, unsigned int max);
  // return a random bool
  static bool getBool();
  // return a uniform random double from the [0,1) interval
  static double getProba();

private:
  static std::mt19937_64 _rng;
  static std::uniform_int_distribution<int> _uniint;
  static std::uniform_real_distribution<double> _uniproba;
};
