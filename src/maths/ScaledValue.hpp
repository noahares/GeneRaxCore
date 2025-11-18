#pragma once

#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>

#define JS_SCALE_FACTOR                                                        \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0 /*  2**256 (exactly)  */
#define JS_SCALE_THRESHOLD (1.0 / JS_SCALE_FACTOR)

const int NULL_SCALER = INT_MAX / 2 - 1;

/**
 *  scale function for a general type
 *  (do nothing)
 */
template <class REAL> void scale(REAL &) {}

/**
 *  getLog function for a general type
 *  (apply std::log)
 */
template <class REAL> double getLog(const REAL &v) { return std::log(v); }

/**
 *  Class representing a double value with a high precision.
 *  It stores a double and a scaling integer to represent very
 *  small double values
 *
 *  When the value is null, the scaler is set to NULL_SCALER
 */
class ScaledValue {
public:
  /**
   *  Null value constructor
   */
  ScaledValue() : value(0.0), scaler(NULL_SCALER) {}

  /**
   *  Conversion constructor
   *  @param v value
   */
  explicit ScaledValue(double v) : value(v), scaler(0) {
    assert(value >= 0.0); // negative values not allowed
  }

  /**
   *  Conversion to a double
   */
  operator double() const {
    if (scaler == NULL_SCALER) {
      return 0.0;
    } else if (scaler == 0) {
      return value;
    } else { // the value is almost zero
      return 0.0;
    }
  }

  /**
   *  ScaledValue sum operator
   */
  inline ScaledValue operator+(const ScaledValue &v) const {
    if (v.scaler == scaler) {
      return ScaledValue(v.value + value, scaler);
    } else if (v.scaler < scaler) {
      return v;
    } else {
      return *this;
    }
  }

  /**
   *  ScaledValue sum operator
   */
  inline ScaledValue &operator+=(const ScaledValue &v) {
    if (v.scaler == scaler) {
      value += v.value;
    } else if (v.scaler < scaler) {
      value = v.value;
      scaler = v.scaler;
    }
    return *this;
  }

  /**
   *  ScaledValue minus operator
   */
  inline ScaledValue operator-(const ScaledValue &v) const {
    if (v.scaler == scaler) {
      if (value - v.value < 0.0) {
        if (fabs(value - v.value) < 0.0000000001) {
          return ScaledValue();
        }
        std::cerr.precision(17);
        std::cerr << *this << " - " << v << std::endl;
      }
      assert(value - v.value >= 0);
      auto res = ScaledValue(value - v.value, scaler);
      res.scale();
      return res;
    } else if (v.scaler < scaler) {
      std::cerr << *this << " - " << v << std::endl;
      assert(false); // negative values not allowed
      return v;
    } else {
      return *this;
    }
  }

  /**
   *  ScaledValue multiplication operator
   */
  inline ScaledValue operator*(const ScaledValue &v) const {
    auto res = ScaledValue(v.value * value, v.scaler + scaler);
    return res;
  }

  /**
   *  ScaledValue multiplication operator
   */
  inline ScaledValue &operator*=(const ScaledValue &v) {
    value *= v.value;
    scaler += v.scaler;
    return *this;
  }

  /**
   *  double multiplication operator
   */
  inline ScaledValue operator*(double v) const {
    auto res = ScaledValue(v * value, scaler);
    return res;
  }

  /**
   *  double multiplication operator
   */
  inline ScaledValue &operator*=(double v) {
    value *= v;
    return *this;
  }

  /**
   *  double division operator
   */
  inline ScaledValue operator/(double v) const {
    auto res = ScaledValue(value / v, scaler);
    return res;
  }

  /**
   *  double division operator
   */
  inline ScaledValue &operator/=(double v) {
    value /= v;
    return *this;
  }

  /**
   *  @return true if the value is 0
   */
  inline bool isNull() const { return value == 0.0; }

  /**
   *  Comparison with ScaledValue operators
   */
  inline bool operator<(const ScaledValue &v) const {
    if (isNull()) {
      return !v.isNull();
    }
    if (scaler != v.scaler) {
      return scaler > v.scaler;
    }
    return value < v.value;
  }

  inline bool operator>(const ScaledValue &v) const { return !(*this <= v); }

  inline bool operator==(const ScaledValue &v) const {
    if (isNull()) {
      return v.isNull();
    }
    return (scaler == v.scaler) && (value == v.value);
  }

  inline bool operator!=(const ScaledValue &v) const { return !(*this == v); }

  inline bool operator<=(const ScaledValue &v) const {
    if (isNull()) {
      return true;
    }
    if (scaler != v.scaler) {
      return scaler > v.scaler;
    }
    return value <= v.value;
  }

  inline bool operator>=(const ScaledValue &v) const { return !(*this < v); }

  /**
   *  std::ostream << operator
   */
  friend std::ostream &operator<<(std::ostream &os, const ScaledValue &v) {
    os << v.value << "s" << v.scaler;
    return os;
  }

  friend void scale<ScaledValue>(ScaledValue &v);
  friend double getLog<ScaledValue>(const ScaledValue &v);

private:
  /**
   *  General constructor
   *  @param v value
   *  @param s scaler
   */
  ScaledValue(double v, int s) : value(v), scaler(s) {}

  void checkNull() {
    if (value == 0.0) {
      scaler = NULL_SCALER;
    }
  }

  void scale() {
    if (value < JS_SCALE_THRESHOLD) {
      scaler += 1;
      value *= JS_SCALE_FACTOR;
      checkNull();
    }
  }

  double value;
  int scaler;
};

/**
 *  scale function for the ScaledValue type
 *  Should be applied every time when converting from a double
 *  or after a series of multiplication and/or division operations
 */
template <> inline void scale<ScaledValue>(ScaledValue &v) { v.scale(); }

/**
 *  getLog function for the ScaledValue type
 */
template <> inline double getLog<ScaledValue>(const ScaledValue &v) {
  if (v.scaler == NULL_SCALER) {
    return -std::numeric_limits<double>::infinity();
  }
  return std::log(v.value) + v.scaler * std::log(JS_SCALE_THRESHOLD);
}
