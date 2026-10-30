#pragma once

#include <algorithm>
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
 *  getCheckedLog function for a general type: the log, or -inf when the value
 *  is zero or subnormal. Below the normal range a double has lost significant
 *  bits, so -inf tells the caller to recompute the value with ScaledValue.
 */
template <class REAL> double getCheckedLog(const REAL &v) {
  return std::isnormal(v) ? std::log(v)
                          : -std::numeric_limits<double>::infinity();
}

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
   *
   *  A band is 2^256 wide, so two values one scaler apart can have the same
   * magnitude (one at the bottom of its band, the other at the top of the
   * next): they are added at the smaller scaler. Two or more scalers apart the
   * smaller value is below 2^-256 of the larger one, far below double's
   * rounding, and is dropped exactly as before.
   */
  inline ScaledValue operator+(const ScaledValue &v) const {
    int scaler_diff = scaler - v.scaler;
    if (scaler_diff == 0) {
      return ScaledValue(v.value + value, scaler);
    } else if (scaler_diff == 1) {
      return ScaledValue(v.value + value * JS_SCALE_THRESHOLD, v.scaler);
    } else if (scaler_diff == -1) {
      return ScaledValue(value + v.value * JS_SCALE_THRESHOLD, scaler);
    } else if (scaler_diff >= 2) {
      return v;
    } else {
      return *this;
    }
  }

  /**
   *  ScaledValue sum operator (see operator+)
   */
  inline ScaledValue &operator+=(const ScaledValue &v) {
    int scaler_diff = scaler - v.scaler;
    if (scaler_diff == 0) {
      value += v.value;
    } else if (scaler_diff == 1) {
      value = v.value + value * JS_SCALE_THRESHOLD;
      scaler = v.scaler;
    } else if (scaler_diff == -1) {
      value += v.value * JS_SCALE_THRESHOLD;
    } else if (scaler_diff >= 2) {
      value = v.value;
      scaler = v.scaler;
    }
    return *this;
  }

  /**
   *  ScaledValue minus operator
   *
   *  Precondition: *this >= v. ScaledValues are probabilities, so a difference
   * below zero is a bug, except for rounding when the two are equal in exact
   * arithmetic (returned as 0). AleRax subtracts only a sum over all species
   *  minus a sum over a subset of the same stored terms (the transfer-forbidden
   * species), both accumulated with +=, which keeps the smallest scaler of its
   * terms; the minuend's scaler therefore never exceeds the subtrahend's there.
   *  Values one scaler apart are subtracted at the smaller scaler (see
   * operator+).
   */
  inline ScaledValue operator-(const ScaledValue &v) const {
    int scaler_diff = scaler - v.scaler;
    double a = value; // both at the smaller scaler s
    double b = v.value;
    int s = scaler;
    if (scaler_diff == 1) {
      a *= JS_SCALE_THRESHOLD;
      s = v.scaler;
    } else if (scaler_diff == -1) {
      b *= JS_SCALE_THRESHOLD;
    } else if (scaler_diff <= -2) {
      return *this; // v is below 2^-256 of *this
    } else if (scaler_diff >= 2) {
      std::cerr << *this << " - " << v << std::endl;
      assert(false); // negative values not allowed
      return ScaledValue();
    }
    const double r = a - b;
    if (r < 0.0) {
      // Rounding between two values equal in exact arithmetic, relative to the
      // operands (a mantissa can be as small as 2^-256, so an absolute
      // tolerance would accept a true negative).
      if (-r <= 0.0000000001 * std::max(a, b)) {
        return ScaledValue();
      }
      std::cerr.precision(17);
      std::cerr << *this << " - " << v << std::endl;
      assert(false); // negative values not allowed
      return ScaledValue();
    }
    auto res = ScaledValue(r, s);
    res.scale();
    res.checkNull();
    return res;
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
   *  Comparison with ScaledValue operators. As for the sum, values one scaler
   * apart can have the same magnitude and are compared at the smaller scaler;
   * two or more apart the scaler decides. Zero is below every positive value.
   */
  inline int compare(const ScaledValue &v) const {
    if (isNull() || v.isNull()) {
      return isNull() ? (v.isNull() ? 0 : -1) : 1;
    }
    int scaler_diff = scaler - v.scaler;
    double a = value;
    double b = v.value;
    if (scaler_diff == 1) {
      a *= JS_SCALE_THRESHOLD;
    } else if (scaler_diff == -1) {
      b *= JS_SCALE_THRESHOLD;
    } else if (scaler_diff != 0) {
      return scaler_diff > 0 ? -1 : 1;
    }
    return a < b ? -1 : (b < a ? 1 : 0);
  }

  inline bool operator<(const ScaledValue &v) const { return compare(v) < 0; }

  inline bool operator>(const ScaledValue &v) const { return compare(v) > 0; }

  inline bool operator==(const ScaledValue &v) const { return compare(v) == 0; }

  inline bool operator!=(const ScaledValue &v) const { return compare(v) != 0; }

  inline bool operator<=(const ScaledValue &v) const { return compare(v) <= 0; }

  inline bool operator>=(const ScaledValue &v) const { return compare(v) >= 0; }

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

/**
 *  getCheckedLog function for the ScaledValue type: a ScaledValue keeps its
 *  precision near zero, so this is getLog
 */
template <> inline double getCheckedLog<ScaledValue>(const ScaledValue &v) {
  return getLog(v);
}
