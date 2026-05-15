#pragma once

#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>

using ScaledValueType = double;
constexpr unsigned SCALE_THRESHOLD =
    std::numeric_limits<ScaledValueType>::digits - 1;
constexpr ScaledValueType JS_SCALE_FACTOR =
    ScaledValueType(1ull << SCALE_THRESHOLD);
constexpr ScaledValueType JS_SCALE_THRESHOLD  = 1.0/JS_SCALE_FACTOR;
constexpr ScaledValueType LOG_2 =
    static_cast<ScaledValueType>(0.693147180559945309417232121458176568);

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
  explicit ScaledValue(ScaledValueType v) : value(v), scaler(0) {
    assert(value >= 0.0); // negative values not allowed
    scale();
  }

  explicit ScaledValue(double v, int s):value(v), scaler(s)  {
    assert(value >= 0.0); // negative values not allowed
    scale();
  }

  /**
   *  Conversion to a double
   */
  operator ScaledValueType() const {
    if (scaler == NULL_SCALER) {
      return 0.0;
    }
    return std::ldexp(value, -scaler);
  }

  /**
   *  ScaledValue sum operator
   */
  inline ScaledValue operator+(const ScaledValue &v) const {
    ScaledValue result{*this};
    result += v;
    return result;
  }

  /**
   *  ScaledValue sum operator
   */
  inline ScaledValue &operator+=(const ScaledValue &v) {
    if (v.isNull()) {
      return *this;
    }
    if (isNull()) {
      value = v.value;
      scaler = v.scaler;
      return *this;
    }

    const int common_scaler = std::min(scaler, v.scaler);
    const double lhs = std::ldexp(value, common_scaler - scaler);
    const double rhs = std::ldexp(v.value, common_scaler - v.scaler);
    value = lhs + rhs;
    scaler = common_scaler;
    scale();
    return *this;
  }

  /**
   *  ScaledValue minus operator
   */
  inline ScaledValue operator-(const ScaledValue &v) const {
    ScaledValue result{*this};
    result += ScaledValue{-v.value, v.scaler};
    return result;
  }

  /**
   *  ScaledValue multiplication operator
   */
  inline ScaledValue operator*(const ScaledValue &v) const {
    ScaledValue res{*this};
    res *= v;
    return res;
  }

  /**
   *  ScaledValue multiplication operator
   */
  inline ScaledValue &operator*=(const ScaledValue &v) {
    if (isNull() || v.isNull()) {
      setNull();
      return *this;
    }
    value *= v.value;
    scaler += v.scaler;
    scale();
    return *this;
  }

  /**
   *  double multiplication operator
   */
  inline ScaledValue operator*(ScaledValueType v) const {
    return *this * ScaledValue{v};
  }

  /**
   *  double multiplication operator
   */
  inline ScaledValue &operator*=(ScaledValueType v) {
    return *this *= ScaledValue{v};
  }

  /**
   *  double division operator
   */
  inline ScaledValue operator/(ScaledValueType v) const {
    ScaledValue res{v};
    res /= v;
    return res;
  }

  /**
   *  double division operator
   */
  inline ScaledValue &operator/=(ScaledValueType v) {
    assert(v != 0.0);
    if (isNull()) {
      *this = ScaledValue{};
      return *this;
    }
    value /= v;
    scale();
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
    if (isNull() && v.isNull()) {
      return false;
    }
    if (isNull()) {
      return v.value > 0.0;
    }
    if (v.isNull()) {
      return value < 0.0;
    }

    assert(value >= 0.0);
    assert(v.value >= 0.0);
    if (scaler != v.scaler) {
      return scaler > v.scaler;
    }
    return value < v.value;
  }

  inline bool operator>(const ScaledValue &v) const { return !(*this <= v); }

  inline bool operator==(const ScaledValue &v) const {
    if (isNull() || v.isNull()) {
      return isNull() && v.isNull();
    }
    return scaler == v.scaler &&
      (std::fabs(v.value - value) <=
       std::numeric_limits<double>::epsilon() *
       std::max({1.0, std::fabs(value), std::fabs(v.value)}));
  }

  inline bool operator!=(const ScaledValue &v) const { return !(*this == v); }

  inline bool operator<=(const ScaledValue &v) const {
    return *this < v || *this == v;
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

  void checkNull() {
    if (value == 0.0) {
      scaler = NULL_SCALER;
    }
  }

  void setNull() {
    value = 0.0;
    scaler = NULL_SCALER;
  }

  void scale() {
    if (value == 0.0) {
      setNull();
      return;
    }

    int exponent = 0;
    value = std::frexp(value, &exponent);
    scaler -= exponent;
    checkNull();
  }

  ScaledValueType value;
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
    if (v.isNull()) {
      return -std::numeric_limits<double>::infinity();
    }
    if (v.value < 0.0) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return std::log(v.value) - static_cast<double>(v.scaler) * LOG_2;
}
