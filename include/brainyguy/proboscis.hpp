#pragma once
#ifndef BRAINYGUY_PROBOSCIS_H
#define BRAINYGUY_PROBOSCIS_H 1

/*
 *   Copyright 2022 Carlos Reyes
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <memory>

// -----------------------------------------------------------------------------
template<typename Scalar, size_t rows_compile_time, size_t cols_compile_time>
class Matrix
{

};

// -----------------------------------------------------------------------------
// https://en.cppreference.com/w/cpp/algorithm/inplace_merge
class SixSigma {
 public:
  MomentSketch()
  : _element_count(), _top_elements()
  {}

  void reset_values() {
    _element_count = 0;
    _top_elements = 1;
    _next_count = calc_next_count(_top_elements);
  }

  void add_element(const double sample) {
    ++_element_count;
    if (_element_count == 1) {
      _top[0] = sample;
      return;
    }

    if (_element_count == _next_count) {
      ++_top_elements;
      _next_count = calc_next_count(_top_elements);
    }

    if (sample > _top[_top_elements-1]) {

    }
  }


 private:
  consteval double six_sigma_quantile = 0.0000001973;
  consteval size_t max_top_elements = (1024/8)-1;
  consteval size_t max_elements = max_top_elements/six_sigma_quantile; // 643'689'812
  std::array<double, max_top_elements+1> _top;
  size_t _element_count;
  size_t _top_elements;
  size_t _next_count; // _element_count value for increasing _top_elements by one

  size_t calc_next_count(const size_t current_top_elements) {
    _next_count = lround((current_top_elements+1) / six_sigma_quantile);
  }
};

// -----------------------------------------------------------------------------
class SummaryStatistics {
 public:
  SummaryStatistics()
  : _count(0), _m1(0.0), _m2(0.0), _m3(0.0), _m4(0.0),
    _minimum(DBL_MAX), _maximum(DBL_MIN)
  {}

  void add_element(const double sample) {
    _count++;
    const double delta = sample - _m1;
    const double delta_n = delta / _count;
    const double delta_n2 = delta_n * delta_n;
    const double term = delta * delta_n * (_count - 1);
    _m1 += delta_n;
    _m4 += term * delta_n2 * (_count * _count - 3 * _count + 3) +
        6 * delta_n2 * _m2 - 4 * delta_n * _m3;
    _m3 += term * delta_n * (_count - 2) - 3 * delta_n * _m2;
    _m2 += term;
  }

  void reset_elements() {
    _count = 0;
    _m1 = _m2 = _m3 = _m4 = 0.0;
    _minimum = DBL_MAX;
    _maximum = DBL_MIN;
  }

  uint64_t samples() const {
    return _count;
  }

  // arithmetic mean
  double mean() const
  {
    return _m1;
  }

  // variance of unknown distribution
  // use 1.0 for gaussian distribution
  double variance() const
  {
    return _m2/(_count-1.5);
  }

  double standard_deviation() const
  {
    return std::sqrt(variance());
  }

  double skewness() const
  {
    return std::sqrt(static_cast<double>(n)) * _m3 / std::pow(_m2, 1.5);
  }

  double kurtosis() const
  {
    return static_cast<double>(n) * _m4 / (_m2*_m2) - 3.0;
  }

  double minimum() const
  {
    return _minimum;
  }

  double maximum() const
  {
    return _maximum;
  }

  double range() const
  {
    return _maximum - _minimum;
  }

  SummaryStatistics &operator+=(const SummaryStatistics &rhs) {
    *this = *this + rhs;
    return *this;
  }

  friend SummaryStatistics operator+(const SummaryStatistics& lhs, const SummaryStatistics &rhs) {
    SummaryStatistics sum;

    const double delta  = rhs._m1 - lhs._m1;
    const double delta2 = delta*delta;
    const double delta3 = delta*delta2;
    const double delta4 = delta2*delta2;

    sum._count = lhs._count + rhs._count;
    sum._m1 = (lhs._count*lhs._m1 + rhs._count*rhs._m1) / sum._count;
    sum._m2 = lhs._m2 + rhs._m2 + delta2 * lhs._count * rhs._count / sum._count;
    sum._m3 = lhs._m3 + rhs._m3 + delta3 * lhs._count * rhs._count * (lhs._count - rhs._count)/(sum._count*sum._count) +
              3.0*delta * (lhs._count*rhs._m2 - rhs._count*lhs._m2) / sum._count;
    sum._m4 = lhs._m4 + rhs._m4 + delta4*lhs._count*rhs._count * (lhs._count*lhs._count - lhs._count*rhs._count + rhs._count*rhs._count) /
              (sum._count*sum._count*sum._count) +
              6.0*delta2 * (lhs._count*lhs._count*rhs._m2 + rhs._count*rhs._count*lhs._m2)/(sum._count*sum._count) +
              4.0*delta*(lhs._count*rhs._m3 - rhs._count*lhs._m3) / sum._count;

    sum._minimum = std::min(lhs._minimum, rhs._minimum);
    sum._maximum = std::max(lhs._maximum, rhs._maximum);

    return sum;
  }

 private:
  uint64_t _count;
  double _m1, _m2, _m3, _m4;
  double _minimum, _maximum;
};

#endif // BRAINYGUY_PROBOSCIS_H
