#pragma once

#include "types.h"

#include <complex>
#include <vector>

namespace dpcpp {

bool is_power_of_two(int value);
void fft(std::vector<std::complex<Real>>& values, bool inverse);
void dct_orthonormal(std::vector<Real>& values);
void idct_orthonormal(std::vector<Real>& values);
void dct2_orthonormal(std::vector<Real>& values, int nx, int ny);
void idct2_orthonormal(std::vector<Real>& values, int nx, int ny);

}  // namespace dpcpp

