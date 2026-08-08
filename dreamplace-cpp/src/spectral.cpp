#include "spectral.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dpcpp {
namespace {
constexpr Real kPi = 3.141592653589793238462643383279502884;
}

bool is_power_of_two(int value) {
    return value > 0 && (value & (value - 1)) == 0;
}

void fft(std::vector<std::complex<Real>>& values, bool inverse) {
    const int n = static_cast<int>(values.size());
    if (!is_power_of_two(n)) throw std::runtime_error("FFT length must be a power of two");
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(values[i], values[j]);
    }
    for (int length = 2; length <= n; length <<= 1) {
        const Real angle = (inverse ? 2.0 : -2.0) * kPi / length;
        const std::complex<Real> root(std::cos(angle), std::sin(angle));
        for (int begin = 0; begin < n; begin += length) {
            std::complex<Real> w(1.0, 0.0);
            for (int j = 0; j < length / 2; ++j) {
                const auto even = values[begin + j];
                const auto odd = values[begin + j + length / 2] * w;
                values[begin + j] = even + odd;
                values[begin + j + length / 2] = even - odd;
                w *= root;
            }
        }
    }
    if (inverse) {
        for (auto& value : values) value /= static_cast<Real>(n);
    }
}

void dct_orthonormal(std::vector<Real>& values) {
    const int n = static_cast<int>(values.size());
    std::vector<std::complex<Real>> extension(2 * n);
    for (int i = 0; i < n; ++i) {
        extension[i] = values[i];
        extension[2 * n - 1 - i] = values[i];
    }
    fft(extension, false);
    for (int k = 0; k < n; ++k) {
        const Real angle = -kPi * k / (2.0 * n);
        const std::complex<Real> phase(std::cos(angle), std::sin(angle));
        const Real sum = 0.5 * std::real(extension[k] * phase);
        values[k] = sum * (k == 0 ? std::sqrt(1.0 / n) : std::sqrt(2.0 / n));
    }
}

void idct_orthonormal(std::vector<Real>& values) {
    const int n = static_cast<int>(values.size());
    std::vector<std::complex<Real>> spectrum(2 * n, {0.0, 0.0});
    for (int k = 0; k < n; ++k) {
        const Real sum = values[k] /
            (k == 0 ? std::sqrt(1.0 / n) : std::sqrt(2.0 / n));
        const Real angle = kPi * k / (2.0 * n);
        spectrum[k] = 2.0 * sum * std::complex<Real>(std::cos(angle), std::sin(angle));
    }
    spectrum[n] = {0.0, 0.0};
    for (int k = 1; k < n; ++k) spectrum[2 * n - k] = std::conj(spectrum[k]);
    fft(spectrum, true);
    for (int i = 0; i < n; ++i) values[i] = spectrum[i].real();
}

void dct2_orthonormal(std::vector<Real>& values, int nx, int ny) {
    if (static_cast<int>(values.size()) != nx * ny) throw std::runtime_error("DCT2 shape mismatch");
    #pragma omp parallel
    {
        std::vector<Real> line(std::max(nx, ny));
        #pragma omp for schedule(static)
        for (int y = 0; y < ny; ++y) {
            line.resize(nx);
            for (int x = 0; x < nx; ++x) line[x] = values[y * nx + x];
            dct_orthonormal(line);
            for (int x = 0; x < nx; ++x) values[y * nx + x] = line[x];
        }
        #pragma omp for schedule(static)
        for (int x = 0; x < nx; ++x) {
            line.resize(ny);
            for (int y = 0; y < ny; ++y) line[y] = values[y * nx + x];
            dct_orthonormal(line);
            for (int y = 0; y < ny; ++y) values[y * nx + x] = line[y];
        }
    }
}

void idct2_orthonormal(std::vector<Real>& values, int nx, int ny) {
    if (static_cast<int>(values.size()) != nx * ny) throw std::runtime_error("IDCT2 shape mismatch");
    #pragma omp parallel
    {
        std::vector<Real> line(std::max(nx, ny));
        #pragma omp for schedule(static)
        for (int y = 0; y < ny; ++y) {
            line.resize(nx);
            for (int x = 0; x < nx; ++x) line[x] = values[y * nx + x];
            idct_orthonormal(line);
            for (int x = 0; x < nx; ++x) values[y * nx + x] = line[x];
        }
        #pragma omp for schedule(static)
        for (int x = 0; x < nx; ++x) {
            line.resize(ny);
            for (int y = 0; y < ny; ++y) line[y] = values[y * nx + x];
            idct_orthonormal(line);
            for (int y = 0; y < ny; ++y) values[y * nx + x] = line[y];
        }
    }
}

namespace {

// Synthesize the sine basis sin((sample + 1/2) k pi / n), k=1..n-1,
// through the existing orthonormal DCT-III implementation.  The identity
// sin((i+1/2)k*pi/n)=(-1)^i cos((i+1/2)(n-k)*pi/n) avoids a second FFT kernel.
void inverse_sine_orthonormal(std::vector<Real>& values) {
    const int n = static_cast<int>(values.size());
    std::vector<Real> reversed(n, 0.0);
    for (int k = 1; k < n; ++k) reversed[n - k] = values[k];
    idct_orthonormal(reversed);
    for (int i = 0; i < n; ++i) values[i] = (i & 1) ? -reversed[i] : reversed[i];
}

}  // namespace

void inverse_mixed_sine_cosine2(std::vector<Real>& values, int nx, int ny,
                                bool sine_x) {
    if (static_cast<int>(values.size()) != nx * ny)
        throw std::runtime_error("mixed inverse transform shape mismatch");
    #pragma omp parallel
    {
        std::vector<Real> line(std::max(nx, ny));
        #pragma omp for schedule(static)
        for (int y = 0; y < ny; ++y) {
            line.resize(nx);
            for (int x = 0; x < nx; ++x) line[x] = values[y * nx + x];
            if (sine_x) inverse_sine_orthonormal(line);
            else idct_orthonormal(line);
            for (int x = 0; x < nx; ++x) values[y * nx + x] = line[x];
        }
        #pragma omp for schedule(static)
        for (int x = 0; x < nx; ++x) {
            line.resize(ny);
            for (int y = 0; y < ny; ++y) line[y] = values[y * nx + x];
            if (sine_x) idct_orthonormal(line);
            else inverse_sine_orthonormal(line);
            for (int y = 0; y < ny; ++y) values[y * nx + x] = line[y];
        }
    }
}

}  // namespace dpcpp
