#include "spectral.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <stdexcept>

namespace dpcpp {
namespace {
constexpr Real kPi = 3.141592653589793238462643383279502884;

struct FftPlan {
    int size = 0;
    std::vector<int> reversed;
    std::vector<std::vector<std::complex<Real>>> forward_twiddles;
    std::vector<std::vector<std::complex<Real>>> inverse_twiddles;
};

struct DctPlan {
    int size = 0;
    std::vector<std::complex<Real>> forward_phase;
    std::vector<std::complex<Real>> inverse_phase;
    std::vector<Real> forward_scale;
    std::vector<Real> inverse_scale;
    std::vector<Real> lee_dct_cos;
    std::vector<Real> lee_idct_cos;
    std::vector<Real> orthonormal_scale;
    std::vector<Real> orthonormal_inverse_input_scale;
};

struct RealTransformWorkspace {
    std::vector<Real> stage;
    std::vector<Real> buffer;
    std::vector<Real> transposed;
};

FftPlan make_fft_plan(int size) {
    if (!is_power_of_two(size))
        throw std::runtime_error("FFT length must be a power of two");
    FftPlan plan;
    plan.size = size;
    plan.reversed.resize(size);
    for (int i = 1, j = 0; i < size; ++i) {
        int bit = size >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        plan.reversed[i] = j;
    }
    for (int length = 2; length <= size; length <<= 1) {
        const Real angle = 2.0 * kPi / length;
        const std::complex<Real> forward_root(std::cos(-angle), std::sin(-angle));
        const std::complex<Real> inverse_root(std::cos(angle), std::sin(angle));
        std::vector<std::complex<Real>> forward(length / 2);
        std::vector<std::complex<Real>> inverse(length / 2);
        std::complex<Real> wf(1.0, 0.0);
        std::complex<Real> wi(1.0, 0.0);
        for (int j = 0; j < length / 2; ++j) {
            forward[j] = wf;
            inverse[j] = wi;
            wf *= forward_root;
            wi *= inverse_root;
        }
        plan.forward_twiddles.push_back(std::move(forward));
        plan.inverse_twiddles.push_back(std::move(inverse));
    }
    return plan;
}

const FftPlan& fft_plan(int size) {
    static std::mutex mutex;
    static std::map<int, FftPlan> plans;
    std::lock_guard<std::mutex> lock(mutex);
    auto found = plans.find(size);
    if (found == plans.end())
        found = plans.emplace(size, make_fft_plan(size)).first;
    return found->second;
}

DctPlan make_dct_plan(int size) {
    if (!is_power_of_two(size))
        throw std::runtime_error("DCT length must be a power of two");
    DctPlan plan;
    plan.size = size;
    plan.forward_phase.resize(size);
    plan.inverse_phase.resize(size);
    plan.forward_scale.resize(size);
    plan.inverse_scale.resize(size);
    plan.lee_dct_cos.resize(std::max(0, size - 1));
    plan.lee_idct_cos.resize(std::max(0, size - 1));
    plan.orthonormal_scale.resize(size);
    plan.orthonormal_inverse_input_scale.resize(size);
    for (int k = 0; k < size; ++k) {
        const Real angle = kPi * k / (2.0 * size);
        plan.forward_phase[k] = {std::cos(-angle), std::sin(-angle)};
        plan.inverse_phase[k] = {std::cos(angle), std::sin(angle)};
        const Real normalization = k == 0
            ? std::sqrt(1.0 / size) : std::sqrt(2.0 / size);
        plan.forward_scale[k] = 0.5 * normalization;
        plan.inverse_scale[k] = 2.0 / normalization;
        plan.orthonormal_scale[k] = normalization;
        plan.orthonormal_inverse_input_scale[k] =
            k == 0 ? 2.0 * normalization : normalization;
    }
    int offset = 0;
    for (int half = size / 2; half > 0; half /= 2) {
        const Real step = 0.5 * kPi / half;
        Real phase = 0.5 * step;
        for (int i = 0; i < half; ++i, phase += step)
            plan.lee_dct_cos[offset + i] = 0.5 / std::cos(phase);
        offset += half;
    }
    offset = 0;
    for (int half = 1; half < size; half *= 2) {
        const Real step = 0.5 * kPi / half;
        Real phase = 0.5 * step;
        for (int i = 0; i < half; ++i, phase += step)
            plan.lee_idct_cos[offset + i] = 0.5 / std::cos(phase);
        offset += half;
    }
    return plan;
}

const DctPlan& dct_plan(int size) {
    static std::mutex mutex;
    static std::map<int, DctPlan> plans;
    std::lock_guard<std::mutex> lock(mutex);
    auto found = plans.find(size);
    if (found == plans.end())
        found = plans.emplace(size, make_dct_plan(size)).first;
    return found->second;
}

void fft_transform(std::vector<std::complex<Real>>& values, bool inverse,
                   const FftPlan& plan) {
    const int n = static_cast<int>(values.size());
    if (n != plan.size) throw std::runtime_error("FFT plan size mismatch");
    for (int i = 1; i < n; ++i) {
        const int j = plan.reversed[i];
        if (i < j) std::swap(values[i], values[j]);
    }
    const auto& stages = inverse ? plan.inverse_twiddles : plan.forward_twiddles;
    int length = 2;
    for (const auto& twiddles : stages) {
        for (int begin = 0; begin < n; begin += length) {
            for (int j = 0; j < length / 2; ++j) {
                const auto even = values[begin + j];
                const auto odd = values[begin + j + length / 2] * twiddles[j];
                values[begin + j] = even + odd;
                values[begin + j + length / 2] = even - odd;
            }
        }
        length <<= 1;
    }
    if (inverse) {
        const Real scale = 1.0 / static_cast<Real>(n);
        for (auto& value : values) value *= scale;
    }
}

void dct_transform(std::vector<Real>& values, const DctPlan& plan,
                   const FftPlan& fft,
                   std::vector<std::complex<Real>>& work) {
    const int n = static_cast<int>(values.size());
    work.resize(2 * n);
    for (int i = 0; i < n; ++i) {
        work[i] = values[i];
        work[2 * n - 1 - i] = values[i];
    }
    fft_transform(work, false, fft);
    for (int k = 0; k < n; ++k) {
        values[k] = std::real(work[k] * plan.forward_phase[k]) *
                    plan.forward_scale[k];
    }
}

void idct_transform(std::vector<Real>& values, const DctPlan& plan,
                    const FftPlan& fft,
                    std::vector<std::complex<Real>>& work) {
    const int n = static_cast<int>(values.size());
    work.resize(2 * n);
    for (int k = 0; k < n; ++k) {
        work[k] = values[k] * plan.inverse_scale[k] * plan.inverse_phase[k];
    }
    work[n] = {0.0, 0.0};
    for (int k = 1; k < n; ++k) work[2 * n - k] = std::conj(work[k]);
    fft_transform(work, true, fft);
    for (int i = 0; i < n; ++i) values[i] = work[i].real();
}

void inverse_sine_transform(std::vector<Real>& values, const DctPlan& plan,
                            const FftPlan& fft,
                            std::vector<Real>& reversed,
                            std::vector<std::complex<Real>>& work) {
    const int n = static_cast<int>(values.size());
    reversed.resize(n);
    std::fill(reversed.begin(), reversed.end(), 0.0);
    for (int k = 1; k < n; ++k) reversed[n - k] = values[k];
    idct_transform(reversed, plan, fft, work);
    for (int i = 0; i < n; ++i)
        values[i] = (i & 1) ? -reversed[i] : reversed[i];
}

void lee_dct_orthonormal(const Real* input, Real* output, Real* buffer,
                         const DctPlan& plan) {
    const int n = plan.size;
    Real* current = output;
    Real* next = buffer;
    std::copy(input, input + n, current);
    int length = n;
    int half = length / 2;
    int cosine_offset = 0;
    while (half) {
        int offset = 0;
        const int steps = n / length;
        for (int k = 0; k < steps; ++k) {
            for (int i = 0; i < half; ++i) {
                next[offset + i] = current[offset + i] +
                                   current[offset + length - i - 1];
                next[offset + half + i] =
                    (current[offset + i] - current[offset + length - i - 1]) *
                    plan.lee_dct_cos[cosine_offset + i];
            }
            offset += length;
        }
        std::swap(current, next);
        cosine_offset += half;
        length = half;
        half /= 2;
    }
    length = 4;
    half = 2;
    while (half < n) {
        int offset = 0;
        const int steps = n / length;
        for (int k = 0; k < steps; ++k) {
            for (int i = 0; i < half - 1; ++i) {
                next[offset + i * 2] = current[offset + i];
                next[offset + i * 2 + 1] =
                    current[offset + half + i] +
                    current[offset + half + i + 1];
            }
            next[offset + length - 2] = current[offset + half - 1];
            next[offset + length - 1] = current[offset + length - 1];
            offset += length;
        }
        std::swap(current, next);
        half = length;
        length *= 2;
    }
    if (current != output) std::copy(current, current + n, output);
    for (int i = 0; i < n; ++i) output[i] *= plan.orthonormal_scale[i];
}

void lee_idct_orthonormal(const Real* input, Real* output, Real* buffer,
                          const DctPlan& plan) {
    const int n = plan.size;
    Real* current = output;
    Real* next = buffer;
    for (int i = 0; i < n; ++i)
        current[i] = input[i] * plan.orthonormal_inverse_input_scale[i];
    current[0] *= 0.5;
    int length = n;
    int half = length / 2;
    while (half) {
        int offset = 0;
        const int steps = n / length;
        for (int k = 0; k < steps; ++k) {
            next[offset] = current[offset];
            next[offset + half] = current[offset + 1];
            for (int i = 1; i < half; ++i) {
                next[offset + i] = current[offset + i * 2];
                next[offset + half + i] =
                    current[offset + i * 2 - 1] +
                    current[offset + i * 2 + 1];
            }
            offset += length;
        }
        std::swap(current, next);
        length = half;
        half /= 2;
    }
    length = 2;
    half = 1;
    int cosine_offset = 0;
    while (half < n) {
        int offset = 0;
        const int steps = n / length;
        for (int k = 0; k < steps; ++k) {
            for (int i = 0; i < half; ++i) {
                const Real g = current[offset + i];
                const Real h = current[offset + half + i] *
                               plan.lee_idct_cos[cosine_offset + i];
                next[offset + i] = g + h;
                next[offset + length - 1 - i] = g - h;
            }
            offset += length;
        }
        std::swap(current, next);
        cosine_offset += half;
        half = length;
        length *= 2;
    }
    if (current != output) std::copy(current, current + n, output);
}

void transpose_blocked(const Real* input, Real* output, int rows, int columns) {
    constexpr int block = 32;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int row_block = 0; row_block < rows; row_block += block) {
        for (int column_block = 0; column_block < columns;
             column_block += block) {
            const int row_end = std::min(rows, row_block + block);
            const int column_end = std::min(columns, column_block + block);
            for (int row = row_block; row < row_end; ++row) {
                for (int column = column_block; column < column_end; ++column)
                    output[column * rows + row] = input[row * columns + column];
            }
        }
    }
}

void lee_transform2(std::vector<Real>& values, int nx, int ny, bool inverse) {
    const DctPlan& plan_x = dct_plan(nx);
    const DctPlan& plan_y = dct_plan(ny);
    static RealTransformWorkspace workspace;
    const std::size_t count = static_cast<std::size_t>(nx) * ny;
    workspace.stage.resize(count);
    workspace.buffer.resize(count);
    workspace.transposed.resize(count);
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < ny; ++y) {
        if (inverse) {
            lee_idct_orthonormal(values.data() + y * nx,
                                 workspace.stage.data() + y * nx,
                                 workspace.buffer.data() + y * nx, plan_x);
        } else {
            lee_dct_orthonormal(values.data() + y * nx,
                                workspace.stage.data() + y * nx,
                                workspace.buffer.data() + y * nx, plan_x);
        }
    }
    transpose_blocked(workspace.stage.data(), workspace.transposed.data(), ny, nx);
    #pragma omp parallel for schedule(static)
    for (int x = 0; x < nx; ++x) {
        if (inverse) {
            lee_idct_orthonormal(workspace.transposed.data() + x * ny,
                                 workspace.stage.data() + x * ny,
                                 workspace.buffer.data() + x * ny, plan_y);
        } else {
            lee_dct_orthonormal(workspace.transposed.data() + x * ny,
                                workspace.stage.data() + x * ny,
                                workspace.buffer.data() + x * ny, plan_y);
        }
    }
    transpose_blocked(workspace.stage.data(), values.data(), nx, ny);
}

}  // namespace

bool is_power_of_two(int value) {
    return value > 0 && (value & (value - 1)) == 0;
}

void fft(std::vector<std::complex<Real>>& values, bool inverse) {
    fft_transform(values, inverse, fft_plan(static_cast<int>(values.size())));
}

void dct_orthonormal(std::vector<Real>& values) {
    std::vector<std::complex<Real>> work;
    const int n = static_cast<int>(values.size());
    dct_transform(values, dct_plan(n), fft_plan(2 * n), work);
}

void idct_orthonormal(std::vector<Real>& values) {
    std::vector<std::complex<Real>> work;
    const int n = static_cast<int>(values.size());
    idct_transform(values, dct_plan(n), fft_plan(2 * n), work);
}

void dct2_orthonormal(std::vector<Real>& values, int nx, int ny) {
    if (static_cast<int>(values.size()) != nx * ny)
        throw std::runtime_error("DCT2 shape mismatch");
    lee_transform2(values, nx, ny, false);
}

void idct2_orthonormal(std::vector<Real>& values, int nx, int ny) {
    if (static_cast<int>(values.size()) != nx * ny)
        throw std::runtime_error("IDCT2 shape mismatch");
    lee_transform2(values, nx, ny, true);
}

void inverse_mixed_sine_cosine2(std::vector<Real>& values, int nx, int ny,
                                bool sine_x) {
    if (static_cast<int>(values.size()) != nx * ny)
        throw std::runtime_error("mixed inverse transform shape mismatch");
    const DctPlan& plan_x = dct_plan(nx);
    const DctPlan& plan_y = dct_plan(ny);
    const FftPlan& fft_x = fft_plan(2 * nx);
    const FftPlan& fft_y = fft_plan(2 * ny);
    #pragma omp parallel
    {
        std::vector<Real> line(std::max(nx, ny));
        std::vector<Real> reversed(std::max(nx, ny));
        std::vector<std::complex<Real>> work(2 * std::max(nx, ny));
        #pragma omp for schedule(static)
        for (int y = 0; y < ny; ++y) {
            line.resize(nx);
            for (int x = 0; x < nx; ++x) line[x] = values[y * nx + x];
            if (sine_x)
                inverse_sine_transform(line, plan_x, fft_x, reversed, work);
            else idct_transform(line, plan_x, fft_x, work);
            for (int x = 0; x < nx; ++x) values[y * nx + x] = line[x];
        }
        #pragma omp for schedule(static)
        for (int x = 0; x < nx; ++x) {
            line.resize(ny);
            for (int y = 0; y < ny; ++y) line[y] = values[y * nx + x];
            if (sine_x) idct_transform(line, plan_y, fft_y, work);
            else inverse_sine_transform(line, plan_y, fft_y, reversed, work);
            for (int y = 0; y < ny; ++y) values[y * nx + x] = line[y];
        }
    }
}

}  // namespace dpcpp
