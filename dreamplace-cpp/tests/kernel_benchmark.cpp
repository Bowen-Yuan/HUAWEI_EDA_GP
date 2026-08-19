#include "bookshelf.h"
#include "electric.h"
#include "wirelength.h"

#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace dpcpp;

namespace {
double measure_ms(int repeats, const std::function<void()>& operation) {
    operation();
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) operation();
    return std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count() / repeats;
}
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: kernel_benchmark <benchmark-base> [bins] [repeats]\n";
        return 2;
    }
    const int bins = argc > 2 ? std::atoi(argv[2]) : 512;
    const int repeats = argc > 3 ? std::atoi(argv[3]) : 10;
    Database db = read_bookshelf(argv[1]);
    center_gaussian_initialize(db, 1000, 0.001);
    std::vector<Filler> fillers = initialize_fillers(db, 1.0, 1000);
    ElectricDensity density(db, bins, bins, 1.0);
    Database fixed_only = db;
    fixed_only.movable_ids.clear();
    ElectricDensity fixed_density(fixed_only, bins, bins, 1.0);
    const std::vector<Filler> no_fillers;
    std::vector<Real> gx, gy, fgx, fgy;
    volatile Real sink = 0.0;

    const double exact_gradient = measure_ms(repeats, [&]() {
        sink += exact_hpwl_subgradient(db, 100, &gx, &gy);
    });
    const double exact_value = measure_ms(repeats, [&]() {
        sink += exact_hpwl(db, 100);
    });
    const double active_gradient = measure_ms(repeats, [&]() {
        sink += exact_hpwl_active_set_direction(
            db, 100, 125.0, 4.0, &gx, &gy);
    });
    const double active_value = measure_ms(repeats, [&]() {
        sink += exact_hpwl_active_set_direction(
            db, 100, 125.0, 4.0, nullptr, nullptr);
    });
    const double wa_gradient = measure_ms(repeats, [&]() {
        sink += weighted_average_wirelength(db, 1000.0, 100, &gx, &gy);
    });
    const double wa_value = measure_ms(repeats, [&]() {
        sink += weighted_average_wirelength(db, 1000.0, 100, nullptr, nullptr);
    });
    const double density_gradient = measure_ms(repeats, [&]() {
        sink += density.compute(db, fillers, &gx, &gy, &fgx, &fgy).energy;
    });
    const double density_value = measure_ms(repeats, [&]() {
        sink += density.compute(db, fillers, nullptr, nullptr, nullptr, nullptr).energy;
    });
    const double density_nodes_only = measure_ms(repeats, [&]() {
        sink += density.compute(db, no_fillers, nullptr, nullptr, nullptr, nullptr).energy;
    });
    const double density_transform_only = measure_ms(repeats, [&]() {
        sink += fixed_density.compute(
            fixed_only, no_fillers, nullptr, nullptr, nullptr, nullptr).energy;
    });

    std::cout << std::fixed << std::setprecision(3)
              << "exact_gradient_ms=" << exact_gradient << '\n'
              << "exact_value_ms=" << exact_value << '\n'
              << "active_gradient_ms=" << active_gradient << '\n'
              << "active_value_ms=" << active_value << '\n'
              << "wa_gradient_ms=" << wa_gradient << '\n'
              << "wa_value_ms=" << wa_value << '\n'
              << "density_gradient_ms=" << density_gradient << '\n'
              << "density_value_ms=" << density_value << '\n'
              << "density_nodes_only_ms=" << density_nodes_only << '\n'
              << "density_transform_only_ms=" << density_transform_only << '\n'
              << "sink=" << sink << '\n';
    return 0;
}
