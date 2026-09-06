#include "epsilon_active/bookshelf.hpp"
#include "epsilon_active/density.hpp"
#include "epsilon_active/hpwl.hpp"

#include "electric.h"
#include "bookshelf.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {

using ea::Real;

struct Options {
    std::filesystem::path benchmark;
    std::filesystem::path initial_placement;
    std::filesystem::path output;
    std::uint64_t seed = 1000;
    Real sigma_ratio = 0.001;
    int bins = 128;
    int stages = 6;
    int stage_iterations = 40;
    int snapshot_every = 10;
    int threads = 16;
    Real target_density = 1.0;
    Real lambda_overlap = 1.0;
    bool ramp_overlap = false;
    Real hpwl_epsilon = 125.0;
    Real hpwl_active_power = 1.0;
    int hpwl_degree_limit = -1;
    bool adaptive_weights = false;
    bool gradient_balanced_weights = false;
    int adaptive_warmup_iterations = 100;
    int adaptive_update_interval = 10;
    int stability_patience = 3;
    Real stability_tolerance = 1.0e-4;
    Real stable_mu_decay = 0.90;
    Real adaptive_mu_boost = 1.0;
    Real adaptive_mu_max = 16.0;
    Real adaptive_mu_handoff_overflow = 0.60;
    bool adaptive_handoff_on_overflow = false;
    bool adaptive_mu_rebound_guard = false;
    Real adaptive_lambda_step = 2.0;
    Real adaptive_lambda_max = 150.0;
    Real adaptive_lambda_target_overflow = 0.07;
    Real adaptive_lambda_deadband = 0.005;
    bool adaptive_final_lambda_band = false;
    bool select_feasible_hpwl = false;
    bool restore_feasible_hpwl_for_final = false;
    Real overlap_epsilon_bins = 0.0;
    Real mu_start = 4.0;
    Real mu_decay = 0.25;
    Real mu_floor = 0.0;
    bool force_final_lambda_max = false;
    Real anchor_weight = 1.0e-4;
    std::string anchor_mode = "per-node";
    Real net_batch_weight = 0.0;
    int net_batch_degree_limit = 64;
    int net_batch_refresh = 1;
    Real step_fraction = 0.002;
    Real beta1 = 0.9;
    Real beta2 = 0.99;
    Real optimizer_epsilon = 1.0e-8;
    std::string optimizer = "adam";
};

std::string value_after(int& index, int argc, char** argv) {
    if (++index >= argc) throw std::invalid_argument("missing option value");
    return argv[index];
}

Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--benchmark") o.benchmark = value_after(i, argc, argv);
        else if (arg == "--initial-placement") o.initial_placement = value_after(i, argc, argv);
        else if (arg == "--output") o.output = value_after(i, argc, argv);
        else if (arg == "--seed") o.seed = std::stoull(value_after(i, argc, argv));
        else if (arg == "--sigma-ratio") o.sigma_ratio = std::stod(value_after(i, argc, argv));
        else if (arg == "--bins") o.bins = std::stoi(value_after(i, argc, argv));
        else if (arg == "--stages") o.stages = std::stoi(value_after(i, argc, argv));
        else if (arg == "--stage-iterations") o.stage_iterations = std::stoi(value_after(i, argc, argv));
        else if (arg == "--snapshot-every") o.snapshot_every = std::stoi(value_after(i, argc, argv));
        else if (arg == "--threads") o.threads = std::stoi(value_after(i, argc, argv));
        else if (arg == "--lambda-overlap") o.lambda_overlap = std::stod(value_after(i, argc, argv));
        else if (arg == "--ramp-overlap") o.ramp_overlap = true;
        else if (arg == "--hpwl-epsilon") o.hpwl_epsilon = std::stod(value_after(i, argc, argv));
        else if (arg == "--hpwl-active-power") o.hpwl_active_power = std::stod(value_after(i, argc, argv));
        else if (arg == "--hpwl-degree-limit") o.hpwl_degree_limit = std::stoi(value_after(i, argc, argv));
        else if (arg == "--adaptive-weights") o.adaptive_weights = true;
        else if (arg == "--gradient-balanced-weights") o.gradient_balanced_weights = true;
        else if (arg == "--adaptive-warmup-iterations") o.adaptive_warmup_iterations = std::stoi(value_after(i, argc, argv));
        else if (arg == "--adaptive-update-interval") o.adaptive_update_interval = std::stoi(value_after(i, argc, argv));
        else if (arg == "--stability-patience") o.stability_patience = std::stoi(value_after(i, argc, argv));
        else if (arg == "--stability-tolerance") o.stability_tolerance = std::stod(value_after(i, argc, argv));
        else if (arg == "--stable-mu-decay") o.stable_mu_decay = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-mu-boost") o.adaptive_mu_boost = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-mu-max") o.adaptive_mu_max = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-mu-handoff-overflow") o.adaptive_mu_handoff_overflow = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-handoff-on-overflow") o.adaptive_handoff_on_overflow = true;
        else if (arg == "--adaptive-mu-rebound-guard") o.adaptive_mu_rebound_guard = true;
        else if (arg == "--adaptive-lambda-step") o.adaptive_lambda_step = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-lambda-max") o.adaptive_lambda_max = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-lambda-target-overflow") o.adaptive_lambda_target_overflow = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-lambda-deadband") o.adaptive_lambda_deadband = std::stod(value_after(i, argc, argv));
        else if (arg == "--adaptive-final-lambda-band") o.adaptive_final_lambda_band = true;
        else if (arg == "--select-feasible-hpwl") o.select_feasible_hpwl = true;
        else if (arg == "--restore-feasible-hpwl-for-final")
            o.restore_feasible_hpwl_for_final = true;
        else if (arg == "--overlap-epsilon-bins") o.overlap_epsilon_bins = std::stod(value_after(i, argc, argv));
        else if (arg == "--mu-start") o.mu_start = std::stod(value_after(i, argc, argv));
        else if (arg == "--mu-decay") o.mu_decay = std::stod(value_after(i, argc, argv));
        else if (arg == "--mu-floor") o.mu_floor = std::stod(value_after(i, argc, argv));
        else if (arg == "--force-final-lambda-max") o.force_final_lambda_max = true;
        else if (arg == "--anchor-weight") o.anchor_weight = std::stod(value_after(i, argc, argv));
        else if (arg == "--anchor-mode") o.anchor_mode = value_after(i, argc, argv);
        else if (arg == "--net-batch-weight") o.net_batch_weight = std::stod(value_after(i, argc, argv));
        else if (arg == "--net-batch-degree-limit") o.net_batch_degree_limit = std::stoi(value_after(i, argc, argv));
        else if (arg == "--net-batch-refresh") o.net_batch_refresh = std::stoi(value_after(i, argc, argv));
        else if (arg == "--step-fraction") o.step_fraction = std::stod(value_after(i, argc, argv));
        else if (arg == "--optimizer") o.optimizer = value_after(i, argc, argv);
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if (o.benchmark.empty() || o.output.empty() || o.sigma_ratio < 0.0 ||
        o.bins <= 0 || (o.bins & (o.bins - 1)) != 0 || o.stages <= 0 ||
        o.stage_iterations <= 0 || o.threads <= 0 || o.threads > 40 ||
        o.snapshot_every <= 0 ||
        o.lambda_overlap < 0.0 || o.hpwl_epsilon < 0.0 ||
        o.hpwl_active_power <= 0.0 || o.overlap_epsilon_bins != 0.0 ||
        o.mu_start < 0.0 || o.mu_decay <= 0.0 ||
        o.mu_decay >= 1.0 || o.mu_floor < 0.0 || o.anchor_weight < 0.0 ||
        o.net_batch_weight < 0.0 || o.net_batch_degree_limit < 2 ||
        o.net_batch_refresh <= 0 ||
        o.step_fraction <= 0.0 || o.adaptive_warmup_iterations < 0 ||
        o.adaptive_update_interval <= 0 || o.stability_patience <= 0 ||
        o.stability_tolerance < 0.0 || o.stable_mu_decay <= 0.0 ||
        o.stable_mu_decay >= 1.0 || o.adaptive_lambda_step < 0.0 ||
        o.adaptive_mu_boost < 1.0 || o.adaptive_mu_max < o.mu_start ||
        o.adaptive_mu_handoff_overflow < 0.0 ||
        o.adaptive_mu_handoff_overflow >= 1.0 ||
        o.adaptive_lambda_max < 0.0 || o.adaptive_lambda_target_overflow < 0.0 ||
        o.adaptive_lambda_target_overflow >= 1.0 ||
        o.adaptive_lambda_deadband < 0.0 ||
        o.adaptive_lambda_deadband > o.adaptive_lambda_target_overflow ||
        (o.anchor_mode != "per-node" && o.anchor_mode != "centroid" &&
         o.anchor_mode != "initial") ||
        (o.optimizer != "adam" && o.optimizer != "sgd")) {
        throw std::invalid_argument("invalid homotopy options");
    }
    return o;
}

struct ExactOverlapResult {
    Real energy = 0.0;
    Real overflow = 0.0;
    Real max_density = 0.0;
};

class ExactOverlapOracle {
public:
    ExactOverlapOracle(const ea::Database& db, int bins, Real target)
        : db_(db), nx_(bins), ny_(bins), target_(target),
          bin_w_((db.xh - db.xl) / bins), bin_h_((db.yh - db.yl) / bins),
          bin_area_(bin_w_ * bin_h_), fixed_(static_cast<std::size_t>(bins) * bins),
          occupancy_(fixed_.size()) {
        for (int id : db_.fixed_ids) {
            const ea::Node& node = db_.nodes[id];
            if (!node.terminal_ni) deposit(node, node.x, node.y, fixed_);
        }
    }

    ExactOverlapResult evaluate(std::vector<Real>* gx, std::vector<Real>* gy) {
        occupancy_ = fixed_;
        for (int id : db_.movable_ids) {
            const ea::Node& node = db_.nodes[id];
            deposit(node, node.x, node.y, occupancy_);
        }
        ExactOverlapResult result;
        Real excess_area = 0.0;
        for (std::size_t i = 0; i < occupancy_.size(); ++i) {
            const Real density = occupancy_[i] / bin_area_;
            const Real excess = std::max<Real>(density - target_, 0.0);
            result.max_density = std::max(result.max_density, density);
            excess_area += excess * bin_area_;
            result.energy += 0.5 * bin_area_ * excess * excess;
        }
        result.overflow = excess_area / std::max<Real>(db_.movable_area, 1.0e-30);
        if (!gx && !gy) return result;
        if (gx) gx->assign(db_.nodes.size(), 0.0);
        if (gy) gy->assign(db_.nodes.size(), 0.0);
        #pragma omp parallel for schedule(static)
        for (int m = 0; m < static_cast<int>(db_.movable_ids.size()); ++m) {
            const int id = db_.movable_ids[m];
            const ea::Node& node = db_.nodes[id];
            Real xg = 0.0, yg = 0.0;
            const int x0 = lower(node.x - 0.5 * node.width);
            const int x1 = upper(node.x + 0.5 * node.width);
            const int y0 = lower(node.y - 0.5 * node.height);
            const int y1 = upper(node.y + 0.5 * node.height);
            for (int by = y0; by <= y1; ++by) {
                const Real bl = db_.yl + by * bin_h_;
                const Real oy = overlap(node.y, 0.5 * node.height, bl, bl + bin_h_);
                const Real dy = derivative(node.y, 0.5 * node.height, bl, bl + bin_h_);
                for (int bx = x0; bx <= x1; ++bx) {
                    const Real left = db_.xl + bx * bin_w_;
                    const Real ox = overlap(node.x, 0.5 * node.width, left, left + bin_w_);
                    const Real dx = derivative(node.x, 0.5 * node.width, left, left + bin_w_);
                    const Real excess = std::max<Real>(
                        occupancy_[static_cast<std::size_t>(by) * nx_ + bx] /
                            bin_area_ - target_, 0.0);
                    xg += excess * dx * oy;
                    yg += excess * dy * ox;
                }
            }
            if (gx) (*gx)[id] = xg;
            if (gy) (*gy)[id] = yg;
        }
        return result;
    }

private:
    int lower(Real coordinate) const {
        return std::clamp(static_cast<int>(std::floor((coordinate - db_.xl) / bin_w_)), 0, nx_ - 1);
    }
    int upper(Real coordinate) const {
        return std::clamp(static_cast<int>(std::floor(
            std::nextafter(coordinate, coordinate - 1.0) - db_.xl) / bin_w_), 0, nx_ - 1);
    }
    static Real overlap(Real center, Real half, Real low, Real high) {
        return std::max<Real>(0.0, std::min(center + half, high) - std::max(center - half, low));
    }
    static Real derivative(Real center, Real half, Real low, Real high) {
        const Real lo = center - half, hi = center + half;
        if (hi <= low || lo >= high) return 0.0;
        return (hi < high ? 1.0 : 0.0) - (lo > low ? 1.0 : 0.0);
    }
    void deposit(const ea::Node& node, Real x, Real y, std::vector<Real>& map) const {
        const Real left = x - 0.5 * node.width, right = x + 0.5 * node.width;
        const Real bottom = y - 0.5 * node.height, top = y + 0.5 * node.height;
        const int x0 = lower(left), x1 = upper(right);
        const int y0 = std::clamp(static_cast<int>(std::floor((bottom - db_.yl) / bin_h_)), 0, ny_ - 1);
        const int y1 = std::clamp(static_cast<int>(std::floor(
            std::nextafter(top, bottom) - db_.yl) / bin_h_), 0, ny_ - 1);
        for (int by = y0; by <= y1; ++by) {
            const Real bl = db_.yl + by * bin_h_;
            const Real oy = std::max<Real>(0.0, std::min(top, bl + bin_h_) - std::max(bottom, bl));
            for (int bx = x0; bx <= x1; ++bx) {
                const Real left_bin = db_.xl + bx * bin_w_;
                const Real ox = std::max<Real>(0.0, std::min(right, left_bin + bin_w_) - std::max(left, left_bin));
                map[static_cast<std::size_t>(by) * nx_ + bx] += ox * oy;
            }
        }
    }
    const ea::Database& db_;
    int nx_, ny_;
    Real target_, bin_w_, bin_h_, bin_area_;
    std::vector<Real> fixed_, occupancy_;
};

struct ElectricBridge {
    dpcpp::Database db;
    dpcpp::ElectricDensity field;

    ElectricBridge(const ea::Database& source, int bins, Real target)
        : db(make_database(source)), field(db, bins, bins, target,
                       dpcpp::ElectricFieldModel::FiniteDifference) {
    }

    static dpcpp::Database make_database(const ea::Database& source) {
        dpcpp::Database result;
        result.xl = source.xl; result.yl = source.yl;
        result.xh = source.xh; result.yh = source.yh;
        result.nodes.resize(source.nodes.size());
        for (std::size_t i = 0; i < source.nodes.size(); ++i) {
            const ea::Node& n = source.nodes[i];
            result.nodes[i].id = n.id; result.nodes[i].x = n.x; result.nodes[i].y = n.y;
            result.nodes[i].width = n.width; result.nodes[i].height = n.height;
            result.nodes[i].fixed = n.fixed; result.nodes[i].terminal_ni = n.terminal_ni;
        }
        result.movable_ids = source.movable_ids;
        result.fixed_ids = source.fixed_ids;
        return result;
    }

    void sync(const ea::Database& source) {
        for (std::size_t i = 0; i < source.nodes.size(); ++i) {
            db.nodes[i].x = source.nodes[i].x;
            db.nodes[i].y = source.nodes[i].y;
        }
    }
};

std::vector<Real> capture(const ea::Database& db) {
    const std::size_t n = db.movable_ids.size();
    std::vector<Real> result(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        const ea::Node& node = db.nodes[db.movable_ids[i]];
        result[i] = node.x; result[n + i] = node.y;
    }
    return result;
}

void apply_positions(ea::Database& db, const std::vector<Real>& values) {
    const std::size_t n = db.movable_ids.size();
    for (std::size_t i = 0; i < n; ++i) {
        ea::Node& node = db.nodes[db.movable_ids[i]];
        node.x = std::clamp(values[i], db.xl + 0.5 * node.width, db.xh - 0.5 * node.width);
        node.y = std::clamp(values[n + i], db.yl + 0.5 * node.height, db.yh - 0.5 * node.height);
    }
}

Real norm(const std::vector<Real>& values) {
    Real sum = 0.0;
    for (Real value : values) sum += value * value;
    return std::sqrt(sum);
}

void adam_step(std::vector<Real>& values, const std::vector<Real>& gradient,
               std::vector<Real>& first, std::vector<Real>& second,
               int step, Real learning_rate, Real beta1, Real beta2, Real eps,
               const ea::Database& db) {
    const std::size_t n = db.movable_ids.size();
    const Real bias1 = 1.0 - std::pow(beta1, step);
    const Real bias2 = 1.0 - std::pow(beta2, step);
    for (std::size_t i = 0; i < values.size(); ++i) {
        first[i] = beta1 * first[i] + (1.0 - beta1) * gradient[i];
        second[i] = beta2 * second[i] + (1.0 - beta2) * gradient[i] * gradient[i];
        const Real delta = learning_rate * (first[i] / bias1) /
            (std::sqrt(second[i] / bias2) + eps);
        values[i] -= delta;
    }
    (void)n;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        omp_set_dynamic(0);
        omp_set_num_threads(options.threads);
        std::filesystem::create_directories(options.output / "snapshots");

        ea::Database db = ea::read_bookshelf(options.benchmark);
        if (options.initial_placement.empty()) {
            ea::initialize_center_gaussian(db, options.seed, options.sigma_ratio);
        } else {
            std::string initial_name = options.initial_placement.string();
            std::transform(initial_name.begin(), initial_name.end(), initial_name.begin(),
                           [](unsigned char value) {
                               return static_cast<char>(std::tolower(value));
                           });
            const bool lg_placement = initial_name.size() >= 6 &&
                initial_name.compare(initial_name.size() - 6, 6, ".lg.pl") == 0;
            if (initial_name.find("eplace") != std::string::npos || lg_placement) {
                throw std::invalid_argument(
                    "ePlace and .lg.pl initial placements are disabled");
            }
            ea::load_bookshelf_placement(db, options.initial_placement);
        }
        const std::vector<Real> initial = capture(db);
        const Real die_scale = std::min(db.xh - db.xl, db.yh - db.yl);
        ExactOverlapOracle overlap(db, options.bins, options.target_density);
        ea::ExactOverlapDensity overlap_direction(
            db, options.bins, options.bins, options.target_density);
        ea::ExactHpwl hpwl(db);
        ElectricBridge electric(db, options.bins, options.target_density);

        std::vector<Real> ogx, ogy, hgx, hgy, egx, egy;
        ExactOverlapResult exact = overlap.evaluate(nullptr, nullptr);
        // The overlap direction is always the exact piecewise derivative. Any
        // active-set reach belongs exclusively to the HPWL subgradient.
        const Real overlap_epsilon = 0.0;
        overlap_direction.evaluate(overlap_epsilon, 1.0, &ogx, &ogy);
        const Real initial_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        electric.sync(db);
        std::vector<Real> electric_gx, electric_gy;
        std::vector<dpcpp::Filler> no_fillers;
        dpcpp::DensityResult electro = electric.field.compute(
            electric.db, no_fillers, &electric_gx, &electric_gy, nullptr, nullptr);
        const Real hpwl_scale = std::max<Real>(initial_hpwl, 1.0);
        const Real overlap_scale = std::max<Real>(exact.energy, 1.0);
        const Real electro_scale = std::max<Real>(electro.energy, 1.0);
        const Real area_scale = std::max<Real>(db.movable_area, 1.0);
        hpwl.evaluate(options.hpwl_epsilon, options.hpwl_active_power,
                      options.hpwl_degree_limit, &hgx, &hgy);
        std::vector<Real> initial_batch_x(db.nodes.size(), 0.0);
        std::vector<Real> initial_batch_y(db.nodes.size(), 0.0);
        if (options.net_batch_weight > 0.0) {
            std::vector<int> counts(db.nodes.size(), 0);
            for (const ea::Net& net : db.nets) {
                if (net.pin_count < 2 ||
                    static_cast<int>(net.pin_count) > options.net_batch_degree_limit) continue;
                Real mean_x = 0.0, mean_y = 0.0;
                int count = 0;
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int id = db.pins[p].node;
                    if (db.nodes[id].fixed) continue;
                    mean_x += hgx[id]; mean_y += hgy[id]; ++count;
                }
                if (count == 0) continue;
                mean_x /= static_cast<Real>(count);
                mean_y /= static_cast<Real>(count);
                for (std::size_t p = net.pin_begin;
                     p < net.pin_begin + net.pin_count; ++p) {
                    const int id = db.pins[p].node;
                    if (db.nodes[id].fixed) continue;
                    initial_batch_x[id] += mean_x;
                    initial_batch_y[id] += mean_y;
                    ++counts[id];
                }
            }
            for (int id : db.movable_ids) {
                if (counts[id] > 0) {
                    initial_batch_x[id] /= static_cast<Real>(counts[id]);
                    initial_batch_y[id] /= static_cast<Real>(counts[id]);
                }
            }
        }
        Real initial_hpwl_grad_sq = 0.0;
        Real initial_overlap_grad_sq = 0.0;
        Real initial_electro_grad_sq = 0.0;
        const Real anchor_x = db.xl + 0.5 * (db.xh - db.xl);
        const Real anchor_y = db.yl + 0.5 * (db.yh - db.yl);
        const Real anchor_scale = options.anchor_weight / area_scale;
        Real initial_centroid_x = 0.0;
        Real initial_centroid_y = 0.0;
        for (int id : db.movable_ids) {
            initial_centroid_x += db.nodes[id].x;
            initial_centroid_y += db.nodes[id].y;
        }
        const Real inverse_movable_count = 1.0 /
            std::max<std::size_t>(db.movable_ids.size(), 1);
        initial_centroid_x *= inverse_movable_count;
        initial_centroid_y *= inverse_movable_count;
        const Real initial_centroid_gx = anchor_scale *
            (initial_centroid_x - anchor_x) * inverse_movable_count;
        const Real initial_centroid_gy = anchor_scale *
            (initial_centroid_y - anchor_y) * inverse_movable_count;
        for (std::size_t m = 0; m < db.movable_ids.size(); ++m) {
            const int id = db.movable_ids[m];
            const Real batch_weight = std::clamp(options.net_batch_weight, 0.0, 1.0);
            const Real hx = ((1.0 - batch_weight) * hgx[id] +
                batch_weight * initial_batch_x[id]) / hpwl_scale;
            const Real hy = ((1.0 - batch_weight) * hgy[id] +
                batch_weight * initial_batch_y[id]) / hpwl_scale;
            const Real ox = ogx[id] / overlap_scale;
            const Real oy = ogy[id] / overlap_scale;
            const Real anchor_target_x = options.anchor_mode == "initial"
                ? initial[m] : anchor_x;
            const Real anchor_target_y = options.anchor_mode == "initial"
                ? initial[db.movable_ids.size() + m] : anchor_y;
            const Real anchor_gx = options.anchor_mode == "centroid"
                ? initial_centroid_gx
                : anchor_scale * (db.nodes[id].x - anchor_target_x);
            const Real anchor_gy = options.anchor_mode == "centroid"
                ? initial_centroid_gy
                : anchor_scale * (db.nodes[id].y - anchor_target_y);
            const Real ex = electric_gx[id] / electro_scale + anchor_gx;
            const Real ey = electric_gy[id] / electro_scale + anchor_gy;
            initial_hpwl_grad_sq += hx * hx + hy * hy;
            initial_overlap_grad_sq += ox * ox + oy * oy;
            initial_electro_grad_sq += ex * ex + ey * ey;
        }
        const Real initial_hpwl_grad_norm = std::sqrt(initial_hpwl_grad_sq);
        const Real overlap_gradient_scale = options.gradient_balanced_weights
            ? initial_hpwl_grad_norm /
                std::max<Real>(std::sqrt(initial_overlap_grad_sq), 1.0e-30)
            : 1.0;
        const Real electro_gradient_scale = options.gradient_balanced_weights
            ? initial_hpwl_grad_norm /
                std::max<Real>(std::sqrt(initial_electro_grad_sq), 1.0e-30)
            : 1.0;

        std::ofstream metrics(options.output / "homotopy_metrics.csv");
        if (!metrics) throw std::runtime_error("cannot create homotopy_metrics.csv");
        metrics << "iteration,stage,mu,lambda_stage,exact_hpwl,exact_overflow,exact_energy,"
                   "electro_energy,anchor_energy,electro_overflow,objective,grad_norm,seconds,"
                   "adaptive_started,stable_windows,hpwl_term_norm,"
                   "overlap_term_norm,electro_term_norm,overlap_active_fraction,"
                   "mu_control,lambda_control\n";
        metrics << std::setprecision(12);
        const auto started = std::chrono::steady_clock::now();
        std::vector<Real> values = initial;
        std::vector<Real> best_values = values;
        std::vector<Real> continuation_feasible_values;
        std::vector<Real> final_stage_best_values;
        std::vector<Real> final_stage_feasible_values;
        Real best_overflow = exact.overflow;
        Real best_hpwl = initial_hpwl;
        Real continuation_feasible_hpwl = std::numeric_limits<Real>::infinity();
        Real continuation_feasible_overflow = std::numeric_limits<Real>::infinity();
        Real final_stage_best_overflow = std::numeric_limits<Real>::infinity();
        Real final_stage_best_hpwl = std::numeric_limits<Real>::infinity();
        Real final_stage_feasible_overflow = std::numeric_limits<Real>::infinity();
        Real final_stage_feasible_hpwl = std::numeric_limits<Real>::infinity();
        std::vector<Real> first(values.size(), 0.0), second(values.size(), 0.0);
        const Real learning_rate = options.step_fraction * die_scale;
        int optimizer_step = 0;
        int iteration = 0;
        Real adaptive_mu = options.mu_start;
        Real adaptive_lambda = 0.0;
        Real adaptive_reference_overflow = exact.overflow;
        Real adaptive_window_best = exact.overflow;
        int stable_windows = 0;
        bool adaptive_started = false;
        Real adaptive_trigger_overflow = exact.overflow;
        std::vector<Real> batch_x;
        std::vector<Real> batch_y;

        for (int stage = 0; stage < options.stages; ++stage) {
            if (stage == options.stages - 1) {
                // Continue the exact limiting problem from the strongest
                // homotopy checkpoint, but never report that positive-mu
                // checkpoint as the final nonsmooth solution.
                const bool restore_feasible = options.restore_feasible_hpwl_for_final &&
                    !continuation_feasible_values.empty();
                values = restore_feasible ? continuation_feasible_values : best_values;
                apply_positions(db, values);
                // The final stage optimizes only the exact nonsmooth target.
                // Remove electrostatic history from Adam at the homotopy limit.
                std::fill(first.begin(), first.end(), 0.0);
                std::fill(second.begin(), second.end(), 0.0);
                optimizer_step = 0;
                if (options.force_final_lambda_max)
                    adaptive_lambda = options.adaptive_lambda_max;
            }
            Real mu = stage == options.stages - 1 && options.mu_floor == 0.0
                ? 0.0 : std::max(options.mu_floor,
                    options.mu_start * std::pow(options.mu_decay, stage));
            const Real lambda_stage = options.ramp_overlap
                ? options.lambda_overlap * static_cast<Real>(stage) /
                    std::max(1, options.stages - 1)
                : options.lambda_overlap;
            for (int local = 0; local < options.stage_iterations; ++local, ++iteration) {
                Real active_mu_control = options.adaptive_weights ? adaptive_mu : mu;
                if (stage == options.stages - 1) active_mu_control = 0.0;
                const Real active_lambda_control = options.adaptive_weights
                    ? adaptive_lambda : lambda_stage;
                const Real active_mu = active_mu_control * electro_gradient_scale;
                const Real active_lambda =
                    active_lambda_control * overlap_gradient_scale;
                apply_positions(db, values);
                electric.sync(db);
                exact = overlap.evaluate(nullptr, nullptr);
                overlap_direction.evaluate(overlap_epsilon, 1.0, &ogx, &ogy);
                const Real current_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
                // Epsilon-active is a subgradient selection for exact HPWL;
                // it does not alter current_hpwl or any audited objective.
                hpwl.evaluate(options.hpwl_epsilon, options.hpwl_active_power,
                               options.hpwl_degree_limit, &hgx, &hgy);
                std::vector<Real> electric_gx_local, electric_gy_local;
                electro = electric.field.compute(electric.db, no_fillers,
                                                  &electric_gx_local,
                                                  &electric_gy_local, nullptr, nullptr);
                if (exact.overflow < best_overflow ||
                    (std::abs(exact.overflow - best_overflow) < 1.0e-12 &&
                     current_hpwl < best_hpwl)) {
                    best_overflow = exact.overflow;
                    best_hpwl = current_hpwl;
                    best_values = capture(db);
                }
                if (stage != options.stages - 1 &&
                    exact.overflow <= options.adaptive_lambda_target_overflow &&
                    current_hpwl < continuation_feasible_hpwl) {
                    continuation_feasible_hpwl = current_hpwl;
                    continuation_feasible_overflow = exact.overflow;
                    continuation_feasible_values = capture(db);
                }
                if (stage == options.stages - 1 &&
                    (exact.overflow < final_stage_best_overflow ||
                     (std::abs(exact.overflow - final_stage_best_overflow) < 1.0e-12 &&
                      current_hpwl < final_stage_best_hpwl))) {
                    final_stage_best_overflow = exact.overflow;
                    final_stage_best_hpwl = current_hpwl;
                    final_stage_best_values = capture(db);
                }
                if (stage == options.stages - 1 &&
                    exact.overflow <= options.adaptive_lambda_target_overflow &&
                    current_hpwl < final_stage_feasible_hpwl) {
                    final_stage_feasible_overflow = exact.overflow;
                    final_stage_feasible_hpwl = current_hpwl;
                    final_stage_feasible_values = capture(db);
                }
                std::vector<Real> gradient(values.size(), 0.0);
                const std::size_t n = db.movable_ids.size();
                Real hpwl_term_sq = 0.0;
                Real overlap_term_sq = 0.0;
                Real electro_term_sq = 0.0;
                Real centroid_x = 0.0;
                Real centroid_y = 0.0;
                for (int id : db.movable_ids) {
                    centroid_x += db.nodes[id].x;
                    centroid_y += db.nodes[id].y;
                }
                centroid_x *= inverse_movable_count;
                centroid_y *= inverse_movable_count;
                const Real centroid_anchor_gx = anchor_scale *
                    (centroid_x - anchor_x) * inverse_movable_count;
                const Real centroid_anchor_gy = anchor_scale *
                    (centroid_y - anchor_y) * inverse_movable_count;
                Real anchor_energy = options.anchor_mode == "centroid"
                    ? 0.5 * anchor_scale *
                         ((centroid_x - anchor_x) * (centroid_x - anchor_x) +
                          (centroid_y - anchor_y) * (centroid_y - anchor_y))
                    : 0.0;
                std::size_t overlap_active = 0;
                if (options.net_batch_weight > 0.0 &&
                    (batch_x.empty() || iteration % options.net_batch_refresh == 0)) {
                    batch_x.assign(db.nodes.size(), 0.0);
                    batch_y.assign(db.nodes.size(), 0.0);
                    std::vector<int> counts(db.nodes.size(), 0);
                    for (const ea::Net& net : db.nets) {
                        if (net.pin_count < 2 ||
                            static_cast<int>(net.pin_count) > options.net_batch_degree_limit) continue;
                        Real mean_x = 0.0, mean_y = 0.0;
                        int count = 0;
                        for (std::size_t p = net.pin_begin;
                             p < net.pin_begin + net.pin_count; ++p) {
                            const int pin_id = db.pins[p].node;
                            if (db.nodes[pin_id].fixed) continue;
                            mean_x += hgx[pin_id]; mean_y += hgy[pin_id]; ++count;
                        }
                        if (count == 0) continue;
                        mean_x /= static_cast<Real>(count);
                        mean_y /= static_cast<Real>(count);
                        for (std::size_t p = net.pin_begin;
                             p < net.pin_begin + net.pin_count; ++p) {
                            const int pin_id = db.pins[p].node;
                            if (db.nodes[pin_id].fixed) continue;
                            batch_x[pin_id] += mean_x;
                            batch_y[pin_id] += mean_y;
                            ++counts[pin_id];
                        }
                    }
                    for (int pin_id : db.movable_ids) {
                        if (counts[pin_id] > 0) {
                            batch_x[pin_id] /= static_cast<Real>(counts[pin_id]);
                            batch_y[pin_id] /= static_cast<Real>(counts[pin_id]);
                        }
                    }
                }
                for (std::size_t m = 0; m < n; ++m) {
                    const int id = db.movable_ids[m];
                const Real batch_weight = std::clamp(options.net_batch_weight, 0.0, 1.0);
                const Real hx = ((1.0 - batch_weight) * hgx[id] + batch_weight *
                    (batch_x.empty() ? 0.0 : batch_x[id])) / hpwl_scale;
                const Real hy = ((1.0 - batch_weight) * hgy[id] + batch_weight *
                    (batch_y.empty() ? 0.0 : batch_y[id])) / hpwl_scale;
                    const Real ox = active_lambda * ogx[id] / overlap_scale;
                    const Real oy = active_lambda * ogy[id] / overlap_scale;
                    const Real anchor_target_x = options.anchor_mode == "initial"
                        ? initial[m] : anchor_x;
                    const Real anchor_target_y = options.anchor_mode == "initial"
                        ? initial[n + m] : anchor_y;
                    const Real anchor_gx = options.anchor_mode == "centroid"
                        ? centroid_anchor_gx
                        : anchor_scale * (db.nodes[id].x - anchor_target_x);
                    const Real anchor_gy = options.anchor_mode == "centroid"
                        ? centroid_anchor_gy
                        : anchor_scale * (db.nodes[id].y - anchor_target_y);
                    if (options.anchor_mode != "centroid") {
                        anchor_energy += 0.5 * anchor_scale *
                            ((db.nodes[id].x - anchor_target_x) *
                                 (db.nodes[id].x - anchor_target_x) +
                             (db.nodes[id].y - anchor_target_y) *
                                 (db.nodes[id].y - anchor_target_y));
                    }
                    const Real ex = active_mu *
                        (electric_gx_local[id] / electro_scale + anchor_gx);
                    const Real ey = active_mu *
                        (electric_gy_local[id] / electro_scale + anchor_gy);
                    gradient[m] = hx + ox + ex;
                    gradient[n + m] = hy + oy + ey;
                    hpwl_term_sq += hx * hx + hy * hy;
                    overlap_term_sq += ox * ox + oy * oy;
                    electro_term_sq += ex * ex + ey * ey;
                    overlap_active += (ogx[id] != 0.0 || ogy[id] != 0.0) ? 1 : 0;
                }
                // Metrics describe the state before this optimizer update.
                // Capture the same coordinates so iter_N.pl is a valid,
                // reproducible continuation point for that metrics row.
                if (iteration % options.snapshot_every == 0) {
                    ea::write_bookshelf_placement(db,
                        options.output / "snapshots" /
                        ("iter_" + std::to_string(iteration) + ".pl"));
                }
                ++optimizer_step;
                if (options.optimizer == "adam") {
                    adam_step(values, gradient, first, second, optimizer_step,
                              learning_rate, options.beta1, options.beta2,
                              options.optimizer_epsilon, db);
                } else {
                    for (std::size_t i = 0; i < values.size(); ++i)
                        values[i] -= learning_rate * gradient[i];
                }
                apply_positions(db, values);
                const double seconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - started).count();
                const Real objective = current_hpwl / hpwl_scale +
                    active_lambda * exact.energy / overlap_scale +
                    active_mu * (electro.energy / electro_scale + anchor_energy);
                metrics << iteration << ',' << stage << ',' << active_mu << ','
                        << active_lambda << ','
                        << current_hpwl << ',' << exact.overflow << ','
                        << exact.energy << ',' << electro.energy << ','
                        << anchor_energy << ',' << electro.overflow << ',' << objective << ','
                        << norm(gradient) << ',' << seconds << ','
                        << (adaptive_started ? 1 : 0) << ',' << stable_windows << ','
                        << std::sqrt(hpwl_term_sq) << ','
                        << std::sqrt(overlap_term_sq) << ','
                        << std::sqrt(electro_term_sq) << ','
                        << static_cast<Real>(overlap_active) /
                            std::max<std::size_t>(n, 1) << ','
                        << active_mu_control << ',' << active_lambda_control << '\n';
                if (iteration % 10 == 0 || local + 1 == options.stage_iterations) {
                    std::cout << "[Homotopy] iter=" << iteration
                              << " stage=" << stage << " mu=" << active_mu
                              << " mu_control=" << active_mu_control
                              << " lambda=" << active_lambda
                              << " lambda_control=" << active_lambda_control
                              << " hpwl=" << current_hpwl
                              << " exact_overflow=" << exact.overflow
                              << " electro_energy=" << electro.energy
                              << " grad=" << norm(gradient)
                              << " seconds=" << seconds << '\n';
                }
                if (options.adaptive_weights) {
                    adaptive_window_best = std::min(adaptive_window_best, exact.overflow);
                    const bool update = iteration + 1 >= options.adaptive_warmup_iterations &&
                        (iteration + 1 - options.adaptive_warmup_iterations) %
                            options.adaptive_update_interval == 0;
                    if (update) {
                        const Real improvement = adaptive_reference_overflow - adaptive_window_best;
                        const bool stable = improvement <= options.stability_tolerance;
                        stable_windows = stable ? stable_windows + 1 : 0;
                        if (!adaptive_started &&
                            options.adaptive_handoff_on_overflow &&
                            adaptive_window_best <=
                                options.adaptive_mu_handoff_overflow) {
                            adaptive_started = true;
                            adaptive_trigger_overflow = adaptive_window_best;
                            stable_windows = 0;
                        }
                        if (!adaptive_started &&
                            stable_windows >= options.stability_patience) {
                            const bool high_overflow = adaptive_window_best >
                                options.adaptive_mu_handoff_overflow;
                            const bool can_boost = options.adaptive_mu_boost > 1.0 &&
                                adaptive_mu < options.adaptive_mu_max;
                            if (high_overflow && can_boost) {
                                adaptive_mu = std::min(options.adaptive_mu_max,
                                    adaptive_mu * options.adaptive_mu_boost);
                                stable_windows = 0;
                            } else {
                                adaptive_started = true;
                                adaptive_trigger_overflow = adaptive_window_best;
                            }
                        }
                        if (options.adaptive_final_lambda_band &&
                            stage == options.stages - 1) {
                            const Real lower = std::max<Real>(0.0,
                                options.adaptive_lambda_target_overflow -
                                    options.adaptive_lambda_deadband);
                            if (exact.overflow >
                                options.adaptive_lambda_target_overflow) {
                                adaptive_lambda = std::min(
                                    options.adaptive_lambda_max,
                                    adaptive_lambda + options.adaptive_lambda_step);
                            } else if (exact.overflow < lower) {
                                adaptive_lambda = std::max<Real>(0.0,
                                    adaptive_lambda - options.adaptive_lambda_step);
                            }
                        } else if (adaptive_started &&
                                   stage != options.stages - 1) {
                            const bool density_rebounded =
                                options.adaptive_mu_rebound_guard &&
                                improvement < -options.stability_tolerance;
                            adaptive_mu = density_rebounded
                                ? std::min(options.adaptive_mu_max,
                                    adaptive_mu / options.stable_mu_decay)
                                : std::max(options.mu_floor,
                                    adaptive_mu * options.stable_mu_decay);
                            const Real denominator = std::max<Real>(1.0e-12,
                                adaptive_trigger_overflow -
                                    options.adaptive_lambda_target_overflow);
                            const Real normalized_error = std::clamp(
                                (adaptive_window_best -
                                    options.adaptive_lambda_target_overflow) /
                                    denominator,
                                0.0, 1.0);
                            adaptive_lambda = std::min(options.adaptive_lambda_max,
                                adaptive_lambda + options.adaptive_lambda_step *
                                    normalized_error);
                        }
                        adaptive_reference_overflow = adaptive_window_best;
                        adaptive_window_best = exact.overflow;
                    }
                }
            }
        }

        if (final_stage_best_values.empty())
            throw std::runtime_error("final mu-zero stage produced no checkpoint");
        const bool selected_feasible_hpwl = options.select_feasible_hpwl &&
            !final_stage_feasible_values.empty();
        apply_positions(db, selected_feasible_hpwl
            ? final_stage_feasible_values : final_stage_best_values);
        exact = overlap.evaluate(nullptr, nullptr);
        const Real final_hpwl = hpwl.evaluate(0.0, 1.0, -1, nullptr, nullptr);
        ea::write_bookshelf_placement(db, options.output / "best.pl");
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        std::ofstream summary(options.output / "summary.txt");
        summary << std::setprecision(12)
                << "mode=exact_hpwl_plus_exact_overlap_plus_electrostatic_homotopy\n"
                << "initialization_mode="
                << (options.initial_placement.empty()
                        ? "fresh-center-gaussian" : "explicit-stage-continuation")
                << '\n'
                << "initial_placement=" << options.initial_placement.string() << '\n'
                << "initialization_seed=" << options.seed << '\n'
                << "initialization_sigma_ratio=" << options.sigma_ratio << '\n'
                << "benchmark_raw_pl=" << db.raw_pl_path << '\n'
                << "optimizer=" << options.optimizer << '\n'
                << "bins=" << options.bins << 'x' << options.bins << '\n'
                << "stages=" << options.stages << '\n'
                << "stage_iterations=" << options.stage_iterations << '\n'
                << "snapshot_every=" << options.snapshot_every << '\n'
                << "mu_start=" << options.mu_start << '\n'
                << "mu_decay=" << options.mu_decay << '\n'
                << "lambda_overlap=" << options.lambda_overlap << '\n'
                << "ramp_overlap=" << (options.ramp_overlap ? "true" : "false") << '\n'
                << "hpwl_epsilon=" << options.hpwl_epsilon << '\n'
                << "hpwl_active_power=" << options.hpwl_active_power << '\n'
                << "hpwl_degree_limit=" << options.hpwl_degree_limit << '\n'
                << "net_batch_weight=" << options.net_batch_weight << '\n'
                << "net_batch_degree_limit=" << options.net_batch_degree_limit << '\n'
                << "net_batch_refresh=" << options.net_batch_refresh << '\n'
                << "adaptive_weights=" << (options.adaptive_weights ? "true" : "false") << '\n'
                << "gradient_balanced_weights="
                << (options.gradient_balanced_weights ? "true" : "false") << '\n'
                << "overlap_gradient_scale=" << overlap_gradient_scale << '\n'
                << "electro_gradient_scale=" << electro_gradient_scale << '\n'
                << "adaptive_warmup_iterations=" << options.adaptive_warmup_iterations << '\n'
                << "adaptive_update_interval=" << options.adaptive_update_interval << '\n'
                << "stability_patience=" << options.stability_patience << '\n'
                << "stability_tolerance=" << options.stability_tolerance << '\n'
                << "stable_mu_decay=" << options.stable_mu_decay << '\n'
                << "adaptive_mu_boost=" << options.adaptive_mu_boost << '\n'
                << "adaptive_mu_max=" << options.adaptive_mu_max << '\n'
                << "adaptive_mu_handoff_overflow="
                << options.adaptive_mu_handoff_overflow << '\n'
                << "adaptive_handoff_on_overflow="
                << (options.adaptive_handoff_on_overflow ? "true" : "false") << '\n'
                << "adaptive_mu_rebound_guard="
                << (options.adaptive_mu_rebound_guard ? "true" : "false") << '\n'
                << "adaptive_lambda_step=" << options.adaptive_lambda_step << '\n'
                << "adaptive_lambda_max=" << options.adaptive_lambda_max << '\n'
                << "adaptive_lambda_target_overflow="
                << options.adaptive_lambda_target_overflow << '\n'
                << "adaptive_lambda_deadband="
                << options.adaptive_lambda_deadband << '\n'
                << "adaptive_final_lambda_band="
                << (options.adaptive_final_lambda_band ? "true" : "false") << '\n'
                << "select_feasible_hpwl="
                << (options.select_feasible_hpwl ? "true" : "false") << '\n'
                << "restore_feasible_hpwl_for_final="
                << (options.restore_feasible_hpwl_for_final ? "true" : "false") << '\n'
                << "final_mu=0 (forced in final stage)\n"
                << "force_final_lambda_max="
                << (options.force_final_lambda_max ? "true" : "false") << '\n'
                << "overlap_epsilon_bins=" << options.overlap_epsilon_bins << '\n'
                << "anchor_weight=" << options.anchor_weight << '\n'
                << "anchor_mode=" << options.anchor_mode << '\n'
                << "step_fraction=" << options.step_fraction << '\n'
                << "final_optimizer_state_reset=true\n"
                << "initial_hpwl=" << initial_hpwl << '\n'
                << "continuation_best_hpwl=" << best_hpwl << '\n'
                << "continuation_best_overflow=" << best_overflow << '\n'
                << "continuation_feasible_hpwl=" << continuation_feasible_hpwl << '\n'
                << "continuation_feasible_overflow="
                << continuation_feasible_overflow << '\n'
                << "final_stage_best_hpwl=" << final_stage_best_hpwl << '\n'
                << "final_stage_best_overflow=" << final_stage_best_overflow << '\n'
                << "final_stage_feasible_hpwl=" << final_stage_feasible_hpwl << '\n'
                << "final_stage_feasible_overflow="
                << final_stage_feasible_overflow << '\n'
                << "selected_feasible_hpwl="
                << (selected_feasible_hpwl ? "true" : "false") << '\n'
                << "final_hpwl=" << final_hpwl << '\n'
                << "final_exact_overflow=" << exact.overflow << '\n'
                << "wall_seconds=" << seconds << '\n'
                << "fixed_macro_capacity=true\n"
                << "legalization=disabled\n"
                << "exact_overlap=true\n"
                << "strict_final_stage_only=true\n"
                << "electrostatic_dct=true\n";
        std::cout << "[Result] hpwl=" << final_hpwl
                  << " exact_overflow=" << exact.overflow
                  << " best_overflow=" << best_overflow
                  << " seconds=" << seconds << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "homotopy failure: " << error.what() << '\n';
        return 1;
    }
}
