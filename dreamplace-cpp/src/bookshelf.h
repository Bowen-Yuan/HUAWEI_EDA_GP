#pragma once

#include "types.h"

#include <cstdint>
#include <string>

namespace dpcpp {

Database read_bookshelf(const std::string& benchmark_base);
void center_gaussian_initialize(Database& db, std::uint64_t seed,
                                Real sigma_ratio);
void write_bookshelf_pl(const Database& db, const std::string& path);
std::string basename_of(const std::string& path);

}  // namespace dpcpp

