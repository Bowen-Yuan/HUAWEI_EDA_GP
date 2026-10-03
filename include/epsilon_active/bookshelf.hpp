#pragma once

#include "epsilon_active/types.hpp"

#include <cstdint>
#include <filesystem>

namespace ea {

Database read_bookshelf(const std::filesystem::path& benchmark);
void load_bookshelf_placement(Database& db, const std::filesystem::path& path);
void initialize_center_gaussian(Database& db, std::uint64_t seed,
                                Real sigma_ratio);
void write_bookshelf_placement(const Database& db,
                               const std::filesystem::path& path);

}  // namespace ea
