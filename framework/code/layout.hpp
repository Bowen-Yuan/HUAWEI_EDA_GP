#pragma once
#include "common.hpp"
#include <vector>
namespace nsgp { struct Layout { std::vector<Real> x,y; std::vector<unsigned char> orientation; std::uint64_t revision=0; void touch(){++revision;} }; }
