#pragma once

#include "optimizer.h"

namespace dpcpp {

GlobalPlaceResult nonsmooth_global_place(Database& db,
                                         std::vector<Filler>& fillers,
                                         const GlobalPlaceConfig& config);

}  // namespace dpcpp
