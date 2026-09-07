#include "microkernel.hpp"
#include "epsilon_active/bookshelf.hpp"

namespace nsgp::modules {
void register_layout_init(ModuleRegistry& registry) {
    registry.add("layout_init", [](StageContext& context, const Json& config) {
        const std::string mode = config.value("mode", std::string("raw"));
        if (mode == "center_gaussian") {
            ea::initialize_center_gaussian(
                context.db, config.value("seed", std::uint64_t(219)),
                config.value("sigma_ratio", .001));
        } else if (mode != "raw") {
            throw std::invalid_argument("layout_init mode must be raw or center_gaussian");
        }
        return StageStats{};
    });
}
}
