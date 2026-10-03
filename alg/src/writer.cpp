#include "types.h"
#include <fstream>
#include <iomanip>
#include <cstdio>

void write_placement(const std::string& out_file,
                     const std::vector<Cell>& cells,
                     float_t chip_xl, float_t chip_yl,
                     float_t chip_xh, float_t chip_yh)
{
    std::ofstream fout(out_file);
    if (!fout) {
        printf("[Writer] Cannot open %s for writing\n", out_file.c_str());
        return;
    }

    fout << "UCLA pl 1.0\n";
    fout << "# Non-smooth Global Placement Result\n";
    fout << "# Chip: [" << chip_xl << ", " << chip_xh << "] x ["
         << chip_yl << ", " << chip_yh << "]\n";

    fout << std::fixed << std::setprecision(0);
    for (const auto& c : cells) {
        fout << "o" << c.id << "\t" << c.x << "\t" << c.y;
        if (c.is_terminal)
            fout << "\t: N";
        fout << "\n";
    }

    fout.close();
    printf("[Writer] Wrote %zu cells to %s\n", cells.size(), out_file.c_str());
}
