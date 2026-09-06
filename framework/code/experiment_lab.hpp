#pragma once

// V3 metrics-only global-view experiment entry point.  Kept separate from
// main.cpp so the legacy replay CLI remains reproducible, while still linking
// exactly the same exact HPWL/density core and optimizer implementation.
int run_global_view_lab(int argc, char** argv);
