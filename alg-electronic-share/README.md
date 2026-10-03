# NSP Electronic Placement

This package contains the source code, build configuration, focused tests, and
optional visualisation code for the `alg-electronic` placer. Benchmark datasets,
compiled binaries, experiment logs, placement outputs, and exported snapshots
are intentionally excluded.

## Directory layout

```text
.
├── src/                 C++ placer and snapshot-export implementation
├── tests/               Focused C++ tests
├── Makefile             Build and test commands
├── visualize_nsp.py     Offline GIF renderer for exported snapshots
├── ELECTROSTATIC.md     Electrostatic density-model notes
└── ispd2005/            Dataset location (create locally; not included)
```

## Requirements

- C++17 compiler with OpenMP support. The supplied Windows Makefile defaults to
  `D:/MingGW/ucrt64/bin/g++.exe`; override `MINGW_HOME` if MinGW is elsewhere.
- GNU Make (`mingw32-make` on a typical Windows MinGW installation).
- Python 3.10+ with `numpy`, `matplotlib`, and `pillow` for GIF rendering:

```powershell
python -m pip install numpy matplotlib pillow
```

## Dataset placement

Obtain the complete ISPD 2005 Bookshelf benchmark separately. Create an
`ispd2005` directory at this repository root and preserve the following layout:

```text
ispd2005/
├── adaptec1/
│   ├── adaptec1.aux
│   ├── adaptec1.nodes
│   ├── adaptec1.nets
│   ├── adaptec1.pl
│   ├── adaptec1.scl
│   ├── adaptec1.wts
│   └── adaptec1.eplace-ip.pl
├── adaptec2/
├── adaptec3/
├── adaptec4/
├── bigblue1/
├── bigblue2/
├── bigblue3/
└── bigblue4/
```

Run the commands below from the repository root. The first argument is the
Bookshelf file base path, without any extension. For example,
`ispd2005/adaptec1/adaptec1` resolves the files shown above. The default
initial-placement mode uses `<benchmark>.eplace-ip.pl`; retain that file for
the default `--init eplace` flow, or choose another `--init` mode.

## Build

Windows (MinGW):

```powershell
mingw32-make -j4
```

If needed, specify the MinGW root explicitly:

```powershell
mingw32-make MINGW_HOME=D:/path/to/mingw64 -j4
```

macOS (Homebrew LLVM/OpenMP path expected by the Makefile):

```bash
make -j4
```

The executable is `nsp_placer.exe` on Windows and `nsp_placer` elsewhere.

## Run

Run one benchmark with the default NSP configuration:

```powershell
.\nsp_placer.exe ispd2005/adaptec1/adaptec1
```

Useful optional flags include `--init eplace|uniform|random|row|scaled|gaussian|quadrant`,
`--cwtr`, `--spread`, `--diffuse`, `--max-iters N`, and `--adam`. Run the
executable without arguments to see the complete current option list.

Run all eight ISPD 2005 benchmarks in PowerShell:

```powershell
$benches = 'adaptec1','adaptec2','adaptec3','adaptec4','bigblue1','bigblue2','bigblue3','bigblue4'
foreach ($bench in $benches) {
  .\nsp_placer.exe "ispd2005/$bench/$bench"
}
```

## Visualisation

The C++ program only writes snapshots when `--visualize` is supplied. Export a
snapshot every 10 iterations into a dedicated generated-output directory:

```powershell
.\nsp_placer.exe ispd2005/adaptec1/adaptec1 `
  --visualize visualizations/adaptec1 `
  --visualize-every 10
```

The directory will contain `metadata.txt`, `metrics.csv`, and
`snapshot_*.bin`. Render these files offline as a GIF:

```powershell
python .\visualize_nsp.py `
  --snapshots visualizations/adaptec1 `
  --benchmark ispd2005/adaptec1/adaptec1 `
  --output visualizations/adaptec1/optimization_animation.gif `
  --fps 5
```

`visualize_nsp.py` reads the original `<benchmark>.nodes` file to obtain cell
dimensions and fixed objects, so the dataset must remain at the path passed to
`--benchmark` while rendering. Use `--max-points N` to cap the number of
movable cells shown in dense benchmarks.

## Tests

```powershell
mingw32-make test-electric
mingw32-make test-lambda
mingw32-make test-trajectory
```

The generic `mingw32-make test` also requires `ispd2005/adaptec1` in the
directory layout above.
