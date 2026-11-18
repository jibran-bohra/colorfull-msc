# colorfull-msc

Colour-algebra calculations for my MSc thesis, *The Jet/CGC Correspondence: A
Perspective through Conformal Transformations* (UCT, 2023), redone with
[ColorFull](https://colorfull.hepforge.org/) 1.1 by M. Sjödahl.

Each program in `examples/` compares ColorFull's output with a result you can
work out by hand or one from the old Mathematica notebooks in
`msc-qcd-journal`.

| Example | What it does |
|---|---|
| [`01_basic_color_factors`](examples/01_basic_color_factors.cc) | Textbook SU(Nc) colour factors: q q̄ and q q̄ g norms, Tr(tᵃtᵇ), a q q̄ colour correlator, f^abc f^abc, and the q q̄ g g trace basis |
| [`02_one_loop_renormalization`](examples/02_one_loop_renormalization.cc) | Colour factors of the one-loop QCD self-energy and vertex diagrams, and β₀ |
| [`03_stereographic_projection`](examples/03_stereographic_projection.cc) | Stereographic map from `notes-CFT/stereographic.nb`; the BMS kernel on the sphere equals 2 × the JIMWLK dipole kernel on the plane; ColorFull colour correlators T_i·T_j for q q̄ and q g q̄ |

## Running an example

ColorFull has to be built first (see below). Then, from the repository root:

```sh
make run EX=03_stereographic_projection   # compile and run
make EX=03_stereographic_projection       # compile only (binary goes in bin/)
make clean                                # remove compiled programs
```

Programs run inside `results/`, so any files they write end up there.

## Building ColorFull

ColorFull 1.1 is GPL-2.0-only and is not included in this repository. Get
`ColorFull-1.1.tar.gz` from [colorfull.hepforge.org](https://colorfull.hepforge.org/)
(in a browser; the site blocks scripted downloads) or from the openSUSE source
RPM at download.opensuse.org/repositories/science/openSUSE_Tumbleweed/src/.
Locally it is kept in `downloads/rpm/`, which git ignores. To build it into
`install/`:

```sh
tar xzf downloads/rpm/ColorFull-1.1.tar.gz
cd ColorFull-1.1
./configure --prefix="$PWD/../install"
make && make install
cd ..
```

The project `Makefile` compiles with `clang++ -std=c++11` against
`install/include` and `install/lib/ColorFull`.

## Layout

```
examples/    one .cc file per calculation
Makefile     build/run helper
```

Ignored by git: `ColorFull-1.1/` (unpacked source), `install/`, `bin/`,
`results/`, logs, and the local-only `downloads/` (sources) and `lib/`
(thesis PDF).
