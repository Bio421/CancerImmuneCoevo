# CancerImmuneCoevo

Implementation of SSA (Gillespie algorithm) for the Cancer-Immune Coevolution model
in our paper "Cancer-immune coevolution dictated by antigenic mutation accumulation".

## Installation

To install the package, you can install with julia built-in package manager:

```julia
using Pkg; Pkg.add(url="https://github.com/Bio421/CancerImmuneCoevo.git")
```

or with in `pkg` mode of REPL:

```julia
pkg> add https://github.com/Bio421/CancerImmuneCoevo.git
```

## Usage

The entry point of SSA is function `ssa`, you can use it to run the simulation:

```julia
using CancerImmuneCoevo; ssa();
```

The return value of `ssa` is quite large. Be careful when running it in REPL.

The augments and how to access the system state are described in the
docstring of `ssa`, which can be accessed by:

```julia
help?> ssa
```

## License

This package is released under the MIT License, see LICENSE.md for more information.
