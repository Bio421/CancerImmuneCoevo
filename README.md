# CancerImmuneCoevo

Implementation of SSA (Gillespie algorithm) for the Cancer-Immune Coevolution model
in our paper "Cancer-immune coevolution dictated by antigenic mutation accumulation".

## Installation

This package is not registered in the official Julia package registry,
so you need to install it from the GitHub repository or clone the repository.

To install the package, you can install with julia built-in package manager:

```julia
using Pkg;Pkg.add(url="https://github.com/Bio421/CancerImmuneCoevo.git")
```

or within `pkg` mode of REPL:

```julia
pkg> add https://github.com/Bio421/CancerImmuneCoevo.git
```

**NOTE**: Julia 1.10 or later is recommended, lower versions may not work.

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

## Implementation

The biggest challenges in implementing the model is the large number of reactions
(about 10^5 - 10^6), which means we need to calculate the propensity of each reaction
after each reaction. This is very time-consuming and inefficient.

To solve this problem, we can avoid recalculating the propensity of reactions
that are not affected by the last reaction. This can be done by maintaining a
dependency graph of reactions. However, the dependency graph dynamically changes
during the evolution of the system, which makes it difficult to maintain.
To maintain the dependency graph, we maintain a list of substances that
are affected by the reaction, and for each substance, we maintain a list of reactions
that depend on it. This way, we can quickly find the reactions that need to be updated
after each reaction.

## License

This package is released under the MIT License, see `LICENSE.md` for more information.
