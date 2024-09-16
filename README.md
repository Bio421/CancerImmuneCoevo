# CancerImmuneCoevo

Implementation of SSA (Gillespie algorithm) for the Cancer-Immune Coevolution model
in our paper "Cancer-immune coevolution dictated by antigenic mutation accumulation".

## Requirements

This package is implemented in Julia, a high-performance dynamic programming language.
So you need to have Julia installed on your system to run the simulation,
see [Julia's download page](https://julialang.org/downloads/) for more information.

**NOTE**: Julia 1.10 or later is recommended, lower versions may not work.

All the dependencies are listed in the `Project.toml` file and will be installed
automatically when you instantiate the package.

## Installation

You can download the zip file and extract it to a directory.
Then activate and instantiate the package in the Julia REPL:

```julia
using Pkg; Pkg.activate("path/to/CancerImmuneCoevo"); Pkg.instantiate()
```

Alternatively, you can install this package directly with Julia's package manager:

```julia
using Pkg; Pkg.add(url="https://github.com/Bio421/CancerImmuneCoevo.git")
```

## Usage

The entry point of SSA is function `ssa`, you can use it to run the simulation:

```julia
using CancerImmuneCoevo; ssa();
```

The return value of `ssa` is quite large. Be careful when running it in REPL.

The augments and how to access the system state are described in the
documentation of `ssa`. You can access the documentation in the REPL:

```julia
help?> ssa
```

Here is an example of how to run the simulation and record mutation information and
effector (CTL) cell number:

```julia
using CancerImmuneCoevo

# Wrap the simulation in a let block to avoid to use global variables
let
    ts = Float64[]
    antigenic_mutations = Vector{Int}[] # record the number of antigenic mutations of each cancer cell
    immunogenic_mutations = Vector{Int}[] # record the number of immunogenic mutations of each cancer cell
    num_effector_cells = Vector{Int}[] # record the number of effector cells of each type
    antigenicity = Vector{Float64}[] # record the antigenicity of each cancer cell
    immunogenicity = Vector{Float64}[] # record the immunogenicity of each cancer cell

    δt = 0.1 # record the state every 0.1 time unit
    last_t = Ref(-Inf) # record the last time the state was recorded, box the value to avoid type instability
    # parameters are passed as keyword arguments
    # the system state can be accessed through passed a function to ssa (with do syntax)
    # when the ssa finishes, some information can be returned
    t, epoch, reaction, reason = ssa(; α_0 = 1e-3, β_0 = 1e-3) do t, epoch, reaction_nt
        t - last_t[] < δt && return

        last_t[] = t
        push!(ts, t)

        antigenic_mutations_epoch = Int[]
        immunogenic_mutations_epoch = Int[]
        # each cancer cell division reaction corresponds to an unique cancer cell
        foreach(reaction_nt.cancer_cell_division) do reaction
            cancer = children(reaction)[1]
            # record the number of antigenic and immunogenic mutations
            # without length to get detailed mutation information
            push!(antigenic_mutations_epoch, length(antigenic(cancer)))
            push!(immunogenic_mutations_epoch, length(immunogenic(cancer)))
        end
        push!(antigenic_mutations, antigenic_mutations_epoch)
        push!(immunogenic_mutations, immunogenic_mutations_epoch)

        num_effector_cells_epoch = Int[]
        antigenicity_epoch = Float64[]
        immunogenicity_epoch = Float64[]
        # each effector cell migration reaction corresponds to a type of effector cell
        foreach(reaction_nt.clt_migration) do reaction
            effector = children(reaction)
            # ignore the effector cell if it is not antigenic
            # which is an template of the effector cell type, has no effect
            iszero(antigenicity(effector)) && return

            # record the number of effector cells of each type
            push!(num_effector_cells_epoch, e_num(effector))
            # record the antigenicity and immunogenicity
            push!(antigenicity_epoch, antigenicity(effector))
            push!(immunogenicity_epoch, immunogenicity(effector))
            return
        end
        push!(num_effector_cells, num_effector_cells_epoch)

        return
    end

    # analyze the data or save it to a file (e.g. julia serialization, JLD2, etc.)
end
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

## Directory Structure

```
.
├── src
│   ├── CancerImmuneCoevo.jl # Main module
│   ├── cells.jl             # Cell types
│   ├── effects.jl           # Effects when a reaction happens
│   ├── reactionlist.jl      # A data structure to store reactions
│   ├── reactions.jl         # Reactions
│   ├── response.jl          # Functional response functions
│   ├── ssa.jl               # main function ssa
│   └── utils.jl             # Utility functions
├── LICENSE             # License
├── Manifest.toml       # Dependencies
├── Project.toml        # Dependencies
└── README.md           # This file
```

## License

This package is released under the MIT License, see `LICENSE.md` for more information.
