module CancerImmuneCoevo

export ssa
export children, mutations, antigenic, passenger, antigenicity, immunogenicity, e_num, c_num

include("utils.jl")

# Concrete type definitions
include("cells.jl")

include("reactionlist.jl")

include("reactions.jl")

# Functional response functions
include("response.jl")

# Reaction effects
include("effects.jl")

# main ssa function
include("ssa.jl")

end # module CancerImmune
