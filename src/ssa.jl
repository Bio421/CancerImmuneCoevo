# This file contains the ssa function that is used to simulate the model

using Random
using Distributions

# do nothing
default_callback(::Any...) = nothing

"""
    ssa(callback; ...) -> (t, epoch, reaction_nt, terminate_reason)

Simulate the cancer-immune coevolution model using the stochastic simulation algorithm (SSA).

## Arguments

### Callback Function

Callback function is the only positional argument of the `ssa` function,
which is executed after each epoch, used to record the state of the system.

The signature of the callback function should be:

```julia
function callback(t::Real, epoch::Int, reaction_nt::NamedTuple)
    # do something
end
```

The `t` is the current time of the simulation, `epoch` is the current epoch number,
`reaction_nt` is a named tuple containing the reaction lists of this system.

### Keyword Arguments

#### Termination Conditions

- `t_max::Real`: maximum time of the simulation, after which the simulation will be terminated, default to 16.0;
- `epoch_max::Integer`: maximum number of epochs (one epoch means one reaction event), default to 100,000,000;
- `pop_max::Integer`: maximum number of cells, after which the simulation will be terminated, default to 30,000;

#### Model Parameters

##### Reaction Coefficients

- `b_c::Real`: cancer cell birth rate, default to 1.0;
- `d_c::Real`: cancer cell death rate, default to 0.1;
- `b_e::Real`: passive recuitment rate of effector cells, default to 0.2;
- `d_e::Real`: effector cell death rate, default to 0.1;
- `α_0::Real`: activation coefficient, default to `1e-3`;
- `β_0::Real`: killing coefficient, default to `1e-3`;
- `γ_0::Real`: inhibition coefficient, default to `1e-3`;
- `h_α::Real`: activation handling time (by default, the activation rate is type II functional response), default to 1e-2;
- `fα`: activation rate function, default to `FunctionalResponse.TypeII(α_0, h_α)`, override `α_0` and `h_α`;
- `fβ`: killing rate function, default to `FunctionalResponse.TypeI(β_0)`, override `β_0`;
- `fγ`: inhibition rate function, default to `FunctionalResponse.TypeI(γ_0)`, override `γ_0`;

##### Mutation Parameters

- `λ::Real`: the mean number of mutation sites during cell division (Poisson distribution), default to 1.0;
- `p_a::Real`: the probability of a mutation site results in a new antigenic epitope, default to 0.075;
- `p_e::Real`: the probability of a mutation site results in a new immunogenic epitope, default to 1e-6;
- `m_i::Real`: the mean of the exponential distribution of immunogenicity, default to 1.0;
- `m_a::Real`: the mean of the exponential distribution of antigenicity, default to 1.0;
- `dist_i::Exponential`: immunogenicity distribution, default to `Exponential(m_i)`, override `m_i`;
- `dist_a::Exponential`: antigenicity distribution, default to `Exponential(m_a)`, override `m_a`;

#### Other Parameters

- `rng::Random.AbstractRNG`: random number generator used in the simulation, default to `Random.default_rng()`, you can specify a seed to reproduce the result;

## Return

The return value is the tuple `(t, epoch, reaction_nt, terminate_reason)`, which
is similar to signature of the callback function but contains the final state of the simulation.

The `terminate_reason` is a symbol indicating the reason of termination, which can be one of the following:
- `:time_limit`: the simulation reaches the maximum time `t_max`;
- `:epoch_limit`: the simulation reaches the maximum number of epochs `epoch_max`;
- `:extinction`: the cancer cells are extinct;
- `:population_limit`: the total number of cells exceeds the maximum number of cells `pop_max`.
"""
function ssa(
    callback = default_callback; # callback function, executed after each epoch
    t_max::Real = 16.0, # maximum time
    rng::Random.AbstractRNG = Random.default_rng(), # random number generator
    epoch_max::Integer = 100_000_000, # maximum number of epochs
    pop_max::Integer = 30_000, # maximum number of cells
    # mutation parameters
    λ::Real = 1.0, # the mean mutation site number during cell division
    p_a::Real = 0.075, # the probability of a mutation site results in a new antigenic epitope
    p_e::Real = 1e-6, # the probability of a mutation site results in a new immunogenic epitope
    m_i::Real = 1.0, # the mean of exponential distribution of immunogenicity
    m_a::Real = 1.0, # the mean of exponential distribution of antigenicity
    dist_i = nothing, # immunogenicity distribution
    dist_a = nothing, # antigenicity distribution
    # cancer cell intrinsic parameters
    b_c::Real = 1.0, # cancer cell birth rate
    d_c::Real = 0.1, # cancer cell death rate
    # effector cell intrinsic parameters
    b_e::Real = 0.2, # effector cell birth rate
    d_e::Real = 0.1, # effector cell death rate
    # immune interaction parameters
    α_0::Real = 1e-3, # activation coefficient
    β_0::Real = 1e-3, # inhibition coefficient
    γ_0::Real = 1e-3, # immune escape (cancer kills immune cell) coefficient
    h_α::Real = 1e-2, # activation handling time
    fα = nothing, # activation rate function
    fβ = nothing, # inhibition rate function
    fγ = nothing, # escape rate function
)
    # Start with underscore to avoid confliction with the variable name of arguments.
    # The confliction will cause type instability, despite the compiler can optimize it.
    # Use the underscore to avoid the confliction to make @code_warntype happy.
    _t_max, _λ, _p_a, _p_e, _m_i, _m_a, _b_c, _d_c, _b_e, _d_e, _α_0, _β_0, _γ_0, _h_α =
        promote(t_max, λ, p_a, p_e, m_i, m_a, b_c, d_c, b_e, d_e, α_0, β_0, γ_0, h_α)

    _epoch_max, _pop_max = promote(epoch_max, pop_max)

    ssa_typed(
        callback;
        rng,
        t_max = _t_max,
        epoch_max = _epoch_max,
        pop_max = _pop_max,
        λ = _λ,
        p_a = _p_a,
        p_e = _p_e,
        b_c = _b_c,
        d_c = _d_c,
        b_e = _b_e,
        d_e = _d_e,
        dist_i = isnothing(dist_i) ? Exponential(_m_i) : dist_i,
        dist_a = isnothing(dist_a) ? Exponential(_m_a) : dist_a,
        fα = isnothing(fα) ? FunctionalResponse.TypeII(_α_0, _h_α) : fα,
        fβ = isnothing(fβ) ? FunctionalResponse.TypeI(_β_0) : fβ,
        fγ = isnothing(fγ) ? FunctionalResponse.TypeI(_γ_0) : fγ,
    )
end

# type stable version of ssa
function ssa_typed(
    callback; # callback function, executed after each epoch
    rng::Random.AbstractRNG,
    t_max::F,
    epoch_max::I,
    pop_max::I,
    # mutation parameters
    λ::F,
    p_a::F,
    p_e::F,
    b_c::F,
    d_c::F,
    b_e::F,
    d_e::F,
    dist_i,
    dist_a,
    fα,
    fβ,
    fγ,
) where {F<:AbstractFloat,I<:Integer}
    #======================== type definitions ========================#
    FA = typeof(fα)
    FI = typeof(fβ)
    FE = typeof(fγ)

    M = LazyVCat{I,Vector{UnitRange{I}}}

    IDX = Tuple{I,I}
    IDCS = Vector{IDX}
    C = CancerCell{
        M,
        @NamedTuple{
            cancer_cell_division::IDX,
            cancer_cell_death::IDX,
            immune_inhibition::IDX,
        }
    }
    E = CLTPopulation{
        I,
        F,
        @NamedTuple{
            clt_migration::IDX,
            clt_death::IDX,
            immune_activation::IDX,
            immune_inhibition::IDCS,
            immune_escape::IDX,
        }
    }
    CE = Tuple{C,Vector{E}}

    CB = CancerCellDivision{F}
    CD = CancerCellDeath{F}
    EM = CLTMigration{F}
    ED = CLTDeath{F}
    IA = Activation{F,FA}
    II = ImmumeInhition{F,FI}
    IE = ImmuneEscape{F,FE}

    CBL = ReactionList{F,CB,CE}
    CDL = ReactionList{F,CD,CE}
    EML = ReactionList{F,EM,E}
    EDL = ReactionList{F,ED,E}
    IAL = ReactionList{F,IA,E}
    IIL = ReactionList{F,II,CE}
    IEL = ReactionList{F,IE,E}

    #====================== initial reaction lists ======================#
    cancer_cell_division_list = CBL()
    cancer_cell_death_list = CDL()
    clt_migration_list = EML()
    clt_death_list = EDL()
    immune_activation_list = IAL()
    immune_inhibition_list = IIL()
    immune_escape_list = IEL()

    #====================== initial reaction nodes ======================#
    # cancer cell intrinsic reactions
    cancer_cell_division = CB(b_c)
    cancer_cell_division_slot = popslot!(cancer_cell_division_list)
    cancer_cell_death = CD(d_c)
    cancer_cell_death_slot = popslot!(cancer_cell_death_list)
    # effector cell intrinsic reactions
    clt_migration = EM(b_e)
    clt_migration_slot = popslot!(clt_migration_list)
    clt_death = ED(zero(F), d_e)
    clt_death_slot = popslot!(clt_death_list)
    # immune interaction
    activation = IA(zero(F), fα)
    activation_slot = popslot!(immune_activation_list)
    immune_inhibition = II(zero(F), zero(F), fβ)
    immune_inhibition_slot = popslot!(immune_inhibition_list)
    immune_escape = IE(zero(F), fγ)
    immune_escape_slot = popslot!(immune_escape_list)

    #======================== initial cell nodes ========================#
    ID_counter = Ref(one(I)) # the id of mutation

    # a cancer cell without any mutation
    # the activation and immune escape are not CLTs related to the cancer cell
    # because the cancer cell number of a specific mutation is stored in effector cell nodes
    cancer_cell = C((
        cancer_cell_division = cancer_cell_division_slot,
        cancer_cell_death = cancer_cell_death_slot,
        immune_inhibition = immune_inhibition_slot,
    ))
    # the effector cell related to mutation 0 with zero immunogenicity and antigenicity
    # this effector cell is not related to any cancer cell just used as a template for effector cells
    # any new effector cell will be duplicated from this template
    # during our analysis, this effector cell will be excluded
    effector_cell = E(
        zero(F),
        zero(F),
        (
            clt_migration = clt_migration_slot,
            clt_death = clt_death_slot,
            immune_activation = activation_slot,
            immune_inhibition = [immune_inhibition_slot],
            immune_escape = immune_escape_slot,
        ),
    )

    #=================== insert reaction and children ===================#
    effectors = E[]
    unsafe_insert!(
        cancer_cell_division_list,
        cancer_cell_division_slot,
        cancer_cell_division,
        (cancer_cell, effectors),
    )
    unsafe_insert!(
        cancer_cell_death_list,
        cancer_cell_death_slot,
        cancer_cell_death,
        (cancer_cell, effectors),
    )
    unsafe_insert!(clt_migration_list, clt_migration_slot, clt_migration, effector_cell)
    unsafe_insert!(clt_death_list, clt_death_slot, clt_death, effector_cell)
    unsafe_insert!(immune_activation_list, activation_slot, activation, effector_cell)
    unsafe_insert!(
        immune_inhibition_list,
        immune_inhibition_slot,
        immune_inhibition,
        (cancer_cell, effectors),
    )
    unsafe_insert!(immune_escape_list, immune_escape_slot, immune_escape, effector_cell)

    #========================== reaction lists ==========================#

    reaction_tp = (
        cancer_cell_division_list,
        cancer_cell_death_list,
        clt_migration_list,
        clt_death_list,
        immune_activation_list,
        immune_inhibition_list,
        immune_escape_list,
    )

    reaction_nt = (
        cancer_cell_division = cancer_cell_division_list,
        cancer_cell_death = cancer_cell_death_list,
        clt_migration = clt_migration_list,
        clt_death = clt_death_list,
        immune_activation = immune_activation_list,
        immune_inhibition = immune_inhibition_list,
        immune_escape = immune_escape_list,
    )

    #=========================== main SSA loop ==========================#
    t = zero(F)
    epoch = zero(epoch_max)

    total_sum = sum(l -> l.sum[], reaction_tp)
    total_cancer_cell = cancer_cell_division_list.sum[] / b_c
    _pop_max = convert(F, pop_max)

    # for recording initial state
    callback(t, epoch, reaction_nt)

    terminate_reason = :time_limit

    while t < t_max
        r = rand(rng) * total_sum
        τ = -log(rand(rng)) / total_sum

        if r < cancer_cell_division_list.sum[]
            cb = sample(cancer_cell_division_list, r)
            effect!(
                cb,
                reaction_nt,
                effector_cell,
                rng,
                ID_counter,
                λ,
                p_e,
                p_a,
                dist_i,
                dist_a,
            )
        elseif (r -= cancer_cell_division_list.sum[]) < cancer_cell_death_list.sum[]
            cd = sample(cancer_cell_death_list, r)
            effect!(cd, reaction_nt)
        elseif (r -= cancer_cell_death_list.sum[]) < clt_migration_list.sum[]
            em = sample(clt_migration_list, r)
            effect!(em, reaction_nt)
        elseif (r -= clt_migration_list.sum[]) < clt_death_list.sum[]
            ed = sample(clt_death_list, r)
            effect!(ed, reaction_nt)
        elseif (r -= clt_death_list.sum[]) < immune_activation_list.sum[]
            ia = sample(immune_activation_list, r)
            effect!(ia, reaction_nt)
        elseif (r -= immune_activation_list.sum[]) < immune_inhibition_list.sum[]
            ii = sample(immune_inhibition_list, r)
            effect!(ii, reaction_nt)
        elseif (r -= immune_inhibition_list.sum[]) < immune_escape_list.sum[]
            ie = sample(immune_escape_list, r)
            effect!(ie, reaction_nt)
        else
            # This should not happen if the total sum is calculated correctly.
            # If this happens, it is a bug (or numerical error) in the code.
            # In our real simulation, this never happens.
            @error "Invalid reaction, $total_sum, $r"
            break
        end

        t += τ
        epoch += 1
        total_sum = sum(l -> l.sum[], reaction_tp)
        total_cancer_cell = cancer_cell_division_list.sum[] / b_c

        callback(t, epoch, reaction_nt)

        epoch >= epoch_max && (terminate_reason = :epoch_limit; break)
        total_cancer_cell < 0.001 && (terminate_reason = :extinction; break)
        total_cancer_cell > _pop_max && (terminate_reason = :population_limit; break)
    end

    return t, epoch, reaction_nt, terminate_reason
end
