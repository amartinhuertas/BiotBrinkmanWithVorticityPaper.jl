module BiotBrinkmanWithVorticityPaper
    using Gridap
    using LinearOperators
    using LinearAlgebra
    using Krylov
    using Printf

    include("PreconditioningTools.jl")
    include("BiotBrinkmanWithVorticityTools.jl")

    export generate_model3d
    export assemble_biotbrinkman
    export compute_errors_biotbrinkman
    export default_u_ex, default_p_ex, default_v_ex
end
