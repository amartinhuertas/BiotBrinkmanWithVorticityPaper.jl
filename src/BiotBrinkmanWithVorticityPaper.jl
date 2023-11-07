module BiotBrinkmanWithVorticityPaper
    using Gridap
    using LinearOperators
    using LinearAlgebra
    using Krylov
    using Printf
    # using GridapPardiso
    # using SparseMatricesCSR

    include("PreconditioningTools.jl")
    include("BiotBrinkmanWithVorticity2DTools.jl")
    include("BiotBrinkmanWithVorticity3DTools.jl")

    export generate_2D_model
    export assemble_2D
    export solve_2D_riesz_mapping_preconditioner_blocks
    export compute_errors_2D
    export default_2D_u_ex, default_2D_p_ex, default_2D_v_ex

    export generate_3D_model
    export assemble_3D
    export solve_3D_riesz_mapping_preconditioner_blocks
    export compute_errors_3D
    export compute_B3_error_norms
    export default_3D_u_ex, default_3D_p_ex, default_3D_v_ex
end
