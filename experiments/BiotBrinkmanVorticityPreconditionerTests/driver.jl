module Driver
  # Add here whichever modules/packages your
  # driver function may need
  using Gridap
  using BiotBrinkmanWithVorticityPaper

  function driver(nk, order, μ, λ, ν, κ, α, c_0, prec_variant)
    model=generate_3D_model(nk)
    Ω  = Triangulation(model)
    dΩ = Measure(Ω,2*(order+2))
    Λ  = SkeletonTriangulation(model)
    Σ = BoundaryTriangulation(model,tags = "Sigma")
   
    dΛ = Measure(Λ,2*(order+2)-1)
    dΣ = Measure(Σ,2*(order+2)-1)

    h_e = CellField(get_array(∫(1)dΛ),Λ)
    h_e_Σ = CellField(get_array(∫(1)dΣ), Σ)

    # We may have here one more than one output
    # (as many as we like)
    output=Dict()
    rtol=1.0e-06

    u_ex=default_3D_u_ex
    p_ex=default_3D_p_ex
    v_ex=default_3D_v_ex

    op=assemble_3D(model, order, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

    xh,stats,rel_residual_convergence=solve_3D_riesz_mapping_preconditioner_blocks(
                                                     op, dΩ, dΛ, dΣ, h_e, h_e_Σ,
                                                     μ, λ, ν, κ, α, c_0;
                                                     verbose=true,
                                                     rtol=rtol,
                                                     prec_variant=prec_variant)
    
    error=
        compute_errors_3D(xh, dΩ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

    output["hh"]              = sqrt(3/2)/(2^(nk))
    output["||r_i||/||r_0||"] = rel_residual_convergence
    output["error"]           = error
    output["ndofs"]           = length(get_free_dof_values(xh))
    output["niters"]          = stats.niter
    output["rtol"]            = rtol
    output["prec_variant"]    = prec_variant
    output
  end
end
