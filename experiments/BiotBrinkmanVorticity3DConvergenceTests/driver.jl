module Driver
  # Add here whichever modules/packages your
  # driver function may need
  using Gridap
  #using GridapPardiso
  using BiotBrinkmanWithVorticityPaper

  function driver(nk, order, μ, λ, ν, κ, α, c_0)
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
    #ps = PardisoSolver(mtype=GridapPardiso.MTYPE_REAL_SYMMETRIC_INDEFINITE)
    #	 	       #msglvl=GridapPardiso.MSGLVL_VERBOSE)
    #solver = LinearFESolver(ps)
    #xh=solve(solver,op)
    xh=solve(op)
    
    eu,ev,eω,eφp=
       compute_B3_error_norms(xh, op, 
                              dΩ, dΛ, dΣ, h_e, h_e_Σ, 
                              μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

    output["hh"]        = sqrt(3.0/2.0)/(2^(nk))
    output["eu"]        = eu
    output["ev"]        = ev
    output["eω"]        = eω
    output["eφp"]       = eφp
    output["ndofs"]     = length(get_free_dof_values(xh))
    output
  end
end
