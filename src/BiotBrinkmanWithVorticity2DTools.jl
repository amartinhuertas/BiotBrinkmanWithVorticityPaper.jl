# manufactured solutions
default_2D_u_ex(x)  = VectorValue(sin(π*(x[1]+x[2])), cos(π*(x[1]^2+x[2]^2)))
default_2D_p_ex(x)  = sin(π*x[1]+x[2])*sin(π*x[2])
default_2D_v_ex(x)  = VectorValue(sin(π*x[1])*sin(π*x[2]),cos(π*x[1])*cos(2*π*x[2]))  

function build_2D_analytical_functions(μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
  Id = TensorValue(1.0,0.0,0.0,1.0)
  φ_ex(x) = α*p_ex(x)-λ*(∇⋅u_ex)(x)
  ω_ex(x) = ν/sqrt(κ)*(∇×v_ex)(x)
  σ_ex(x) = 2*μ*ε(u_ex)(x)- φ_ex(x)*Id
  b_ex(x) = -(∇⋅σ_ex)(x)
  f_ex(x) = ν/κ*v_ex(x) + ν/sqrt(κ)*(∇×ω_ex)(x)-ν^2/κ*∇(∇⋅v_ex)(x)+∇(p_ex)(x)
  g_ex(x) = -(c_0 + α^2/λ)*p_ex(x) + α/λ*φ_ex(x) - (∇⋅v_ex)(x)
  return φ_ex, ω_ex, σ_ex, b_ex, f_ex, g_ex
end

function generate_2D_model(nk)
    domain =(0,1,0,1)
    n      = 2^nk + 1
    # Discrete model
    partition = (n,n)
    CartesianDiscreteModel(domain, partition) |> simplexify
end

function Gridap.Fields.grad2curl(∇u::VectorValue{2})
    VectorValue(∇u[2], -∇u[1])
  end 

function assemble_2D(model, k, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
  φ_ex, ω_ex, σ_ex, b_ex, f_ex, g_ex = 
    build_2D_analytical_functions(μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

   # Define reference FE (P2/P0/P1/P0/P0)
   reffe_u = ReferenceFE(lagrangian,VectorValue{2,Float64},k+2)
   reffe_v = ReferenceFE(raviart_thomas,Float64,k)
   reffe_ω = ReferenceFE(lagrangian,Float64,k+1) # ND0 becomes P1 in 2D
   reffe_φ = ReferenceFE(lagrangian,Float64,k)
   reffe_p = ReferenceFE(lagrangian,Float64,k)

   # FESpaces
   Uh_ = TestFESpace(model,reffe_u,dirichlet_tags="boundary",conformity=:H1)
   Vh_ = TestFESpace(model,reffe_v,dirichlet_tags="boundary",conformity=:HDiv)
   Wh_ = TestFESpace(model,reffe_ω,dirichlet_tags="boundary",conformity=:H1)
   Zh_ = TestFESpace(model,reffe_φ,conformity=:L2)
   Qh_ = TestFESpace(model,reffe_p,conformity=:L2)
   Lh_ = ConstantFESpace(model)

   Uh = TrialFESpace(Uh_,u_ex)
   Vh = TrialFESpace(Vh_,v_ex)
   Wh = TrialFESpace(Wh_,ω_ex)
   Zh = TrialFESpace(Zh_)
   Qh = TrialFESpace(Qh_)
   Lh = TrialFESpace(Lh_)

   Yh = MultiFieldFESpace([Uh_,Vh_,Wh_,Zh_,Qh_,Lh_,Lh_])
   Xh = MultiFieldFESpace([Uh,Vh,Wh,Zh,Qh,Lh,Lh])

   Ω = Triangulation(model)
   dΩ = Measure(Ω,2*(k+4))

   # Bilinear form
   lhs((u,v,ω,φ,p,ρ1,ρ2),(γ,ζ,θ,ψ,q,ξ1,ξ2)) =  ∫(2*μ*ε(γ) ⊙ ε(u) -
                                          (∇⋅γ)*φ +
                                          ν/κ*v⋅ζ +
                                          ν/sqrt(κ)*(∇×ω)⋅ζ +
                                          ν^2/κ*(∇⋅v)*(∇⋅ζ) -
                                          (∇⋅ζ)*p - ω*θ +
                                          ν/sqrt(κ)*(∇×θ)⋅v -
                                          1/λ*φ*ψ +
                                          α/λ*ψ*p -
                                          (∇⋅u)*ψ -
                                          (c_0+α^2/λ)*p*q +
                                          α/λ*φ*q -
                                          (∇⋅v)*q +
                                          φ*ξ1 +
                                          p*ξ2 +
                                          ρ1*ψ +
                                          ρ2*q)dΩ 
                
   # Linear functional
   rhs((γ,ζ,θ,ψ,q,ξ1,ξ2)) = ∫(b_ex⋅γ + f_ex⋅ζ + g_ex*q + φ_ex*ξ1 + p_ex*ξ2)dΩ
                      
   # Build affine FE operator 
   op = AffineFEOperator(lhs,rhs,Xh,Yh)
end

function compute_errors_2D(xh, dΩ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
  uh, vh, ωh, φh, ph, ρ1h, ρ2h = xh

  Xh = Gridap.FESpaces.get_fe_space(xh)
  _, _, _, Zh, _, _, _ = Xh

  φ_ex, ω_ex, σ_ex, b_ex, f_ex, g_ex = 
    build_2D_analytical_functions(μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

  # project the mass loss onto Zh, and take the l_infty norm of that quantity...
  a(u,v) = ∫( u*v )*dΩ
  l(v) = ∫(v*(-(c_0+α^2/λ)*ph+α/λ*φh-(∇⋅vh)-g_ex))*dΩ
  projection = AffineFEOperator(a,l,Zh,Zh)
  loss = solve(projection)
  lossh=norm(get_free_dof_values(loss),Inf)

  # compute errors in the non-weighted norms
  error_u = sqrt(sum(∫((u_ex-uh)⋅(u_ex-uh))*dΩ +∫(∇(u_ex-uh)⊙∇(u_ex - uh))*dΩ))
  error_v = sqrt(sum(∫((v_ex-vh)⋅(v_ex-vh))*dΩ + ∫((∇⋅(v_ex-vh))*(∇⋅(v_ex-vh)))*dΩ))
  error_ω = sqrt(sum(∫((ω_ex-ωh)*(ω_ex-ωh))*dΩ + ∫((∇×(ω_ex-ωh))⋅(∇×(ω_ex-ωh)))*dΩ))
  error_φ = sqrt(sum(∫((φ_ex-φh)*(φ_ex-φh))*dΩ))
  error_p = sqrt(sum(∫((p_ex-ph)*(p_ex-ph))*dΩ))
  error_u,error_v,error_ω,error_φ,error_p, lossh
end
