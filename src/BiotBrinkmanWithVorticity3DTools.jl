# default manufactured solutions
default_3D_u_ex(x)  = VectorValue(0.1*sin(π*(x[1]+x[2]+x[3])),
                                  0.1*cos(π*(x[1]^2+x[2]^2+x[3]^2)),
                                  0.1*sin(π*(x[1]+x[2]+x[3]))*cos(π*(x[1]+x[2]+x[3])))

default_3D_p_ex(x)  = sin(π*x[1])*cos(π*x[2])*sin(π*x[3])

default_3D_v_ex(x)  = VectorValue((sin(π*x[1]))^2*sin(π*x[2])*sin(2*π*x[3]), 
                                   sin(π*x[1])*(sin(π*x[2]))^2*sin(2*π*x[3]),
                                   sin(π*x[1])*sin(2*π*x[2])*sin(π*x[3]))

function build_3D_analytical_functions(μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
  Id = TensorValue(1,0,0,0,1,0,0,0,1)
  φ_ex(x) = α*p_ex(x)-λ*(∇⋅u_ex)(x)
  ω_ex(x) = sqrt(ν/κ)*(∇×v_ex)(x)
  σ_ex(x) = 2*μ*ε(u_ex)(x)-φ_ex(x)*Id
  b_ex(x) = -(∇⋅σ_ex)(x)
  f_ex(x) = 1.0/κ*v_ex(x) + sqrt(ν/κ)*(∇×ω_ex)(x)-ν/κ*∇(∇⋅v_ex)(x)+∇(p_ex)(x)
  g_ex(x) = -(c_0 + α^2/λ)*p_ex(x) + α/λ*φ_ex(x) - (∇⋅v_ex)(x)
  return φ_ex, ω_ex, σ_ex, b_ex, f_ex, g_ex
end 

function setup_3D_model_labels!(model)
    labels = get_face_labeling(model)
    add_tag!(labels,"Gamma",[1,2,4,5,6,8,21,23,26,9,10,11,13,14,16,17,18,20])
    add_tag!(labels,"Sigma",[22,24,25])
end

function generate_3D_model(nk)
  domain =(0,1,0,1,0,1)
  n      = 2^nk+1
  partition = (n,n,n)
  model=CartesianDiscreteModel(domain, partition) |> simplexify
  setup_3D_model_labels!(model)
  model
end

function assemble_3D(model, k, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
   φ_ex, ω_ex, σ_ex, b_ex, f_ex, g_ex =
     build_3D_analytical_functions(μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

   Dc = num_cell_dims(model)

   # Define reference FE (P2/RT0/ND0/P1/P0)
   reffe_u = ReferenceFE(lagrangian,VectorValue{Dc,Float64},k+2)
   reffe_v = ReferenceFE(raviart_thomas,Float64,k)
   reffe_ω = ReferenceFE(nedelec,Float64,k)
   reffe_φ = ReferenceFE(lagrangian,Float64,k+1)
   reffe_p = ReferenceFE(lagrangian,Float64,k)

   # FESpaces
   Uh_ = TestFESpace(model,reffe_u,dirichlet_tags="Gamma",conformity=:H1)
   Vh_ = TestFESpace(model,reffe_v,dirichlet_tags="Gamma",conformity=:HDiv)
   Wh_ = TestFESpace(model,reffe_ω,dirichlet_tags="Gamma",conformity=:HCurl)
   Zh_ = TestFESpace(model,reffe_φ,conformity=:H1)
   Qh_ = TestFESpace(model,reffe_p,conformity=:L2)

   Uh = TrialFESpace(Uh_,u_ex)
   Vh = TrialFESpace(Vh_,v_ex)
   Wh = TrialFESpace(Wh_,ω_ex)
   Zh = TrialFESpace(Zh_)
   Qh = TrialFESpace(Qh_)

   Yh = MultiFieldFESpace([Uh_,Vh_,Wh_,Zh_,Qh_])
   Xh = MultiFieldFESpace([Uh,Vh,Wh,Zh,Qh])

   # Set up triangulations and measures 
   # corresponding to different integration domains 
   Ω = Triangulation(model); dΩ = Measure(Ω,2*(k+2))
   Σ = BoundaryTriangulation(model,tags = "Sigma"); dΣ = Measure(Σ,2*(k+2)-1)
   
   # Obtained unit outer normal vector to the Neumann boundary
   n_Σ = get_normal_vector(Σ)


   # bilinear form
   lhs((u,v,ω,φ,p),(γ,ζ,θ,ψ,q)) = ∫(2*μ*ε(γ) ⊙ ε(u) - (∇⋅γ)*φ +
                                      1.0/κ*v⋅ζ + sqrt(ν/κ)*(∇×ω)⋅ζ + ν/κ*(∇⋅v)*(∇⋅ζ) - (∇⋅ζ)*p - 
                                      ω⋅θ + sqrt(ν/κ)*(∇×θ)⋅v -
                                      1/λ*φ*ψ + α/λ*ψ*p - (∇⋅u)*ψ -
                                      (c_0+α^2/λ)*p*q + α/λ*φ*q - (∇⋅v)*q)dΩ

   # linear functional
   rhs((γ,ζ,θ,ψ,q)) = ∫(b_ex⋅γ + f_ex⋅ζ + g_ex*q)dΩ +
                        ∫(γ⋅(σ_ex⋅n_Σ))dΣ +
                        ∫(ν/κ*(ζ⋅n_Σ)*(∇⋅v_ex))dΣ -
                        ∫(p_ex*(ζ⋅n_Σ))dΣ -
                        ∫(sqrt(ν/κ)*(v_ex⋅(θ×n_Σ)))dΣ 

   # Build affine FE operator
   op = AffineFEOperator(lhs,rhs,Xh,Yh)
end

function assemble_3D_riesz_mapping_preconditioner_blocks(op, dΩ, dΛ, dΣ,
                                                         h_e, h_e_Σ,
                                                         μ, λ, ν, κ, α, c_0;
                                                         prec_variant=:B1)
  Y1,Y2,Y3,Y4,Y5=op.test
  X1,X2,X3,X4,X5=op.trial
  @assert prec_variant in (:B1,:B2,:B3)
  if (prec_variant==:B1)
    a11_B1(u,γ)=∫(2.0*μ*ε(γ)⊙ε(u))dΩ
    a22_B1(v,ζ)=∫((1.0/κ)*v⋅ζ+1.0/κ*(∇⋅v)*(∇⋅ζ))dΩ
    a33_B1(ω,θ)=∫((ω⋅θ)+ν*(∇×ω)⋅(∇×θ))dΩ
    a44_B1(φ,ψ)=∫((1.0/λ+0.5/μ)*φ*ψ)dΩ
    a55_B1(p,q)=∫((c_0+α^2/λ+κ)*p*q)dΩ
    A11,A22,A33,A44,A55=assemble_matrix(a11_B1,X1,Y1),assemble_matrix(a22_B1,X2,Y2),
                        assemble_matrix(a33_B1,X3,Y3),assemble_matrix(a44_B1,X4,Y4),
                        assemble_matrix(a55_B1,X5,Y5)
    return A11,A22,A33,A44,A55
  elseif (prec_variant==:B2)
    a11_B2(u,γ)=∫(2.0*μ*ε(γ)⊙ε(u))dΩ
    a22_B2(v,ζ)=∫((1.0/κ)*v⋅ζ+ν/κ*(∇⋅v)*(∇⋅ζ))dΩ
    a33_B2(ω,θ)=∫((ω⋅θ)+sqrt(ν/κ)*(∇×ω)⋅(∇×θ))dΩ
    a44_B2(φ,ψ)=∫((1.0/λ+0.5/μ)*φ*ψ)dΩ
    a55_B2(p,q)=∫((c_0+α^2/λ)*p*q)dΩ + ∫((κ/ν)*(∇(p))⋅∇(q))dΩ + 
                ∫((κ/ν)*1.0/h_e*jump(p)*jump(q))dΛ + ∫((κ/ν)*1.0/h_e_Σ*p*q)dΣ
    A11,A22,A33,A44,A55=assemble_matrix(a11_B2,X1,Y1),assemble_matrix(a22_B2,X2,Y2),
                        assemble_matrix(a33_B2,X3,Y3),assemble_matrix(a44_B2,X4,Y4),
                        assemble_matrix(a55_B2,X5,Y5)
    return A11,A22,A33,A44,A55       
  elseif (prec_variant==:B3)
    a11_B3(u,γ)=∫(2.0*μ*ε(γ)⊙ε(u))dΩ
    a22_B3(v,ζ)=∫((1.0/κ)*v⋅ζ+(1.0+ν/κ)*(∇⋅v)*(∇⋅ζ))dΩ
    a33_B3(ω,θ)=∫((ω⋅θ)+ν*(∇×ω)⋅(∇×θ))dΩ
    A11,A22,A33=assemble_matrix(a11_B3,X1,Y1),assemble_matrix(a22_B3,X2,Y2),assemble_matrix(a33_B3,X3,Y3)
    Y45 = MultiFieldFESpace([Y4, Y5])
    X45 = MultiFieldFESpace([X4, X5])
    a44_B3a((φ,p),(ψ,q))= ∫((0.5/μ)*φ*ψ)dΩ + ∫((1.0/λ)*(φ+α*p)*(ψ+α*q))dΩ + ∫((1.0+c_0)*p*q)dΩ
    a44_B3b((φ,p),(ψ,q))= ∫((0.5/μ)*φ*ψ)dΩ + ∫((1.0/λ)*(φ+α*p)*(ψ+α*q))dΩ + ∫((c_0)*p*q)dΩ + 
                         ∫(κ*(∇(p))⋅∇(q))dΩ + ∫(κ/(h_e)*jump(p)*jump(q))dΛ + ∫(κ/(h_e_Σ)*p*q)dΣ

    A44a=assemble_matrix(a44_B3a,X45,Y45)
    A44b=assemble_matrix(a44_B3b,X45,Y45)
    return A11,A22,A33,A44a,A44b
  end
end 

function compute_operator_norm(B,uh)
  u=get_free_dof_values(uh)
  sqrt(dot(u,B*u))
end

function compute_B3_error_norms(xh, op, dΩ, dΛ, dΣ, h_e, h_e_Σ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
  blocks=assemble_3D_riesz_mapping_preconditioner_blocks(op, dΩ, dΛ, dΣ,
                                                         h_e, h_e_Σ,
                                                         μ, λ, ν, κ, α, c_0;
                                                         prec_variant=:B3)

  A11,A22,A33,A44a,A44b=blocks
  A11 = LinearOperator(GridapLinearSolverPreconditioner(A11))
  A22 = LinearOperator(GridapLinearSolverPreconditioner(A22))
  A33 = LinearOperator(GridapLinearSolverPreconditioner(A33))
  A4455 = LinearOperator(GridapLinearSolverPreconditioner(A44a))+
          LinearOperator(GridapLinearSolverPreconditioner(A44b))
  P=BlockDiagonalOperator(A11,A22,A33,A4455)

  uh, vh, ωh, φh, ph = xh
  φ_ex, ω_ex, _, _, _, _ =
    build_3D_analytical_functions(μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

  eu = u_ex-uh
  ev = v_ex-vh
  eω = ω_ex-ωh
  eφ = φ_ex-φh
  ep = p_ex-ph

  X1,X2,X3,X4,X5=op.trial
  X45 = MultiFieldFESpace([X4, X5])
  euh=interpolate(eu,X1)
  evh=interpolate(ev,X2)
  eωh=interpolate(eω,X3)
  eφph=interpolate([eφ,ep],X45)
  etotal=interpolate([eu,ev,eω,eφ,ep],op.trial)

  return compute_operator_norm(A11,euh),
           compute_operator_norm(A22,evh),
              compute_operator_norm(A33,eωh),
                 compute_operator_norm(A4455,eφph),
                    compute_operator_norm(P,etotal)
end 

function solve_3D_riesz_mapping_preconditioner_blocks(op, dΩ, dΛ, dΣ, h_e, h_e_Σ,
                                                      μ, λ, ν, κ, α, c_0;
                                                      rtol=1.0e-6,
                                                      itmax=500,
                                                      prec_variant=:B1,
                                                      verbose=false,
                                                      diagnostics=false)

  if (verbose)
    println("###############")
    println("μ=$(μ)")
    println("λ=$(λ)")
    println("ν=$(ν)")
    println("κ=$(κ)")
    println("α=$(α)")
    println("c_0=$(c_0)")
    println("prec_variant=$(prec_variant)")
    println("###############")
  end

  blocks=assemble_3D_riesz_mapping_preconditioner_blocks(op, dΩ, dΛ, dΣ,
                                                            h_e, h_e_Σ,
                                                            μ, λ, ν, κ, α, c_0;
                                                            prec_variant=prec_variant)

  function print_condition_on_screen(A,P)
    Adense=Array(A)
    evals=eigvals(P*Adense)
    @printf("λ_max: %16.7e λ_min: %16.7e\n",maximum(broadcast(abs,evals)),minimum(broadcast(abs,evals)))
    @printf("||A-A'||: %16.7e\n",norm(A-transpose(A)))
    @printf("condition number %16.7e\n", maximum(broadcast(abs,evals))/minimum(broadcast(abs,evals)))
    eigs = sort(broadcast(abs,evals))
    println("Five smallest eigenvalues: ", eigs[1:10])
    println("Five largest eigenvalues: ", eigs[length(eigs)-10:length(eigs)])
    eigs_new = eigs[2:length(eigs)]
    @printf("Condition number without 1 smallest %16.7e\n", maximum(eigs_new)/minimum(eigs_new))
    eigs_new = eigs_new[2:length(eigs_new)]
    @printf("Condition number without 2 smallest %16.7e\n", maximum(eigs_new)/minimum(eigs_new))
    eigs_new = eigs_new[3:length(eigs_new)]
    @printf("Condition number without 3 smallest %16.7e\n", maximum(eigs_new)/minimum(eigs_new))
  end

  if (prec_variant==:B3)
    A11,A22,A33,A44a,A44b=blocks
    A11 = LinearOperator(GridapLinearSolverPreconditioner(A11))
    A22 = LinearOperator(GridapLinearSolverPreconditioner(A22))
    A33 = LinearOperator(GridapLinearSolverPreconditioner(A33))
    A4455 = LinearOperator(GridapLinearSolverPreconditioner(A44a))+
          LinearOperator(GridapLinearSolverPreconditioner(A44b))
    P=BlockDiagonalOperator(A11,A22,A33,A4455)
    if (diagnostics)
      N  = num_free_dofs(op.trial)
      Ns = [num_free_dofs(U) for U in op.trial] 
      Np = zeros(Int,length(Ns)+1)
      Np[1]=1
      for i=2:length(Ns)+1
        Np[i]=Np[i-1]+Ns[i-1]
      end
      Pd=zeros(N,N); 
      A11inv=inv(Array(blocks[1])); 
      A22inv=inv(Array(blocks[2]));  
      A33inv=inv(Array(blocks[3])); 
      A44ainv=inv(Array(blocks[4])); 
      A44binv=inv(Array(blocks[5])); 
      range=Np[1]:Np[2]-1
      Pd[range,range]=A11inv; 
      range=Np[2]:Np[3]-1
      Pd[range,range]=A22inv;  
      range=Np[3]:Np[4]-1
      Pd[range,range]=A33inv;  
      range=Np[4]:Np[6]-1
      Pd[range,range]=A44ainv+A44binv;
      @assert norm(P*ones(N)-Pd*ones(N))/(norm(P*ones(N))) < 1.0e-10
      print_condition_on_screen(op.op.matrix,Pd)
    end
  else 
    P=BlockDiagonalOperator(LinearOperator.(GridapLinearSolverPreconditioner.(blocks))...)
  end  
  
  A=op.op.matrix
  b=op.op.vector
  r=copy(b)
  conv_history = Float64[]
  function custom_stopping_condition(solver::KrylovSolver, A, b, r, tol)
    mul!(r, A, solver.x)
    r .-= b                       # r := b - Ax
    bool = norm(r) ≤ tol*norm(b)  # tolerance based on the 2-norm of the residual
    if (verbose)
        @printf("||b-Ax||/||b||: %16.7e\n",norm(r)/norm(b))
    end
    push!(conv_history,norm(r)/norm(b))
    return bool
  end
  minres_callback(solver) = custom_stopping_condition(solver, A, b, r, rtol)
  (xdofs, stats) = minres(A, b, M=P; itmax=itmax, verbose=1,
                        callback=minres_callback,
                        atol=0.0,
                        rtol=0.0,
                        etol=0.0)
  xh=FEFunction(op.trial,xdofs)

  if (verbose)
    println("#iters: $(stats.niter)")
  end
  xh, stats, conv_history
end

function compute_errors_3D(xh, dΩ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
  uh, vh, ωh, φh, ph = xh

  φ_ex, ω_ex, _, _, _, _ =
    build_3D_analytical_functions(μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

  eu = u_ex-uh
  ev = v_ex-vh
  eω = ω_ex-ωh
  eφ = φ_ex-φh
  ep = p_ex-ph

  # error in the weighted norm
  error = sum(∫(2.0*μ*ε(eu)⊙ε(eu) + 
                1.0/κ*ev⋅ev + 
                ν/κ*(divergence(ev)*divergence(ev)) +
                eω⋅eω +
                ν*(curl(eω)⋅curl(eω)) + 
                0.5/μ*eφ*eφ + 
                (c_0+κ/ν)*ep*ep + 
                1.0/λ*(eφ+α*ep)*(eφ+α*ep))dΩ)
    
  error
end
