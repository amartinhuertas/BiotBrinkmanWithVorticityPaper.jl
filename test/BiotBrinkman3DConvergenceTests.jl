using BiotBrinkmanWithVorticityPaper
using Gridap
using Gridap.FESpaces
using Printf
#using GridapPardiso
#using SparseMatricesCSR 
using Test

function convergence_test_3D(order, μ, λ, ν, κ, α, c_0, nkmax; 
                             u_ex=default_3D_u_ex, 
                             p_ex=default_3D_p_ex, 
                             v_ex=default_3D_v_ex)

    eu   = Float64[]
    ev   = Float64[]
    eω   = Float64[]
    eφp  = Float64[]
    etot = Float64[]
    
    hh   = Float64[]
    ru   = Float64[]
    rv   = Float64[]
    rω   = Float64[]
    rφp  = Float64[]
    rtot = Float64[]
    nn = Int[]

    push!(ru,0.)
    push!(rv,0.)
    push!(rω,0.)
    push!(rφp,0.)
    push!(rtot,0.)


    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_3D_model(nk)
        op=assemble_3D(model, order, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

        #ps = PardisoSolver(mtype=GridapPardiso.MTYPE_REAL_SYMMETRIC_INDEFINITE, msglvl=GridapPardiso.MSGLVL_VERBOSE)
        #solver = LinearFESolver(ps)
        #xh=solve(solver,op)
        #uh, vh, ωh, φh, ph = xh
        #ndofs=num_free_dofs(Gridap.FESpaces.get_fe_space(xh))
        xh=solve(op)
        
        Ω = Triangulation(model)
        dΩ = Measure(Ω,2*(order+2))
        Λ  = SkeletonTriangulation(model)
        Σ = BoundaryTriangulation(model,tags = "Sigma")
        dΛ = Measure(Λ,2*(order+2)-1)
        dΣ = Measure(Σ,2*(order+2)-1)    
        h_e = CellField(get_array(∫(1)dΛ),Λ)
        h_e_Σ = CellField(get_array(∫(1)dΣ), Σ)
    
        error=
          compute_errors_3D(xh, dΩ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

        erroru,errorv,errorω,errorφp=
          compute_B3_error_norms(xh, op, dΩ, dΛ, dΣ, h_e, h_e_Σ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

        ndofs=num_free_dofs(Gridap.FESpaces.get_fe_space(xh))+sum(num_dirichlet_dofs.(op.trial))
        push!(nn,ndofs)
        push!(hh,(sqrt(6.0)/2.0)/(2^(nk)+1)) #for regular tetrahedra: twice the circumradius
        println("******** Total DoFs: ", nn[nk])
        push!(eu,erroru)
        push!(ev,errorv)
        push!(eω,errorω)
        push!(eφp,errorφp)
        push!(etot,sqrt(sum([erroru,errorv,errorω,errorφp].^2)))
        if nk>1
            push!(ru , log(eu[nk]/eu[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rv , log(ev[nk]/ev[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rω , log(eω[nk]/eω[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rφp, log(eφp[nk]/eφp[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rtot, log(etot[nk]/etot[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("====================================================================================================")
    println("   DoF  &    h   &   e(u)   &  r(u)   &   e(v)   &  r(v)   &   e(ω)   &  r(ω)   &   e(φp)   &  r(φp)")
    println("====================================================================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f\n", nn[nk], hh[nk], eu[nk], ru[nk], ev[nk], rv[nk], eω[nk], rω[nk], eφp[nk], rφp[nk]);
    end
    println("===================================")

    println("===================================")
    println("   DoF  &    h   &   e(u)   &  r(u)")
    println("===================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f\n", nn[nk], hh[nk], etot[nk], rtot[nk]);
    end
    println("====================================================================================================")
end

# Arbitrary model parameters
const μ   = 10.0
const λ   = 1.0e+02
const ν   = 0.1
const κ   = 1.0e-03
const α   = 0.1
const c_0 = 0.1

nkmax = 5
order = 0
convergence_test_3D(order, μ, λ, ν, κ, α, c_0, nkmax)

order = 1
convergence_test_3D(order, μ, λ, ν, κ, α, c_0, nkmax)


