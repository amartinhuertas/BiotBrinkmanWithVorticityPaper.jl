using BiotBrinkmanWithVorticityPaper
using Gridap
using Gridap.FESpaces
using Printf
using Test

function convergence_test_3D(order, μ, λ, ν, κ, α, c_0, nkmax; 
                             u_ex=default_3D_u_ex, 
                             p_ex=default_3D_p_ex, 
                             v_ex=default_3D_v_ex)

    e = Float64[]
    r = Float64[]
    hh = Float64[]
    nn = Int[]

    push!(r,0.)

    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_3D_model(nk)
        op=assemble_3D(model, order, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
        xh=solve(op)
        uh, vh, ωh, φh, ph = xh
        ndofs=num_free_dofs(Gridap.FESpaces.get_fe_space(xh))
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

        eu,ev,eω,eφp=compute_B3_error_norms(xh, op, dΩ, dΛ, dΣ, h_e, h_e_Σ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
        print("$(eu) $(ev) $(eω) $(eφp)")

        push!(nn,ndofs)
        push!(hh,sqrt(3.0/2.0)/(2^(nk))) #for regular tetrahedra: twice the circumradius
        println("******** Total DoFs: ", nn[nk])
        push!(e,error)
        if nk>1
            push!(r, log(e[nk]/e[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("===================================")
    println("   DoF  &    h   &   e(u)   &  r(u)")
    println("===================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f  \n", nn[nk], hh[nk], e[nk], r[nk]);
    end
    println("===================================")
end

# Arbitrary model parameters
const μ   = 1# 10.0
const λ   = 1#1.0e+02
const ν   = 1#1.0e-01
const κ   = 1#1.0e-03
const α   = 1#0.1
const c_0 = 1#0.1

nkmax = 3
order = 0
convergence_test_3D(order, μ, λ, ν, κ, α, c_0, nkmax)

# order = 1
# convergence_test_3D(order, μ, λ, ν, κ, α, c_0, nkmax)


