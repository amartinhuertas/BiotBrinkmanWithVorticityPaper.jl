using BiotBrinkmanWithVorticityPaper
using Gridap
using Gridap.FESpaces
using Printf
using Test

function convergence_test(order, μ, λ, ν, κ, α, c_0, nkmax; 
                          u_ex=default_u_ex, 
                          p_ex=default_p_ex, 
                          v_ex=default_v_ex)

    eu = Float64[]
    ev = Float64[]
    eω = Float64[]
    ep = Float64[]
    eφ = Float64[]
    ru = Float64[]
    rv = Float64[]
    rω = Float64[]
    rp = Float64[]
    rφ = Float64[]
    lossn = Float64[]
    hh = Float64[]
    nn = Int[]

    push!(ru,0.)
    push!(rv,0.)
    push!(rω,0.)
    push!(rp,0.)
    push!(rφ,0.)

    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_model3d(nk)
        op=assemble_biotbrinkman(model, order, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)
        xh=solve(op)
        uh, vh, ωh, φh, ph = xh
        ndofs=num_free_dofs(Gridap.FESpaces.get_fe_space(xh))
        Ω = Triangulation(model)
        dΩ = Measure(Ω,2*(order+2))

        error_u,error_v,error_ω,error_φ,error_p,lossh=
          compute_errors_biotbrinkman(xh, dΩ, μ, λ, ν, κ, α, c_0, u_ex, p_ex, v_ex)

        push!(nn,ndofs)
        push!(hh,sqrt(3/2)/(2^(nk))) #for regular tetrahedra: twice the circumradius
        println("******** Total DoFs: ", nn[nk])
        push!(eu,error_u)
        push!(ev,error_v)
        push!(eω,error_ω)
        push!(eφ,error_φ)
        push!(ep,error_p)
        push!(lossn,lossh)
        if nk>1
            push!(ru, log(eu[nk]/eu[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rv, log(ev[nk]/ev[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rω, log(eω[nk]/eω[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rφ, log(eφ[nk]/eφ[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rp, log(ep[nk]/ep[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("============================================================================================================================")
    println("   DoF  &    h   &   e(u)   &  r(u) &   e(v)   &  r(v) &   e(ω)   &  r(ω) &   e(φ)   &  r(φ) &   e(p)   &  r(p) &  |loss|   ")
    println("============================================================================================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e \n", nn[nk], hh[nk], eu[nk], ru[nk], ev[nk], rv[nk],  eω[nk], rω[nk], eφ[nk], rφ[nk], ep[nk], rp[nk], lossn[nk]);
    end
    println("============================================================================================================================")
end

# Arbitrary model parameters
const μ   = 1.0
const λ   = 1.0
const ν   = 1.0
const κ   = 1.0
const α   = 1.0
const c_0 = 1.0

const nkmax = 3
const order = 0 

convergence_test(order, μ, λ, ν, κ, α, c_0, nkmax)
