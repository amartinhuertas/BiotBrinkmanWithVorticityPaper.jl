struct GridapLinearSolverPreconditioner{A,B}
  numerical_setup::A
  matrix::B
  function GridapLinearSolverPreconditioner(matrix;ls=LUSolver())
    ss=symbolic_setup(ls,matrix)
    ns=numerical_setup(ss,matrix)
    A=typeof(ns)
    B=typeof(matrix)
    new{A,B}(ns,matrix)
 end
end

function LinearAlgebra.mul!(z::AbstractVector{T},
              M::GridapLinearSolverPreconditioner,
              r::AbstractVector{T}) where T
  # Solve Mz=r ...
  Gridap.solve!(z,M.numerical_setup,r)
end

function LinearOperators.LinearOperator(solver::GridapLinearSolverPreconditioner)
  function apply(res, x, α, β)
    if (β ≈ 0.0)
      mul!(res,solver,x)
    else 
      @assert α ≈ 1.0
      tmp=copy(res)
      mul!(tmp,solver,x)
      res .=  β .* res .+  α .* tmp
    end
  end
  LinearOperator(Float64, size(solver.matrix,1), size(solver.matrix,2), true, true, apply)
end 




