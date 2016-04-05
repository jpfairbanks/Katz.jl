module Katz
using LightGraphs
using GraphMatrices
using IterativeSolvers
export katz, stdbasis, ikatz, ρ

function katz(A, a, v)
    F = lufact(I-a*A)
    return A*(F\v)
end

function ikatz(A::KatzAdjacency, b::AbstractVector, tol::Real=1e-8, maxiter::Integer=100)
    x, convhist = cg(A, b, tol=tol, maxiter=maxiter)
    return x, convhist
end

function ρ(convhist::ConvergenceHistory)
    rs = @show convhist.residuals
    return mean([rs[i+1]/rs[i] for i in 1:length(rs)-1])
end

function stdbasis(n,i)
    x = ones(n)
    x[i] + 1
    return x
end

end # module
