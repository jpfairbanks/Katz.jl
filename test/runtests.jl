using Katz
using Base.Test
using LightGraphs
using GraphMatrices
import GraphMatrices: KatzAdjacency
a = 0.15
g = LightGraphs.Datasets.BullGraph()
function test1(g, a::Number)
    n = nv(g)
    A = sparse(CombinatorialAdjacency(g))
    v = stdbasis(n, 1)
    k = katz(A, a, v)
    k2 = katz(A, a, ones(n))
    @test_approx_eq_eps k2 k 1e-5
end

function test_iterative(g, a)
    K = KatzAdjacency(CombinatorialAdjacency(g), a)
    n = nv(g)
    y = K*ones(n)
    z, convhist = ikatz(K,ones(n))
    x = K.A*z
    @show convhist.isconverged
    x_lu = katz(sparse(K.A), a, ones(n))
    @test norm(x_lu - x) < 1e-6
    @show ρ(convhist)
end
datasets = [:BullGraph,
            :TutteGraph,
            :DesarguesGraph]
for a in [0.1,0.01,0.001]
for s in datasets
    g = eval(:(LightGraphs.Datasets.$s()))
    println("test1 on $s, $a")
    test1(g, a)
    test_iterative(g,a)
end
end

for s in ["Newman/football"]
    g = MDGraph(s,:read)
    test1(g, α)
    test_iterative(g,α)
end
