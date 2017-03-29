using Distances, StatsBase

dist = Euclidean()

tic()
tic()
# file lines iterator object
input = eachline(STDIN)

state = start(input)
# read data size N
(n, state) = next(input, state)
n = parse(Int32, n)

# 2xN node's coordinates matrix
# Data[:, k] = [x_k, y_k]
data = Array{Float64}(2, n)
# NxN node's distances matric
d = Array{Float64}(n, n)

# read file data to Data matrix
while !done(input, state)
    (v, state) = next(input, state)
    (i_v, x_v, y_v) = split(v)
    ii = parse(Int32, i_v)
    xx = parse(Float64, x_v)
    yy = parse(Float64, y_v)
    data[:,ii] = [xx, yy]
end
# calculate nodes distances as NxN matrix
pairwise!(d, dist, data)

# show(d)

function distance(solution::Array{Int32})
    cost = 0
    for i in 1:length(solution)-1
        cost += d[solution[i],solution[i+1]]
    end
    return cost
end

frequencylist = zeros(Float64, n)
frequencylist[1] = Inf

currentsol = Array{Int32}(n+1)
currentsol[1] = currentsol[n+1] = Int32(1)
bestsol = Array{Int32}(n+1)
bestsol[1] = bestsol[n+1] = 1

#
# for i in 2:n
#     currentsol[i] = Int32(i)
# end

for i in 1:n
    d[i,i] = Inf
end

d1 = copy(d)
k = 1
for i in 2:length(currentsol)-1
    currentsol[i] = findfirst(d1[k,:].==min(d1[k,:]...))
    k = currentsol[i]
    d1[k,1] = Inf
    for visited_node in 1:i
        d1[k, currentsol[visited_node]] = Inf
    end
end

currentcost = distance(currentsol)

bestsol = copy(currentsol)
bestcost = currentcost

bestneighborsol = copy(currentsol)
bestneighborcost = currentcost

neighborsol = copy(bestneighborsol)
neighborcost = bestneighborcost

N = [(j, k) for j=2:length(neighborsol)-2 for k=j+1:length(neighborsol)-1]

neighborscount = length(N)

function istabbed(i::Integer, j::Integer, iteration::Integer)
    return (iteration - tabulist[i]) < penalty || (iteration - tabulist[j]) < penalty
end

function isaspiring(i::Integer, j::Integer, current::Float64, base::Float64, iteration::Integer)
    return base > current + current*(iteration - max(tabulist[i], tabulist[j]))/penalty
end

####### SETTINGS #######
const itermax = Integer(floor((10_000_000*sqrt(n))/n^2))
println(string("max iterations:\t", itermax))

const penalty = round(sqrt(n)) # rounds tabbed

const neighborstocheck = sqrt(neighborscount)
########################

tabulist = fill(-penalty, n)
checks = 0

println(string("initial cost:\t", bestcost))
# println(string("initial route: ", bestsol))
println(string("initial time:\t", toq(), "s"))

for iter in 1:itermax
    neighborsol = copy(bestneighborsol)
    neighborcost = bestneighborcost

    checks = 0
    maxchecks = neighborstocheck
    city1 = 0
    city2 = 0

    for (j, k) in shuffle(N)
        c1 = neighborsol[j]
        c2 = neighborsol[k]
        if checks < maxchecks
            currentcost = neighborcost + sum([
                    -d[neighborsol[j-1], neighborsol[j]],
                    -d[neighborsol[j], neighborsol[j+1]],
                     d[neighborsol[j-1], neighborsol[k]],
                     d[neighborsol[k], neighborsol[j+1]],
                    -d[neighborsol[k-1], neighborsol[k]],
                    -d[neighborsol[k], neighborsol[k+1]],
                     d[neighborsol[k-1], neighborsol[j]],
                     d[neighborsol[j], neighborsol[k+1]]
                ])
            if !istabbed(c1, c2, iter) || isaspiring(c1, c2, currentcost, neighborcost, iter)
                # println(currentcost)

                if currentcost <= bestneighborcost
                    city1 = j
                    city2 = k
                    bestneighborcost = currentcost
                end
                checks += 1
            end
        elseif city1 == 0 || city2 == 0
            maxchecks *= 2
        end

        tabulist[c1] = iter
        tabulist[c2] = iter
        # tab(neighborsol[j], neighborsol[k], iter)
    end

    if city1 != 0 && city2 != 0
        bestneighborsol = copy(neighborsol)
        bestneighborsol[city2], bestneighborsol[city1] = bestneighborsol[city1], bestneighborsol[city2]
    else
        break
    end

    if bestneighborcost <= bestcost
        bestsol = copy(bestneighborsol)
        bestcost = bestneighborcost
    end
end

# println(STDOUT, bestcost)
# for n in bestsol
#     print(STDERR, n)
#     print(STDERR, " ")
# end
# println(STDERR)

println(STDOUT, string("best cost:\t", bestcost))
println(string("total time:\t", toq(), "s"))
