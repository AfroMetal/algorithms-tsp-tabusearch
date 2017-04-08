using Distances, StatsBase

const DEBUG = length(ARGS) >=1 && ARGS[1] == "1"

function distance(solution::Array{Int32})
    cost = 0
    for i in 1:length(solution)-1
        cost += d[solution[i],solution[i+1]]
    end
    return cost
end

function tab(i::Integer, j::Integer, iteration::Integer)
    tabulist[i] = iteration
    tabulist[j] = iteration
end

function istabbed(i::Integer, j::Integer, iteration::Integer)
    return (iteration - tabulist[i]) < penalty || (iteration - tabulist[j]) < penalty
end

function isaspiring(i::Integer, j::Integer, current::Float64, base::Float64, itr::Integer)
    return base > current + 0.7(penalty-(itr-min(tabulist[i], tabulist[j])))
end

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

# read file data to Data matrix
while !done(input, state)
    (v, state) = next(input, state)
    (i_v, x_v, y_v) = split(v)
    data[:,parse(Int32, i_v)] = [parse(Float64, x_v), parse(Float64, y_v)]
end
# calculate nodes distances as NxN matrix

d = Array{Float64}(n, n)
pairwise!(d, Euclidean(1e-12), data)

currentsol = Array{Int32}(n+1)
currentsol[1] = currentsol[n+1] = Int32(1)

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

N = [(j, k) for j=2:length(bestsol)-2 for k=j+1:length(bestsol)-1 if k-j>1]

neighborscount = length(N)

####### SETTINGS #######
const itermax = round(Int32, 1.0e7*round(sqrt(n))/n^2)

const neighborstocheck = round(Int32, sqrt(neighborscount)/log(100, n))

const penalty = round(Int32, 0.7sqrt(n)) # rounds tabbed
########################

tabulist = fill(-penalty, n)
checks = 0

if DEBUG
    println(string("max iterations:\t", itermax))
    println(string("initial cost:\t", bestcost))
    # for n in bestsol
    #     println(STDERR, n)
    # end
    println(string("initial time:\t", toq(), "s"))
end

for iter in 1:itermax
    neighborsol = copy(bestneighborsol)
    neighborcost = bestneighborcost

    checks = 0
    maxchecks = neighborstocheck
    index1 = 0
    index2 = 0

    for (j, k) in shuffle(N)
        city1 = neighborsol[j]
        city2 = neighborsol[k]

        if checks < maxchecks
            currentcost = neighborcost - d[neighborsol[j-1], neighborsol[j]] - d[neighborsol[j], neighborsol[j+1]] + d[neighborsol[j-1], neighborsol[k]] + d[neighborsol[k], neighborsol[j+1]] - d[neighborsol[k-1], neighborsol[k]] - d[neighborsol[k], neighborsol[k+1]] + d[neighborsol[k-1], neighborsol[j]] + d[neighborsol[j], neighborsol[k+1]]

            if !istabbed(city1, city2, iter) || isaspiring(city1, city2, currentcost, neighborcost, iter)
                if currentcost <= bestneighborcost || index1 == 0
                    index1 = j
                    index2 = k
                    bestneighborcost = currentcost
                    if currentcost < bestcost
                        # maxchecks = min(maxchecks, checks+0.01neighborstocheck)
                        break
                    end
                end
                checks += 1
            end
        end
        tab(city1, city2, iter)
    end

    if index1 != 0 && index2 != 0
        bestneighborsol = copy(neighborsol)
        bestneighborsol[index2], bestneighborsol[index1] = bestneighborsol[index1], bestneighborsol[index2]

        if bestneighborcost <= bestcost
            bestsol = copy(bestneighborsol)
            bestcost = bestneighborcost
        end
    else
        break
    end
end

if DEBUG
    println(string("best cost:\t", bestcost))
    println(string("total time:\t", toq(), "s"))
else
    println(STDOUT, bestcost)
    for n in bestsol
        print(STDERR, n)
        print(STDERR, " ")
    end
    println(STDERR)
end
