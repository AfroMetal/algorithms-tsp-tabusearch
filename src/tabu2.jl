using Distances
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

####### SETTINGS #######
const itermax = Integer(floor((3_000_000*sqrt(n))/n^2))
println(string("iterations: ", itermax))
const penalty = 10
########################
tabulist = fill(-penalty, n)

frequencylist = zeros(Float64, n)
frequencylist[1] = Inf

currentsol = Array{Int32}(n+1)
currentsol[1] = currentsol[n+1] = Int32(1)
bestsol = Array{Int32}(n+1)
bestsol[1] = bestsol[n+1] = 1

function tab(i::Integer, j::Integer, iteration::Integer)
    tabulist[currentsol[i]] = iteration
    tabulist[currentsol[j]] = iteration
end

function istabbed(i::Integer, j::Integer, iteration::Integer)
    return (iteration - tabulist[currentsol[i]]) < penalty || (iteration - tabulist[currentsol[j]]) < penalty
end

function distance(solution::Array{Int32})
    cost = 0
    for i in 1:length(solution)-1
        cost += d[solution[i],solution[i+1]]
    end
    return cost
end

function distancediff(tour::Array{Int32}, j::Integer, k::Integer)
    return sum([
            -d[tour[j-1], tour[j]],
            -d[tour[j], tour[j+1]],
             d[tour[j-1], tour[k]],
             d[tour[k], tour[j+1]],
            -d[tour[k-1], tour[k]],
            -d[tour[k], tour[k+1]],
             d[tour[k-1], tour[j]],
             d[tour[j], tour[k+1]]
        ])
end

function minfreq(freq::Array{Float64})
    _, i = findmin(freq)
    freq[i] = Inf
    _, j = findmin(freq)

    return (i, j)
end
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

print("initial cost: ")
println(bestcost)
# print("initial route: ")
# println(bestsol)
println(string("initialization time: ", toq(), "s"))

noimprov = 0
checked = 0

neighborscount = 0.5 * (length(currentsol)-3) * (length(currentsol)-2)

for iter in 1:itermax
    neighborsol = copy(bestneighborsol)
    neighborcost = bestneighborcost

    for j in shuffle(2:length(neighborsol)-2), k in shuffle(j+1:length(neighborsol)-1)
        if k-j > 1 && !istabbed(j,k,iter)
            # checked += 1

            currentcost += distancediff(currentsol, j, k)

            currentsol[k], currentsol[j] = currentsol[j], currentsol[k]

            frequencylist[currentsol[j]] += currentcost - bestcost
            frequencylist[currentsol[k]] += currentcost - bestcost

            # TODO: tab every neighbor, not only ones checked
            tab(j, k, iter)

            if currentcost < bestneighborcost
                bestneighborsol = copy(currentsol)
                bestneighborcost = currentcost
            end

            currentsol = copy(neighborsol)
            currentcost = neighborcost
        end
        if (checked+=1) >= neighborscount/10
            checked = 0
            break
        end
    end

    if bestneighborcost < bestcost
        bestsol = copy(bestneighborsol)
        bestcost = bestneighborcost
        noimprov = 0
    else
        if (noimprov+=1) >= itermax/10
            (freq1, freq2) = minfreq(frequencylist)
            tour1 = findfirst(bestneighborsol.==freq1)
            tour2 = findfirst(bestneighborsol.==freq2)

            # bestneighborcost += distancediff(bestneighborsol, tour1, tour2)
            bestneighborsol[tour1], bestneighborsol[tour2] = bestneighborsol[tour2], bestneighborsol[tour1]
            bestneighborcost = distance(bestneighborsol)

            frequencylist = zeros(Float64, n)
            frequencylist[1] = Inf
            frequencylist[bestneighborsol[tour1]] += bestneighborcost - bestcost
            frequencylist[bestneighborsol[tour2]] += bestneighborcost - bestcost

            noimprov = 0

            tab(tour1, tour2, iter)

            if bestneighborcost < bestcost
                bestsol = copy(bestneighborsol)
                bestcost = bestneighborcost
            end
        end
    end
end
println(STDOUT, bestcost)
# for n in bestsol
#     print(STDERR, n)
#     print(STDERR, " ")
# end
# println(STDERR)

println(string("total time: ", toq(), "s"))
