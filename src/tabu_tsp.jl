using Distances
dist = Euclidean()

# file lines iterator object
I = eachline(STDIN)

state = start(I)
# read data size N
(n, state) = next(I, state)
n = parse(Int16, n)

# 2xN node's coordinates matrix
# Data[:, k] = [x_k, y_k]
Data = Array{Integer}(2, n)
# NxN node's distances matric
D = Array{Float64}(n, n)

# read file data to Data matrix
while !done(I, state)
    (v, state) = next(I, state)
    (i, x, y) = split(v)
    i = parse(Int32, i)
    x = parse(Int32, x)
    y = parse(Int32, y)
    Data[:,i] = [x, y]
end
# calculate nodes distances as NxN matrix
pairwise!(D, dist, Data)

# coordinates getting helper
# return x_k
function x(k::Integer)
    return Data[1,k]
end

# coordinates getting helper
# return y_k
function y(k::Integer)
    return Data[2,k]
end

# distance getting helper
# return d = sqrt(sum((x - y) .^ 2))
function d(l::Integer, m::Integer)
    return D[l,m]
end
