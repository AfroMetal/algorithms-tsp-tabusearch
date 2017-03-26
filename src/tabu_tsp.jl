using Distances
dist = Euclidean()

tic()
tic()
# file lines iterator object
input = eachline(STDIN)

state = start(input)
# read data size N
(n, state) = next(input, state)
n = parse(Int16, n)

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

d_orig = copy(d)
dim1 = size(d,1)
dim12 = size(d)

for i in 1:dim1
    d[i,i] = 10.0e6
end

d1 = copy(d)
tour = zeros(dim12)
cost = 0
min_dist = Array{Float32}(dim1)
short_path = Array{Int32}(dim1)
best_nbr_cost = 0
best_nbr = Array{Float64}(dim1)

k = 1
for i in 1:dim1-1
    min_dist[i] = min(d1[k,:]...)
    short_path[i] = findfirst(d1[k,:].==min(d1[k,:]...))
    cost = cost + min_dist[i]
    k = short_path[i]
    d1[k,1] = 10.0e6
    for visited_node in 1:i
        d1[k, short_path[visited_node]] = 10.0e6
    end
end

tour[1, short_path[1]] = 1
for i in 2:dim1-1
    tour[short_path[i-1], short_path[i]] = 1
end

short_path[dim1] = 1
tour[k, short_path[dim1]] = 1
cost = cost + d[k,1]

crnt_tour = copy(short_path)
best_tour = copy(short_path)
best_obj = cost
crnt_tour_cost = cost

print("initial cost: ")
println(best_obj)
println(string("initialization time: ", toq(), "s"))

nbr_cost = Array{Float64}(dim12)

tabu_tenure = zeros(dim12)
max_tabu_tenure = round(sqrt(dim1))

penalty = zeros(1, Integer((dim1-1)*(dim1-2)/2))
frequency = zeros(dim12)
frequency[1,:] = 100_000
frequency[:,1] = 100_000

for i in 1:dim1
    frequency[i,i] = 100_000
end

iter_snc_last_imprv = 0

best_nbr = copy(crnt_tour)
best_i = 0
best_j = 0

for iter in 1:Integer(round(300_000*sqrt(n)/n^2))
    # @printf("Iteration %d\n", iter)
    if crnt_tour_cost < 0
        @printf("ERROR! crnt_tour_cost < 0\n")
        println(crnt_tour)
        return
    end
    nbr_cost = fill(Inf, dim12)
    for i in 1:dim1-2
        for j in i+1:dim1-1
            if i == 1
                if j-1 == 1
                    nbr_cost[crnt_tour[i], crnt_tour[j]] = crnt_tour_cost - d[1,crnt_tour[i]] + d[1,crnt_tour[j]] - d[crnt_tour[j],crnt_tour[j+1]] + d[crnt_tour[i],crnt_tour[j+1]]
                    best_i = i
                    best_j = j
                    best_nbr_cost = nbr_cost[crnt_tour[i],crnt_tour[j]]
                    tabu_node1 = crnt_tour[i]
                    tabu_node2 = crnt_tour[j]
                else
                    nbr_cost[crnt_tour[i], crnt_tour[j]] = crnt_tour_cost - d[1,crnt_tour[i]] + d[1,crnt_tour[j]] - d[crnt_tour[j],crnt_tour[j+1]] + d[crnt_tour[i],crnt_tour[j+1]] - d[crnt_tour[i],crnt_tour[i+1]] + d[crnt_tour[j],crnt_tour[i+1]] - d[crnt_tour[j-1],crnt_tour[j]] + d[crnt_tour[j-1],crnt_tour[i]]
                end
            else
                if j-1 == 1
                    nbr_cost[crnt_tour[i],crnt_tour[j]] = crnt_tour_cost - d[crnt_tour[i-1],crnt_tour[i]] + d[crnt_tour[i-1],crnt_tour[j]] - d[crnt_tour[j],crnt_tour[j+1]] + d[crnt_tour[i],crnt_tour[j+1]]
                else
                    nbr_cost[crnt_tour[i],crnt_tour[j]] = crnt_tour_cost - d[crnt_tour[i-1],crnt_tour[i]] + d[crnt_tour[i-1],crnt_tour[j]] - d[crnt_tour[j],crnt_tour[j+1]] + d[crnt_tour[i],crnt_tour[j+1]] -
                    d[crnt_tour[i],crnt_tour[i+1]] + d[crnt_tour[j],crnt_tour[i+1]] - d[crnt_tour[j-1],crnt_tour[j]] + d[crnt_tour[j-1],crnt_tour[i]]
                end
            end

            if nbr_cost[crnt_tour[i],crnt_tour[j]] < best_nbr_cost
                best_nbr_cost = nbr_cost[crnt_tour[i],crnt_tour[j]]
                best_i = i
                best_j = j
                tabu_node1 = crnt_tour[i]
                tabu_node2 = crnt_tour[j]
                if tabu_node1 > tabu_node2
                    tabu_node1, tabu_node2 = tabu_node2, tabu_node1
                end
            end
        end
    end

    best_nbr[best_i] = crnt_tour[best_j]
    best_nbr[best_j] = crnt_tour[best_i]

    while tabu_tenure[tabu_node1,tabu_node2] > 0
        if best_nbr_cost < best_obj
            break
        else
            nbr_cost[tabu_node1,tabu_node2] = nbr_cost[tabu_node1,tabu_node2] * 1000
            best_nbr_cost_col = minimum(nbr_cost,2)
            best_nbr_cost = minimum(best_nbr_cost_col,1)[1]
            (R,C) = collect(zip(ind2sub(size(nbr_cost),find(nbr_cost.==best_nbr_cost))...))[1]
            tabu_node1 = R
            tabu_node2 = C
        end
    end

    if best_nbr_cost > crnt_tour_cost
        min_d_col = minimum(d,2)
        penal_nbr_cost = nbr_cost + minimum(min_d_col,1)[1]*frequency
        penal_best_nbr_cost_col = minimum(penal_nbr_cost,2)
        penal_best_nbr_cost = minimum(penal_best_nbr_cost_col,1)[1]
        (Rp,Cp) = collect(zip(ind2sub(size(penal_nbr_cost),find(penal_nbr_cost.==penal_best_nbr_cost))...))[1]
        tabu_node1 = Rp
        tabu_node2 = Cp
        best_nbr_cost = nbr_cost[tabu_node1,tabu_node2]
    end

    for row in 1:dim1-1
        for col in row+1:dim1
            if tabu_tenure[row,col] > 0
                tabu_tenure[row,col] = tabu_tenure[row,col] - 1
                tabu_tenure[col,row] = tabu_tenure[row,col]
            end
        end
    end

    tabu_tenure[tabu_node1,tabu_node2] = max_tabu_tenure
    tabu_tenure[tabu_node2,tabu_node1] = tabu_tenure[tabu_node1,tabu_node2]

    frequency[tabu_node1,tabu_node2] = frequency[tabu_node1,tabu_node2] + 1

    crnt_tour = copy(best_nbr)
    crnt_tour_cost = best_nbr_cost

    if crnt_tour_cost < best_obj
        best_obj = crnt_tour_cost
        # @printf("Iteration %d\n", iter)
        # @printf("best obj = %d\n", best_obj)
        best_tour = copy(crnt_tour)
        iter_snc_last_imprv = 0
    else
        iter_snc_last_imprv += 1
        if iter_snc_last_imprv >= 30*log2(n)
            # @printf("Iterations since last improvement exceeded...\n")
            min_freq_col = minimum(frequency,2)
            min_freq = minimum(min_freq_col,1)[1]
            (R,C) = collect(zip(ind2sub(size(frequency),find(frequency.==minimum(minimum(frequency,2),1)[1]))...))[1]
            freq_indx1 = R
            freq_indx2 = C
            indx_in_crnt_tour1 = findfirst(crnt_tour.==R)
            indx_in_crnt_tour2 = findfirst(crnt_tour.==C)
            # Diversify using a move that has the lowest frequency
            crnt_tour[indx_in_crnt_tour1], crnt_tour[indx_in_crnt_tour2] = crnt_tour[indx_in_crnt_tour2], crnt_tour[indx_in_crnt_tour1]
            tabu_tenure = zeros(dim12)
            frequency = zeros(dim12)
            frequency[1,:] = 100_000
            frequency[:,1] = 100_000
            for i in 1:dim1
                frequency[i,i] = 100_000
            end
            tabu_tenure[R,C] = max_tabu_tenure
            tabu_tenure[C,R] = max_tabu_tenure
            frequency[R,C] = frequency[R,C] + 1
            frequency[C,R] = frequency[R,C]
            # Re-calculare crnt tour cost
            crnt_tour_cost = d[1,crnt_tour[1]]
            for i in 1:dim1-1
                crnt_tour_cost = crnt_tour_cost + d[crnt_tour[i],crnt_tour[i+1]]
            end
            iter_snc_last_imprv = 0
            if crnt_tour_cost < best_obj
                best_obj = crnt_tour_cost
                best_tour = copy(crnt_tour)
            end
        end
    end
    # @printf("current tour cost = %d\t", crnt_tour_cost)
    # @printf("best obj = %d\n", best_obj)
end
println(STDOUT, best_obj)
for n in best_tour
    print(STDERR, n)
    print(STDERR, " ")
end
println(STDERR)

println(string("total time: ", toq(), "s"))
