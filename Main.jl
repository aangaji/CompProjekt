using PyPlot
pygui(true)

dim = Nx,Ny = 10,10

#################################
function lattice(latvecs,dim::Tuple{Integer,Integer}) :: Array{typeof(latvecs[1]),2}
    g1,g2 = latvecs
    lat = Array{typeof(g1)}(undef,dim)
    for x=0:dim[1]-1, y=0:dim[2]-1
            lat[x+1,y+1] = x.*g1.+y.*g2
    end
    return lat
end
##################################
##################################
function index2Tuple(i::Integer)
    global Ny
    x,y = divrem(i,Ny).+(1,0)
    if y == 0 x-=1; y+=Ny end
    y,x
end
function tuple2Index(t::Tuple{Integer,Integer})
    global Ny
    y,x = t
    i = (x-1)*Ny+y
end
##################################
##################################
# nearest neighbour environment function for square lattice with obc
function nnEnv(i::Integer)
    global Nx
    global Ny
    nn = Array{Int64}(undef,0)
    if rem(i,Ny) != 0 push!(nn, i+1) end
    if rem(i,Ny) != 1 push!(nn, i-1) end
    if i > Ny push!(nn, i-Ny) end
    if i <= (Nx-1)*Ny push!(nn, i+Ny) end
    nn
end
################################

##############################################################
# single site classical XY-model
## idea calculate energy for given nearest neighbours environment (nn::collection of indices) of single site (i::index) on the spin map (ϕ::collection of Float64)
function Hnn(i::Integer)
    global dim
    global J
    global h
    global ϕ
    global nnEnvs #nearest neighbour environment matrix

    nns = nnEnvs[i]
    -J*sum(cos(ϕ[nn]-ϕ[i]) for nn in nns) - h*cos(ϕ[i])
end
##############################################################



J = 1
h = 0.
ϕ = 2*pi.*rand(dim)

# collect nearest neighbour environments such that nnEnvs[i] is nnEnv @ site i
nnEnvs = Array{Array{Int64,1}}(undef,dim)
for i in 1:prod(dim)
    nnEnvs[i] = nnEnv(i)
end
nnEnvs
# collect all single site Energies
H = Array{Float64}(undef,dim)
for i in 1:prod(dim)
    H[i] = Hnn(i)
end
H
sum(H)
