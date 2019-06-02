using PyPlot

# lattice
function lattice(latvecs,dim::Tuple{Integer,Integer}) :: Array{typeof(latvecs[1]),2}
    g1,g2 = latvecs
    lat = Array{typeof(g1)}(undef,dim)
    for x=0:dim[1]-1, y=0:dim[2]-1
            lat[x+1,y+1] = x.*g1.+y.*g2
    end
    return lat
end

#test
lat= lattice(((1.,0.5),(.5,1.)),(10,10))
figure("Lattice")
for (x,y) in lat
    scatter(x,y,color="black",s=1.)
end

# matrix indices
A = Array{Int64}(undef,10,10)
t=1
for x=1:10, y=1:10
    A[x,y] = t
    t+=1
end

B = Array{Int64}(undef,10,10)
t=1
for x=1:100
    B[x] = t
    t+=1
end
A,B
A[20],B[20]

function index2Tuple(i::Integer,Ny::Integer)
    x,y = divrem(i,Ny).+(1,0)
    if y == 0 x-=1; y+=Ny end
    y,x
end

# test
B = Array{Int64}(undef,10,10)
t=1
for x=1:100
    B[x] = t
    t+=1
end

B[20]
index2Tuple(20,Ny)
B[10,2]

function tuple2Index(t::Tuple{Integer,Integer},Ny::Integer)
    y,x = t
    i = (x-1)*Ny+y
end

tuple2Index(index2Tuple(10,10),10)

function nnEnv(i::Integer,dim::Tuple{Integer,Integer})
    Nx,Ny = dim
    nn = Array{Int64}(undef,0)
    if rem(i,Ny) != 0 push!(nn, i+1) end
    if rem(i,Ny) != 1 push!(nn, i-1) end
    if i > Ny push!(nn, i-Ny) end
    if i <= (Nx-1)*Ny push!(nn, i+Ny) end
    nn
end

nn= nnEnv(21,(10,10))
#
x=1
function y()
    global x
    return 3+x
end
y()
