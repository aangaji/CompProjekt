using PyPlot
pygui(true)
# lattice
function lattice(latvecs) :: Tuple{Array{Float64,2},Array{Float64,2}}
    global dim
    g1,g2 = latvecs
    X,Y = Array{Float64}(undef,dim),Array{Float64}(undef,dim)
    for x=0:dim[1]-1, y=0:dim[2]-1
            X[x+1,y+1],Y[x+1,y+1] = x.*g1.+y.*g2
    end
    return X,Y
end

#test
X,Y= lattice(((1.,0.5),(.5,1.)))
scatter(X,Y);title("Lattice")

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

function Hnn(i::Integer; ϕi = ϕ[i])
    global dim
    global J
    global h
    global ϕ
    global nnEnvs #nearest neighbour environment matrix

    nns = nnEnvs[i]
    -J*sum(cos(ϕ[nn]-ϕi) for nn in nns) - h*cos(ϕi)
end

nnEnvs = Array{Array{Int64,1}}(undef,dim)
for i in 1:prod(dim)
    nnEnvs[i] = nnEnv(i)
end
nnEnvs


Hnn(15)
Hnn(15; ϕi = 2*pi*rand())

function metropolisStep()
    global ϕ
    global β
    i = rand(1:prod(dim))
    ϕi = 2*pi*rand()
    if rand() < min(1,exp(β*(Hnn(i)-Hnn(i,ϕi))))
        ϕ[i] = ϕi
    end
end

metropolisStep(n::Integer) = begin
    n == 1 ? nothing : metropolisStep(n-1)
    println(n)
end

@time metropolisStep(20)
@time for k in 1:20 println(k) end

X,Y = collect(1:10),zeros(Int64,10)
phi = 2*pi.*rand(10)
U,V = cos.(phi),sin.(phi)

begin
    fig, ax = plt.subplots()
    Q = ax.quiver(X, Y, U, V, pivot="mid")
    scatter(X,Y,s=1.,color="red")
    PyPlot.show()
end

g1,g2 = (1.,0.5),(.5,1.)
X,Y= lattice((g1,g2))

SpinX,SpinY = cos.(ϕ),sin.(ϕ)
scatter(X,Y)

begin
    fig, ax = plt.subplots()
    Q = ax.quiver(X, Y, SpinX, SpinY, pivot="mid")
    scatter(X,Y,s=1.,color="red")
    PyPlot.show()
end

function arrowmap()
    global X
    global Y
    global ϕ
    SpinX,SpinY = cos.(ϕ),sin.(ϕ)
    fig, ax = plt.subplots()
    ax.quiver(X, Y, SpinX, SpinY, pivot="mid")
    scatter(X,Y,s=1.,color="red")
    PyPlot.show()
end

g1,g2 = (1.,0.5),(.5,1.)
X,Y= lattice((g1,g2))
ϕ = 2*pi.*rand(dim)
arrowmap()
