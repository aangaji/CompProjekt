using PyPlot
pygui(true)
#################################
function lattice(latvecs) :: Tuple{Array{Float64,2},Array{Float64,2}}
    global dim
    g1,g2 = latvecs
    X,Y = Array{Float64}(undef,dim),Array{Float64}(undef,dim)
    for x=0:dim[1]-1, y=0:dim[2]-1
            X[x+1,y+1],Y[x+1,y+1] = x.*g1.+y.*g2
    end
    return X,Y
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

## !!! carefull: this doesnt return the total energy by summation since couplings then get counted twice. However it is correct for comparing a change on site i !!!
function Hnn(i::Integer; ϕi = ϕ[i])
    global dim
    global J
    global h
    global ϕ
    global nnEnvs #nearest neighbour environment matrix

    nns = nnEnvs[i]
    -J*sum(cos(ϕ[nn]-ϕi) for nn in nns) - h*cos(ϕi)
end
##############################################################


##############################################################
# Metropolis step

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
    metropolisStep()
end
##############################################################

##############################################################
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
##############################################################



dim = Nx,Ny = 10,10

J = 1
h = 0.
ϕ = 2*pi.*rand(dim)

# collect nearest neighbour environments such that nnEnvs[i] is nnEnv @ site i
nnEnvs = Array{Array{Int64,1}}(undef,dim)
for i in 1:prod(dim)
    nnEnvs[i] = nnEnv(i)
end
nnEnvs


g1,g2 = (1.,0.5),(.5,1.)
X,Y= lattice((g1,g2))

arrowmap()
