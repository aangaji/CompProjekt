using PyPlot
pygui(true)
include("neighborhoods.jl")
include("usefullfunctions.jl")
include("model.jl")
include("ising.jl")

# reducing to 1D case by setting Ny=1 for high Nx shows the existence of continuous domain walls
# very small systems show ordered phase due to low entropy

dim = Ny,Nx = 20,40

J = 1
h = 0.
β = .001

# collect nearest neighbour environments such that nnEnvs[i] is nnEnv @ site i
## choose from cubicObc,cubicPbc,triangObc,triangPbc
conditions = cubicObc
nnEnvs = Array{Array{Int64,1}}(undef,dim)
for i in 1:prod(dim)
    nnEnvs[i] = nnEnv(conditions,i)
end
nnEnvs


g1,g2 = (1.,0.),(0.,1.)
X,Y= lattice((g1,g2))

Hi = similar(X)
##############################################################

mz() = begin global ϕ; cos.(ϕ) end
hlocal() = begin global Hi; Hmatrix(Hi); return Hi end

scale = nothing
ϕ = 2*pi.*rand(dim)
color = hlocal
fig,ax=arrowmap(scale=scale,C=color())
arrowmap((fig,ax); scale=scale,C=color())

#temper
βend = 10.
for βp in linspace(.001,βend,200)
    global β = βp
    metropolisStep(1000)
end

metropolisStep(1000)
for βp in β:0.1:50.
    if !update_figure break end
    global β = βp
    metropolisStep(1000)
    arrowmap((fig,ax); scale=scale,C=color())
    sleep(.05)
end

for n in 1:1000
    if !update_figure break end
    metropolisStep(20)
    arrowmap((fig,ax); scale=scale,C=color())
    sleep(.05)
end






##############################################################
# ISING
scale = nothing
ϕ = pi/2.*rand([1,-1],dim)
fig,ax=arrowmap(scale=scale)

metropolisStepISING(1000)
for βp in β:0.1:10.
    if !update_figure break end
    global β = βp
    metropolisStepISING(1000)
    arrowmap((fig,ax); scale=scale)
    sleep(.05)
end
##############################################################
