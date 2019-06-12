#= about System: mutable structure containing fields.
When instantiated it also sets the global var Sys to the current hash.
Thus Sys only redirects to the System object of its current hash, mutating its fields=#

mutable struct System
    #input fields
    dim
    conditions
    latvecs
    #optional
    ϕ

    nnEnvs
    latpoints

    System(dim,conditions,latvecs,ϕ) = begin
        nnEnvsp = Array{Array{Int64,1}}(undef,dim)
        Ny = dim[1]
        for i in 1:prod(dim)
            nnEnvsp[i] = nnEnv(dim,conditions,i)
        end #*
        Xp,Yp= lattice(dim,latvecs)
        global Sys = new(dim,conditions,latvecs,ϕ,nnEnvsp,(Xp,Yp))
    end
    System(dim,conditions,latvecs) = System(dim,conditions,latvecs,2*pi.*rand(dim))
end
#= *collect nearest neighbour environments such that nnEnvs[i] is nnEnv @ site i
choose from cubicObc,cubicPbc,triangObc,triangPbc =#
##############################################################
# single site classical XY-model
## idea calculate energy for given nearest neighbours environment (nn::collection of indices) of single site (i::index) on the spin map (ϕ::collection of Float64)

## !!! carefull: this doesnt return the total energy by summation since couplings then get counted twice. However it is correct for comparing a change on site i !!!
function Hnn(i::Integer; ϕi = Sys.ϕ[i])
    global Sys
    global J
    global h
    ϕ,nnEnvs = Sys.ϕ,Sys.nnEnvs

    nns = nnEnvs[i]
    -J*sum(cos(ϕ[nn]-ϕi) for nn in nns) - h*cos(ϕi)
end
##############################################################
function Hmatrix(H)
    global Sys
    global J
    global h
    dim,ϕ,nnEnvs = Sys.dim,Sys.ϕ,Sys.nnEnvs

    for i in 1:prod(dim)
        H[i] = -J/2*sum(cos(ϕ[nn]-ϕ[i]) for nn in nnEnvs[i]) - h*cos(ϕ[i])
    end
end

##############################################################
# Metropolis step

function metropolisStep()
    global Sys
    dim = Sys.dim
    global β

    i = rand(1:prod(dim))
    ϕi = 2*pi*rand()
    if rand() < exp(β*(Hnn(i)-Hnn(i,ϕi=ϕi)))
        Sys.ϕ[i] = ϕi
    end
end

metropolisStep(n::Integer) = begin
    n == 1 ? nothing : metropolisStep(n-1)
    metropolisStep()
end
##############################################################
;
