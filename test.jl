mutable struct System
    #input fields
    dim
    conditions
    latvecs

    ϕ

    nnEnvs
    latpoints

    System(dim,conditions,latvecs,ϕ) = begin
        nnEnvsp = Array{Array{Int64,1}}(undef,dim)
        Ny = dim[1]
        for i in 1:prod(dim)
            nnEnvsp[i] = nnEnv(dim,conditions,i)
        end
        Xp,Yp= lattice(dim,latvecs)
        global Sys = new(dim,conditions,latvecs,ϕ,nnEnvsp,(Xp,Yp))
    end
    System(dim,conditions,latvecs) = System(dim,conditions,latvecs,2*pi.*rand(dim))
end

###########################################################################
cubicSys = System((20,40),cubicObc,((1.,0.),(0.,1.)))
triangSys = System((30,30),triangPbc,((1.,0.),(0.5,-1.)))

Sys
Sys = cubicSys

dim,ϕ,nnEnvs = Sys.dim,Sys.ϕ,Sys.nnEnvs
ϕ === Sys.ϕ
dim === Sys.dim

ϕ = rand(dim)

Sys.ϕ = rand(dim)

Sys === cubicSys
###########################################################################
J = 1.
h = 0.
β = .001

cubicSys = System((20,40),cubicObc,((1.,0.),(0.,1.)))

##############################################################

mz() = begin global Sys; cos.(Sys.ϕ) end
Hi = similar(Sys.ϕ)
hlocal() = begin global Hi; Hmatrix(Hi); return Hi end

scale = nothing
color = mz
fig,ax=arrowmap(scale=scale,C=color())
arrowmap((fig,ax); scale=scale,C=color())

metropolisStep()

Sys === cubicSys
Sys === triangSys

metropolisStep(1000)
for βp in β:0.1:50.
    if !update_figure break end
    global β = βp
    metropolisStep(1000)
    arrowmap((fig,ax); scale=scale,C=color())
    sleep(.05)
end
