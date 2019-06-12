include("neighborhoods.jl")
include("usefullfunctions.jl")
include("model.jl")
include("ising.jl")


###########################################################################
J = 1.
h = 0.
β = .001

# reducing to 1D case by setting Ny=1 for high Nx shows the existence of continuous domain walls
# very small systems show ordered phase due to low entropy
cubicSys = System((20,40),cubicObc,((1.,0.),(0.,1.)))
triangSys = System((30,30),triangPbc,((1.,0.),(0.5,-1.)))
# Sys = cubicSys # All functions see mutable global var Sys
##############################################################


mz() = begin global Sys; cos.(Sys.ϕ) end
Hi = similar(Sys.ϕ)
hlocal() = begin global Hi; Hmatrix(Hi); return Hi end

scale = nothing
color = mz
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
Sys.ϕ = pi/2.*rand([1,-1],dim)
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
