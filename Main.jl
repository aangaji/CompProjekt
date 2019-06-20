using PyPlot
pygui(true)

include("neighborhoods.jl")
include("usefullfunctions.jl")
include("model.jl")


###########################################################################
J = 1.
h = 0.
βcrit = π*J/2
β = .001

# reducing to 1D case by setting Ny=1 for high Nx shows the existence of continuous domain walls
# very small systems show ordered phase due to low entropy
linearSys = System((1,20),cubicPbc,((1.,0.),(0.,1.)))
cubicSys = System((150,150),cubicObc,((1.,0.),(0.,1.)))
triangSys = System((150,150),triangPbc,((1.,0.),(0.5,-1.)))
Sys = cubicSys # All functions see mutable global var Sys
##############################################################


mz() = begin global Sys; sin.(Sys.ϕ) end
hlocal() = begin Hmatrix(); return Hi end

scale = nothing
color = mz
fig,ax=arrowmap(scale=scale,C=color())
arrowmap((fig,ax); scale=scale,C=color())

#temper
function temper(βend)
    for βp in linspace(.01,βend,300)
        global β = βp
        metropolisStep(2000)
    end
end

temper(50.)
arrowmap((fig,ax); scale=scale,C=color(),showgrid=false)
#visualize
for Tp in 1/β:0.01:10.
    if !update_figure break end
    global β = 1/Tp
    metropolisStep(2000)
    arrowmap((fig,ax); scale=scale,C=color(),showgrid=false)
    sleep(.02)
end

################################################
#Measurements
begin
    color = hlocal
    Ts = collect(linspace(.1,3.,200))
    Es,dEs = similar(Ts),similar(Ts)
    Cvs,χs = similar(Ts),similar(Ts)
    Ms,dMs = similar(Ts),similar(Ts)
    βs = 1./Ts
    temper(βs[1])
    arrowmap((fig,ax); scale=scale,C=color())
    It = 150
    Et = Array{Float64}(It)
    Mt = Array{Float64}(It)
end
@time for t=1:length(βs)
    global β = βs[t]
    print(round(1/β,digits=2)," ")
    metropolisStep(1000)
    for i=1:It
        metropolisStep(1500)
        Et[i] = energy()
        Mt[i] = magnetization()
    end
    Es[t],dEs[t] = mean(Et),std(Et)/It
    Cvs[t] = Cv(Et)
    Ms[t],dMs[t] = mean(Mt),std(Mt)/It
    χs[t] = χ(Mt)
end

arrowmap((fig,ax); scale=scale,C=color())

measfig,axs = plt.subplots(2, 2)

begin
    axs[1,1].errorbar(Ts*βcrit,Es, yerr=sqrt.(dEs)); axs[1,1].set_title("Energy")
    axs[1,2].errorbar(Ts*βcrit,Ms, yerr=sqrt.(dMs)); axs[1,2].set_title("Magnetization")
    axs[2,1].plot(Ts*βcrit,Cvs); axs[2,1].set_title("C_v")
    axs[2,2].plot(Ts*βcrit,χs, label="L=$(prod(Sys.dim))"); axs[2,2].set_title("χ")
    legend()
    measfig.tight_layout()
end

for axi in axs
    axi.set_xscale("log")
    axi.set_xticks([.1,1,10])
    axi.set_xticklabels([.1,1,10])
end

for t=1:length(βs)-10
    Cvs[t] = mean(Cvs[t:t+5])
    χs[t] = mean(χs[t:t+5])
end






##############################################################
# ISING
include("ising.jl")

Sys = cubicSys
scale = nothing
color=hlocal
Sys.ϕ = pi/2.*rand([1,-1],Sys.dim)
fig,ax=arrowmap(; scale=scale,C=color())

β=.01
J = 1.
h = 10.

for βp in linspace(.001,5.,200)
    global β = βp
    metropolisStepISING(1000)
end
arrowmap((fig,ax); scale=scale,C=color())

metropolisStepISING(1000)
for βp in β:-0.05:.1
    if !update_figure break end
    global β = βp
    metropolisStepISING(1000)
    arrowmap((fig,ax); scale=scale,C=color())
    sleep(.05)
end
##############################################################
