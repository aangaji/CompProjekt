include("neighborhoods.jl")
include("usefullfunctions.jl")
include("model.jl")
include("ising.jl")


###########################################################################
J = 1.
h = 0.
βcrit = π*J/2
β = .001

# reducing to 1D case by setting Ny=1 for high Nx shows the existence of continuous domain walls
# very small systems show ordered phase due to low entropy
cubicSys = System((8,8),cubicObc,((1.,0.),(0.,1.)))
triangSys = System((30,30),triangPbc,((1.,0.),(0.5,-1.)))
Sys = triangSys # All functions see mutable global var Sys
##############################################################


mz() = begin global Sys; cos.(Sys.ϕ) end
hlocal() = begin Hmatrix(); return Hi end

scale = nothing
color = hlocal
fig,ax=arrowmap(scale=scale,C=color())
arrowmap((fig,ax); scale=scale,C=color())

#temper
function temper(βend)
    for βp in linspace(.001,βend,200)
        global β = βp
        metropolisStep(1000)
    end
end

temper(20.)
arrowmap((fig,ax); scale=scale,C=color())
#visualize
for βp in β:-0.1:.1
    if !update_figure break end
    global β = βp
    metropolisStep(2000)
    arrowmap((fig,ax); scale=scale,C=color())
    sleep(.02)
end

################################################
#Measurements
begin
    color = hlocal
    Ts = linspace(.001,10.,200)
    Es,dEs = similar(Ts),similar(Ts)
    Cvs,χs = similar(Ts),similar(Ts)
    Ms,dMs = similar(Ts),similar(Ts)
    βs = 1./Ts
    temper(βs[1])
    arrowmap((fig,ax); scale=scale,C=color())
    It = 80
    Et = Array{Float64}(It)
    Mt = Array{Float64}(It)
end
@time for t=1:length(βs)
    global β = βs[t]
    metropolisStep(1000)
    for i=1:It
        metropolisStep(1000)
        Et[i] = energy()
        Mt[i] = magnetization()
    end
    Es[t],dEs[t] = mean(Et),var(Et)
    Cvs[t] = Cv(Et)
    Ms[t],dMs[t] = mean(Mt),var(Mt)
    χs[t] = χ(Mt)
end
arrowmap((fig,ax); scale=scale,C=color())

begin
    measfig,axs = plt.subplots(2, 2)
    axs[1,1].errorbar(Ts*βcrit,Es, yerr=sqrt.(dEs)); axs[1,1].set_title("Energy")
    axs[1,2].errorbar(Ts*βcrit,Ms, yerr=sqrt.(dMs)); axs[1,2].set_title("Magnetization")
    axs[2,1].plot(Ts*βcrit,Cvs); axs[2,1].set_title("C_v")
    axs[2,2].plot(Ts*βcrit,χs); axs[2,2].set_title("χ")
    for axi in axs
        axi.set_xscale("log")
        axi.set_xticks([.001,.01,.1,1,10])
        axi.set_xticklabels([.001,.01,.1,1,10])
    end
    measfig.tight_layout()
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
