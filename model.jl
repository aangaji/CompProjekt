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
function Hmatrix(H)
    global dim
    global J
    global h
    global ϕ
    global nnEnvs

    for i in 1:prod(dim)
        H[i] = -J/2*sum(cos(ϕ[nn]-ϕ[i]) for nn in nnEnvs[i]) - h*cos(ϕ[i])
    end
end

##############################################################
# Metropolis step

function metropolisStep()
    global ϕ
    global β
    i = rand(1:prod(dim))
    ϕi = 2*pi*rand()
    if rand() < exp(β*(Hnn(i)-Hnn(i,ϕi=ϕi)))
        ϕ[i] = ϕi
    end
end

metropolisStep(n::Integer) = begin
    n == 1 ? nothing : metropolisStep(n-1)
    metropolisStep()
end
##############################################################
;
