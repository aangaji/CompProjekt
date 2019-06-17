function HnnISING(i::Integer; ϕi = Sys.ϕ[i])
    global Sys
    global J
    global h
    ϕ,nnEnvs = Sys.ϕ,Sys.nnEnvs

    nns = nnEnvs[i]
    -J/4*sum(sign(ϕ[nn])*sign(ϕi) for nn in nns) - h/2*sign(ϕi)
end

function metropolisStepISING()
    global Sys
    global β
    ϕ,dim = Sys.ϕ,Sys.dim
    i = rand(1:prod(dim))
    ϕi = -ϕ[i]
    if rand() < exp(β*(HnnISING(i)-HnnISING(i,ϕi=ϕi)))
        ϕ[i] = ϕi
    end
end

metropolisStepISING(n::Integer) = begin
    n == 1 ? nothing : metropolisStepISING(n-1)
    metropolisStepISING()
end
;
