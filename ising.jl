function HnnISING(i::Integer; ϕi = ϕ[i])
    global dim
    global J
    global h
    global ϕ
    global nnEnvs #nearest neighbour environment matrix

    nns = nnEnvs[i]
    -J/4*sum(sign(ϕ[nn])*sign(ϕi) for nn in nns) - h/2*sign(ϕi)
end

function metropolisStepISING()
    global ϕ
    global β
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
