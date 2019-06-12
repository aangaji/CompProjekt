##############################################################
function nnEnv(dim,conditions,i::Integer)
    Ny,Nx = dim
    nn = Array{Int64}(undef,0)
    notBottom :: Bool = (rem(i,Ny) != 0 && Ny != 1)
    notTop:: Bool = (rem(i,Ny) != 1 && Ny != 1)
    notRight :: Bool = (i <= (Nx-1)*Ny && Nx != 1)
    notLeft :: Bool = (i > Ny && Nx != 1)
    nots = (notTop,notBottom,notRight,notLeft)
    conditions(i,nn,(Nx,Ny),nots)
    nn
end

function cubicObc(i::Integer,nn::Array{Int64,1},(Nx,Ny)::Tuple{Integer,Integer},nots::Tuple{Bool,Bool,Bool,Bool})
    notTop,notBottom,notRight,notLeft = nots
    if notBottom push!(nn, i+1) end
    if notTop push!(nn, i-1) end
    if notLeft push!(nn, i-Ny) end
    if notRight push!(nn, i+Ny) end
end

function cubicPbc(i::Integer,nn::Array{Int64,1},(Nx,Ny)::Tuple{Integer,Integer},nots::Tuple{Bool,Bool,Bool,Bool})
    notTop,notBottom,notRight,notLeft = nots
    if notBottom push!(nn, i+1) else push!(nn, i-Ny+1) end
    if notTop push!(nn, i-1) else push!(nn, i+Ny-1) end
    if notLeft push!(nn, i-Ny) else push!(nn, i+(Nx-1)*Ny) end
    if notRight push!(nn, i+Ny) else push!(nn, i-(Nx-1)*Ny) end
end

function triangObc(i::Integer,nn::Array{Int64,1},(Nx,Ny)::Tuple{Integer,Integer},nots::Tuple{Bool,Bool,Bool,Bool})
    notTop,notBottom,notRight,notLeft = nots
    if notTop push!(nn, i-1)
        if notRight push!(nn, i+Ny-1) end
    end
    if notBottom push!(nn, i+1)
        if notLeft push!(nn,i-Ny+1) end
    end
    if notLeft push!(nn, i-Ny) end
    if notRight push!(nn, i+Ny) end
end

function triangPbc(i::Integer,nn::Array{Int64,1},(Nx,Ny)::Tuple{Integer,Integer},nots::Tuple{Bool,Bool,Bool,Bool})
    notTop,notBottom,notRight,notLeft = nots
    if notTop && notBottom && notRight && notLeft
        push!(nn, i-Ny,i-Ny+1,i-1,i+1,i+Ny-1,i+Ny)
    else
        if i==1 push!(nn, 2,Ny+1,Ny,2*Ny,(Nx-1)*Ny+1,(Nx-1)*Ny+2)
        elseif i==Ny push!(nn, 1,Ny-1,2*Ny-1,2*Ny,Ny*Nx,(Nx-1)*Ny+1)
        elseif i==(Nx-1)*Ny+1 push!(nn, i-Ny,i-Ny+1,i+1,1,Nx*Ny,Ny)
        elseif i==Nx*Ny push!(nn, i-Ny,i-1,i-Ny+1,Ny,Ny-1,i-2*Ny+1)
        elseif !notLeft push!(nn, i-1,i+1,i+Ny,i+Ny-1,i+(Nx-1)*Ny,i+(Nx-1)*Ny+1)
        elseif !notRight push!(nn, i-1,i+1,i-Ny,i-Ny+1,i-(Nx-1)*Ny,i-(Nx-1)*Ny-1)
        elseif !notTop push!(nn, i-Ny,i-Ny+1,i+1,i+Ny,i+Ny-1,i+2*Ny-1)
        elseif !notBottom push!(nn, i-Ny,i-Ny+1,i-1,i+Ny,i+Ny-1,i-2*Ny+1)
        end
    end
end
##############################################################
;
