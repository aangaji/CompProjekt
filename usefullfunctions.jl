##############################################################
update_figure = true
function handle_close(event)
    global update_figure
    update_figure = false
end

function arrowmap(;scale= nothing, showgrid=true, C = zeros(Float64,Sys.dim))
    global update_figure = true
    fig = figure("System")
    ax = plt.subplot()
    arrowmap((fig,ax),scale=scale, C=C)
    fig,ax
end
function arrowmap((fig,ax); scale= nothing, showgrid=true, C = zeros(Float64,Sys.dim))
    fig[:canvas][:mpl_connect]("close_event", handle_close)
    if update_figure
        global Sys
        global β
        X,Y = Sys.latpoints
        ϕ = Sys.ϕ
        U,V = ones(Float64,length(X)),ones(Float64,length(X))
        ax.cla()
        ax.quiver(X, Y, U, V, C, cmap = "brg", pivot="mid", scale=scale, angles = 360.*ϕ/(2*pi), alpha=1.)
        showgrid ? ax.scatter(X,Y,s=3.,color="black") : nothing
        ax.set_title("System T=$(round(1/β,digits=3))")
        PyPlot.show()
    end
end
##############################################################

##############################################################
function lattice(dim,latvecs) :: Tuple{Array{Float64,2},Array{Float64,2}}
    g1,g2 = latvecs
    X,Y = Array{Float64}(undef,dim),Array{Float64}(undef,dim)
    for x=0:dim[2]-1, y=0:dim[1]-1
            X[y+1,x+1],Y[y+1,x+1] = x.*g1.+y.*g2
    end
    return X,Y
end
##############################################################
;
