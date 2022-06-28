

function plotSol(sol,p; t_step = 0.1)
    
    ratio = p.Nx/p.Ny
    fig = figure(figsize=(10,10/ratio))
    t_plot = sol.t[1]:t_step:sol.t[end]

    vmin = minimum(log.(sol .+1))
    vmax = maximum(log.(sol .+1))

    # t_plot = [0.0,200.0,410.0,425.0,500.0,1000.0]
    for t in t_plot
        fig.clear()
        imshow(log.(sol(t).+1)',origin="lower",
            vmin = vmin, vmax = vmax, extent=(0.0,p.Lx,0.0,p.Ly))
        xticks(fontsize=15)
        yticks(fontsize=15)

        # title("T = $t")
        # colorbar()
        savefig("pics/vieEmilia$t.eps")
        pause(0.01)
    end
end

function plotSurfSol(sol,p; t_step = 0.1)
    ratio = p.Nx/p.Ny
    fig = figure(figsize=(10,10/ratio))
    t_plot = sol.t[1]:t_step:sol.t[end]

    zmin = 0.0
    zmax = maximum(log.(max.(sol.+1,0.0)))
    X = repeat(p.x,1,p.Ny)
    Y = repeat(p.y',p.Nx,1)
    for t in t_plot
        fig.clear()
        ax = fig.add_subplot(111, projection="3d")
        surf(X,Y,log.(sol(t).+1),rstride=5,cstride=5,facecolors=get_cmap("jet")(log.(sol(t).+1)))
        ax.axes.set_zlim3d(bottom=zmin, top=zmax) 
        # pause(0.1)

        # title("T = $t")
        # savefig("pics/2AD$t.eps")
    end
end

function plotMass(sol,p)
    @unpack_AD p 

    minX = 3.5
    maxX = 4.5 
    minY = y[1]
    maxY = y[Ny]

    mini = findfirst(x -> x >= minX,p.x)
    maxi = findfirst(x -> x >= maxX,p.x)
    minj = findfirst(x -> x >= minY,p.y)
    maxj = findfirst(x -> x >= maxY,p.y)
    

    times = sol.t
    mass = zeros(length(times))
    for (i,t) in enumerate(times)
        mass[i] = sum(sol(t)[mini:maxi,minj:maxj]) * Δx^2
    end
    plot(times,mass)
    xlabel("Period")
    ylabel("Mass")
    grid()
    
end
