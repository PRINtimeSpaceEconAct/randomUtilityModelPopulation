

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
        title("T = $t")

        colorbar()
        # savefig("pics/vieEmilia$t.eps")
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

        title("T = $t")
        # savefig("pics/2AD$t.eps")
    end
end

function plotFunction(I,p)
    ratio = p.Nx/p.Ny
    fig = figure(figsize=(10,10/ratio))

    imshow(I',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    colorbar()

end

function plotAll(sol,p; t = p.T_end)
    ratio = p.Nx/p.Ny
    fig = figure(figsize=(15,10/ratio))

    l = sol(t)
    Wₕᴾl = imfilter(l,p.WₕᴾM,Fill(0,l))           # Wₕᴾ * l 
    Al = p.GM .* Wₕᴾl                             # local technological progress
    w =  Al .* abs.(l).^(p.β-1) .* (l .> 0 )      # wages
    Y =  Al .* abs.(l).^p.β                       # production
    AEN =  (p.τ * abs.(Y)).^p.φ .- p.γA * l
    Utility = p.γw * w + p.γEN * AEN

    subplot(231)
    imshow(l',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    colorbar()
    title("l")

    subplot(232)
    imshow(w',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    colorbar()
    title("w")

    subplot(233)
    imshow(Y',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    colorbar()
    title("Total income")

    subplot(234)
    imshow(AEN',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    colorbar()
    title("A endogenous")

    subplot(235)
    imshow(Al',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    colorbar()
    title("Al")

    subplot(236)
    imshow(Utility',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    colorbar()
    title("Utility")

    tight_layout()
end

