

function plotInitialProjection(p)
    subplot(121)
    plot(p.y,vec(sum(p.AES,dims=1))/maximum(vec(sum(p.AES,dims=1))),label="AES")
    plot(p.y,vec(sum(p.GM,dims=1))/maximum(vec(sum(p.GM,dims=1))),label="G")
    plot(p.y,vec(sum(p.u₀,dims=1))/maximum(vec(sum(p.u₀,dims=1))),label="u0")
    plot([p.hₚ,p.hₚ],[0,1])
    grid()
    legend()

    subplot(122)
    plot(p.x,vec(sum(p.AES,dims=2))/maximum(vec(sum(p.AES,dims=2))),label="AES")
    plot(p.x,vec(sum(p.GM,dims=2))/maximum(vec(sum(p.GM,dims=2))),label="G")
    plot(p.x,vec(sum(p.u₀,dims=2))/maximum(vec(sum(p.u₀,dims=2))),label="u0")
    plot([p.hₚ,p.hₚ],[0,1])
    grid()
    legend()

end

function plotFunction(I,p; surfPlot=false)
    ratio = p.Nx/p.Ny
    fig = figure(figsize=(10,10))
    if surfPlot == false
        ct = contour(I',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
        clabel(ct,inline=1);
        # colorbar()
    else
        xGrid = repeat(p.x,1,p.Ny)
        yGrid = repeat(p.y',p.Nx,1)
        surf(xGrid,yGrid,I,rstride=1,cstride=1,facecolors=get_cmap("jet")(I))
        ax = gca()
        ax.set_yticks(0:p.Ly)
        ax.set_box_aspect((p.Lx,p.Ly,4))
    end
    # tight_layout()
end

function plotSol(sol,p; t_step = 0.1)
    
    t_plot = sol.t[1]:t_step:sol.t[end]

    global fig = figure(figsize=(15,max(10/(p.Nx/p.Ny),3.0)))
    for t in t_plot
        plotAll(sol,p; t = t)
        if !isdir(p.folder_name)
            mkdir(p.folder_name)
        end
        savefig("$(p.folder_name)/RUM$(lpad(string(round(t,digits=1)),5,"0")).png")
    end
end

function plotSurfSol(sol,p; t_step = 0.1)

    t_plot = sol.t[1]:t_step:sol.t[end]

    global fig = figure(figsize=(15,10))
    for t in t_plot
        plotAllSurf(sol,p; t = t)
        if !isdir(p.folder_name)
            mkdir(p.folder_name)
        end
        savefig("$(p.folder_name)/RUMSurf$(lpad(string(round(t,digits=1)),5,"0")).png")
    end
end


function plotAllSurf(sol,p; t = p.T_end, saveFig = false)
    subplots = [231,232,233,234,235,236]

    if t > 0
        for plt in subplots
            subplot(plt)
            fig = gca()
            lines = fig.get_lines()
            for line in lines
                line.remove()
            end
        end
    end

    xGrid = repeat(p.x,1,p.Ny)
    yGrid = repeat(p.y',p.Nx,1)

    l = sol(t)
    Wₕᴾl = imfilter(l,p.WₕᴾM,Fill(0,l))           # Wₕᴾ * l 
    Wₕᴾl[Wₕᴾl .< p.ϵTol] .= 0

    Al = p.GM .* Wₕᴾl                                           # local technological progress
    w =  Al .* fw.(abs.(l),p.ϵY,p.β)                            # wages
    Y =  Al .* fY.(abs.(l),p.ϵY,p.β)                            # production
    AEN =  p.GM .* ((p.τ^p.φ) * fY.(Y,p.ϵY,p.φ) .- p.γA * l)    # endogenous amenities 
    Utility = p.γw * w + p.γEN * AEN

    suptitle("t = $(round(t,digits=3))")

    ax = fig.add_subplot(231, projection="3d")
    surf(xGrid,yGrid,l,rstride=1,cstride=1,facecolors=get_cmap("jet")(l))
    # colorbar(shrink=0.7)
    title("Workers")

    ax = fig.add_subplot(232, projection="3d")
    surf(xGrid,yGrid,w,rstride=1,cstride=1,facecolors=get_cmap("jet")(w))
    # colorbar(shrink=0.7)
    title("Wages")

    ax = fig.add_subplot(233, projection="3d")
    surf(xGrid,yGrid,Y,rstride=1,cstride=1,facecolors=get_cmap("jet")(Y))
    # colorbar(shrink=0.7)
    title("Total income")

    ax = fig.add_subplot(234, projection="3d")
    surf(xGrid,yGrid,AEN,rstride=1,cstride=1,facecolors=get_cmap("jet")(AEN))
    # colorbar(shrink=0.7)
    title("Endogenous amenities")

    ax = fig.add_subplot(235, projection="3d")
    surf(xGrid,yGrid,Al,rstride=1,cstride=1,facecolors=get_cmap("jet")(Al))
    # colorbar(shrink=0.7)
    title("Technological progress")

    ax = fig.add_subplot(236, projection="3d")
    surf(xGrid,yGrid,Utility,rstride=1,cstride=1,facecolors=get_cmap("jet")(Utility))
    # colorbar(shrink=0.7)
    title("Individual utility")

    # tight_layout()
    if !isdir(p.folder_name)
        mkdir(p.folder_name)
    end
    pause(0.01)

    if saveFig
        savefig("$(p.folder_name)/RUMSurf$(lpad(string(round(t,digits=1)),5,"0")).png")
    end
end


function plotAll(sol,p; t = p.T_end, saveFig = false)
    # fig = figure(figsize=(15,10/(p.Nx/p.Ny)))
    # fig.clear()
    
    subplots = [231,232,233,234,235,236]

    if t > 0
        for plt in subplots
            subplot(plt)
            fig = gca()
            lines = fig.get_lines()
            for line in lines
                line.remove()
            end
        end
    end

    l = sol(t)
    Wₕᴾl = imfilter(l,p.WₕᴾM,Fill(0,l))           # Wₕᴾ * l 
    Wₕᴾl[Wₕᴾl .< p.ϵTol] .= 0

    Al = p.GM .* Wₕᴾl                                           # local technological progress
    w =  Al .* fw.(abs.(l),p.ϵY,p.β)                            # wages
    Y =  Al .* fY.(abs.(l),p.ϵY,p.β)                            # production
    AEN =  p.GM .* ((p.τ^p.φ) * fY.(Y,p.ϵY,p.φ) .- p.γA * l)    # endogenous amenities 
    Utility = p.γw * w + p.γEN * AEN

    # suptitle("t = $(round(t,digits=3))")

    subplot(231)
    im = imshow((l)',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    # colorbar(shrink=0.7)
    title("Workers")

    subplot(232)
    imshow(w',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    # colorbar(shrink=0.7)
    title("Wages")

    subplot(233)
    imshow((Y)',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    # colorbar(shrink=0.7)
    title("Total income")

    subplot(234)
    imshow((AEN)',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    # colorbar(shrink=0.7)
    title("Endogenous amenities")

    subplot(235)
    imshow((Al)',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    # colorbar(shrink=0.7)
    title("Technological progress")

    subplot(236)
    imshow((Utility)',origin="lower",extent=(0.0,p.Lx,0.0,p.Ly))
    # colorbar(shrink=0.7)
    title("Individual utility")

    # tight_layout()
    if !isdir(p.folder_name)
        mkdir(p.folder_name)
    end
    pause(0.01)

    if saveFig
        savefig("$(p.folder_name)/RUM$(lpad(string(round(t,digits=1)),5,"0")).png")
    end
end

function saveAllPaper()
    # @load "Lx=4,Ly=4,sigma005,h039/solp.jld2"
    
    global fig = figure(figsize=(15,max(10/(p.Nx/p.Ny),3.0)))
    plotAll(sol,p,t=185; saveFig=true)
    tight_layout()
    savefig("allt=185h=03.eps")
    # close()

    # plotFunction(sol(1),p;surfPlot =true)
    # savefig("single0.eps")
    # close()
    # plotFunction(sol(0.5),p;surfPlot =true)
    # savefig("single05.eps")
    # close()
    # plotFunction(sol(1),p;surfPlot =true)
    # savefig("single1.eps")
    # close()
    # plotFunction(sol(5),p;surfPlot =true)
    # savefig("single5.eps")
    # close()
    # plotFunction(sol(10),p;surfPlot =true)
    # savefig("single10.eps")
    # close()
    # plotFunction(sol(20),p;surfPlot =true)
    # savefig("single20.eps")
    # close()

    plotFunction(sol(181),p;surfPlot =true)
    savefig("t=181h=0.3.eps")
    close()
    plotFunction(sol(182),p;surfPlot =true)
    savefig("t=182h=0.3.eps")
    close()
    plotFunction(sol(183),p;surfPlot =true)
    savefig("t=183h=0.3.eps")
    close()
    plotFunction(sol(184),p;surfPlot =true)
    savefig("t=184h=0.3.eps")
    close()
    plotFunction(sol(185),p;surfPlot =true)
    savefig("t=185h=0.3.eps")
    close()

    plotFunction(sol(0),p;surfPlot =true)
    savefig("t=0ViaEmilia.eps")
    close()
    plotFunction(sol(10),p;surfPlot =true)
    savefig("t=10ViaEmilia.eps")
    close()
    plotFunction(sol(20),p;surfPlot =true)
    savefig("t=20ViaEmilia.eps")
    close()
    plotFunction(sol(100),p;surfPlot =true)
    savefig("t=100ViaEmilia.eps")
    close()

end






