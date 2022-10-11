

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

function plotSol(sol,p; t_step = 0.1)
    
    t_plot = sol.t[1]:t_step:sol.t[end]

    fig = figure(figsize=(15,10/(p.Nx/p.Ny)))
    T_plot = t_step 
    for t in t_plot
        plotAll(sol,p; t = t)
        if !isdir(p.folder_name)
            mkdir(p.folder_name)
        end
        savefig("$(p.folder_name)/RUM$(lpad(string(round(t,digits=1)),4,"0")).png")
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

    Al = p.GM .* Wₕᴾl                                   # local technological progress
    w =  Al .* fw.(abs.(l),p.ϵY,p.β)                    # wages
    Y =  Al .* fY.(abs.(l),p.ϵY,p.β)                    # production
    AEN =  (p.τ^p.φ) * fY.(Y,p.ϵY,p.φ) .- p.γA * l      # endogenous amenities 
    Utility = p.γw * w + p.γEN * AEN

    suptitle("t = $(round(t,digits=3))")

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

    tight_layout()
    if !isdir(p.folder_name)
        mkdir(p.folder_name)
    end
    savefig("$(p.folder_name)/RUM$(lpad(string(round(t,digits=1)),4,"0")).png")
    pause(0.01)
end

function histSignificative(l,Img,p; lCut = 1e-2, ICut = 0.0, nb = 20, titleStr=" ", pesi=false)
    # l = matrix of worker distribution
    # I = matrix of distribution to hist
    # lCut = lower level where to cut l
    # ICut = lower level where to cut further Img

    i = findall(l .> lCut)
    print("keeping $(sum(l[i])*p.Δx^2) total mass\n")
    l = l[i]
    Img = Img[i]

    i = findall(Img .> ICut)
    l = l[i]
    Img = Img[i]

    if pesi == true
        minl = minimum(l)
        weigh = round.(Int,l/minl)
        ImgTmp = [Img[i]*ones(weigh[i]) for i in 1:length(weigh)]
        Img = vcat(ImgTmp...)
    end
    
    
    hist(Img,density=true,nb)
    density = kde(Img)
    i = findall((density.x .> minimum(Img)) .& (density.x .< maximum(Img)))
    xval = density.x[i]
    val = density.density[i]
    plot(xval,val)

    
    title(titleStr)
    tight_layout() 
end

function spatMean(l,p)
    l1 = vec(sum(l,dims=2)*p.Δx)
    m1 = sum(l1 .* p.x) * p.Δx

    l2 = vec(sum(l,dims=1)*p.Δx)
    m2 = sum(l2 .* p.y) * p.Δx

    return [m1,m2]
end

function moment(l,p,order)
    n = order
    m = spatMean(l,p)

    l1 = vec(sum(l,dims=2)*p.Δx)
    l2 = vec(sum(l,dims=1)*p.Δx)

    mn1 = sum(l1 .* (p.x .- m[1]).^n ) * p.Δx
    mn2 = sum(l2 .* (p.y .- m[2]).^n ) * p.Δx

    return [mn1,mn2]
end


function plotPane(sol,p; t_save = [0.01,0.1,1.0,20.0,40.0],preStr="")
    # plotPane(sol,p,preStr="beta05")
    ioff()
    for t in t_save
        plotAll(sol,p; t = t)
        savefig(preStr * "t" * string(t) * ".pdf")
    end
end

function plotHist(sol,p;preStr="")
    ioff()

    l = sol[end]

    Wₕᴾl = imfilter(l,p.WₕᴾM,Fill(0,l))           # Wₕᴾ * l 
    Al = p.GM .* Wₕᴾl                             # local technological progress
    w =  Al .* abs.(l).^(p.β-1) .* (l .> 0 )      # wages
    Y =  Al .* abs.(l).^p.β                       # production
    AEN =  (p.τ * abs.(Y)).^p.φ .- p.γA * l
    Utility = p.γw * w + p.γEN * AEN

    fig = figure()
    histSignificative(l,l,p, titleStr="Workers")
    savefig(preStr * "L.pdf")
    close(fig)

    fig = figure()
    histSignificative(l,AEN,p, titleStr="Endogenous amenities")
    savefig(preStr * "AEN.pdf")
    close(fig)

    fig = figure()
    histSignificative(l,Al,p, titleStr="Technological progress")
    savefig(preStr * "Al.pdf")
    close(fig)

    fig = figure()
    histSignificative(l,Y,p, titleStr="Total income")
    savefig(preStr * "totalIncome.pdf")
    close(fig)

    fig = figure()
    histSignificative(l,w,p, titleStr="Wages"; pesi = true)
    savefig(preStr * "wages.pdf")
    close(fig)

    fig = figure()
    histSignificative(l,Utility,p, ICut = -Inf, titleStr="Utility"; pesi = true)
    savefig(preStr * "utility.pdf")
    close(fig)

end

function relevantValuesSummary(sol,p; lCut = 1e-2)
    l = sol[end]
    
    Wₕᴾl = imfilter(l,p.WₕᴾM,Fill(0,l))           # Wₕᴾ * l 
    Al = p.GM .* Wₕᴾl                             # local technological progress
    w =  Al .* abs.(l).^(p.β-1) .* (l .> 0 )      # wages
    Y =  Al .* abs.(l).^p.β                       # production
    AEN =  (p.τ * abs.(Y)).^p.φ .- p.γA * l
    Utility = p.γw * w + p.γEN * AEN

    i = findall(l .> lCut)
    println("keeping $(sum(l[i])*p.Δx^2) total mass")
    Area = length(i)*p.Δx^2
    l = l[i]
    Al = Al[i]
    w = w[i]
    Y = Y[i]
    AEN = AEN[i]
    Utility = Utility[i]

    println("Area = $(round(Area,digits=3))")
    println("Densità = $(round(1/Area,digits=3))")

    println("mean w = $(round(mean(w,StatsBase.weights(l)),digits=3))")
    println("sum Y = $(round(sum(Y),digits=3))")
    println("mean Al = $(round(mean(Al),digits=3))")
    println("mean AEN = $(round(mean(AEN),digits=3))")
    println("mean Utility = $(round(mean(Utility,StatsBase.weights(l)),digits=3))")

    println("std w = $(round(std(w,StatsBase.weights(l)),digits=3))")
    println("std Y = $(round(std(Y),digits=3))")
    println("std Al = $(round(std(Al),digits=3))")
    println("std AEN = $(round(std(AEN),digits=3))")
    println("std Utility = $(round(std(Utility,StatsBase.weights(l)),digits=3))")

    println("skewness w = $(round(StatsBase.skewness(w,StatsBase.weights(l)),digits=3))")
    println("skewness Y = $(round(StatsBase.skewness(Y),digits=3))")
    println("skewness Al = $(round(StatsBase.skewness(Al),digits=3))")
    println("skewness AEN = $(round(StatsBase.skewness(AEN),digits=3))")
    println("skewness Utility = $(round(StatsBase.skewness(Utility,StatsBase.weights(l)),digits=3))")

    println("kurtosis w = $(round(StatsBase.kurtosis(w,StatsBase.weights(l)),digits=3))")
    println("kurtosis Y = $(round(StatsBase.kurtosis(Y),digits=3))")
    println("kurtosis Al = $(round(StatsBase.kurtosis(Al),digits=3))")
    println("kurtosis AEN = $(round(StatsBase.kurtosis(AEN),digits=3))")
    println("kurtosis Utility = $(round(StatsBase.kurtosis(Utility,StatsBase.weights(l)),digits=3))")

    println("90/100 w = $(round(quantile(w,StatsBase.weights(l),0.9)/quantile(w,StatsBase.weights(l),0.1),digits=3))")
    println("90/100 Y = $(round(quantile(Y,0.9)/quantile(Y,0.1),digits=3))")
    println("90/100 Al = $(round(quantile(Al,0.9)/quantile(Al,0.1),digits=3))")
    println("90/100 AEN = $(round(quantile(AEN,0.9)/quantile(AEN,0.1),digits=3))")
    println("90/100 Utility = $(round(quantile(Utility,StatsBase.weights(l),0.9)/quantile(Utility,StatsBase.weights(l),0.1),digits=3))")
    
end



