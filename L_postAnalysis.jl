

function totalUtility(sol,p)


    times = 0.0:1.0:p.T_end
    utilityVec = zeros(length(times))
    @showprogress for (i,t) in enumerate(times)
        l = sol(t)
        Wₕᴾl = imfilter(l,p.WₕᴾM,Fill(0,l))           # Wₕᴾ * l 
        Wₕᴾl[Wₕᴾl .< p.ϵTol] .= 0
    
        Al = p.GM .* Wₕᴾl                                           # local technological progress
        w =  Al .* fw.(abs.(l),p.ϵY,p.β)                            # wages
        Y =  Al .* fY.(abs.(l),p.ϵY,p.β)                            # production
        AEN =  p.GM .* ((p.τ^p.φ) * fY.(Y,p.ϵY,p.φ) .- p.γA * l)    # endogenous amenities 
        Utility = p.γw * w + p.γEN * AEN
        utilityVec[i] = sum(Utility .* l) * p.Δx^2
    end
    plot(times,utilityVec)
    title("Total utility")
    grid()
    savefig("totalUtility" * p.folder_name * ".png")
end


# @load "Lx=6,Ly=6,GCross,u0GaussCross/solp.jld2"
# @load "Lx=4,Ly=4,sigma005,T240(metastability)/solp.jld2" sol
