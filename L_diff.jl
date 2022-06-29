
function df!(dl,l,p,t)
    
    Wₕᴾl = imfilter(l,p.WₕᴾM,Fill(0,l))           # Wₕᴾ * l 
    Al = p.GM .* Wₕᴾl                             # local technological progress
    w =  Al .* abs.(l).^(p.β-1) .* (l .> 0 )      # wages
    Y =  Al .* abs.(l).^p.β                       # production
    AEN =  (p.τ * abs.(Y)).^p.φ .- p.γA * l       # endogenous amenities 
    
    dl .= ( p.σ^2/2 * Δ(l,p) 
        - p.γw * (∂x(l .* ∂x(w,p),p) .+  ∂y(l .* ∂y(w,p),p))
        - p.γEN * (∂x(l .* ∂x(AEN,p),p) .+  ∂y(l .* ∂y(AEN,p),p))
        - p.γES * (∂x(l .* p.∂xAES,p) .+  ∂y(l .* p.∂yAES,p)) 
        + p.n * l )
    
    # @show t

    # if sum(isfinite.(dl) .== 0) != 0
    #     @infiltrate 
        # error("error in dl")
    # end

    # if sum(l .< 0) != 0 
    #     @show sum(l .< 0)
    # end
end

function Δ(u,p) return Δ(u,p.Nx,p.Ny,p.Δx) end

function Δ(u,Nx,Ny,Δx)

    Δu = similar(u)
    @tturbo for i in 2:Nx-1, j in 2:Ny-1
        Δu[i,j] = 1/(Δx)^2 * (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])
    end

    # dirichlet 0
    @tturbo for i in 2:Nx-1   
        Δu[i,1] = 1/(Δx)^2 * (u[i-1,1] + u[i+1,1] + u[i,2] - 4*u[i,1])
        Δu[i,Ny] = 1/(Δx)^2 * (u[i-1,Ny] + u[i+1,Ny] + u[i,Ny-1] - 4*u[i,Ny])
    end
    @tturbo for j = 2:Ny-1
        Δu[1,j] = 1/(Δx)^2 * (u[2,j] + u[1,j-1] + u[1,j+1] - 4*u[1,j])
        Δu[Nx,j] = 1/(Δx)^2 * (u[Nx-1,j] + u[Nx,j-1] + u[Nx,j+1] - 4*u[Nx,j])
    end
    Δu[1,1] = 1/(Δx)^2 * (u[2,1] + u[1,2] - 4*u[1,1])
    Δu[1,Ny] = 1/(Δx)^2 * (u[2,Ny] + u[1,Ny-1] - 4*u[1,Ny])
    Δu[Nx,1] = 1/(Δx)^2 * (u[Nx-1,1] + u[Nx,2] - 4*u[Nx,1])
    Δu[Nx,Ny] = 1/(Δx)^2 * (u[Nx-1,Ny] + u[Nx,Ny-1] - 4*u[Nx,Ny])

    return Δu
end

function ∂x(u,p) return ∂x(u,p.Nx,p.Ny,p.Δx) end

function ∂x(u,Nx,Ny,Δx)

    ∂x = similar(u)
    @tturbo for i in 1:Nx, j in 2:Ny-1
        ∂x[i,j] = 1/(2*Δx) * (u[i,j+1] - u[i,j-1])
    end

    # dirichlet 0
    @tturbo for i in 1:Nx
        ∂x[i,1] =  1/(2*Δx) * u[i,2] 
        ∂x[i,Ny] = - 1/(2*Δx) * u[i,Ny-1]
    end
    return ∂x
end

function ∂y(u,p) return ∂y(u,p.Nx,p.Ny,p.Δx) end

function ∂y(u,Nx,Ny,Δx)

    ∂y = similar(u)
    @tturbo for i in 2:Nx-1, j in 1:Ny
        ∂y[i,j] = 1/(2*Δx) * (u[i-1,j] - u[i+1,j])
    end

    # dirichlet 0
    @tturbo for j in 1:Ny
        ∂y[1,j] =  - 1/(2*Δx) * u[2,j] 
        ∂y[Nx,j] = 1/(2*Δx) * u[Nx-1,j]
    end
    return ∂y
end




