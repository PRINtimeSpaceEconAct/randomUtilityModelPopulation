
@with_kw struct RUM{T}

    # domain
    Lx::T = 6.0
    Ly::T = 6.0
    T_end::T = 120.0
    borderLength::T = 0.2

    # numerical
    ϵTol::T = 1e-12                # threshold to consider zero 
    Δx::T = 1e-2
    Nx::Int = Int(Lx/Δx)
    Ny::Int = Int(Ly/Δx)
    x::LinRange{T, Int64} = LinRange(0,Lx,Nx)
    y::LinRange{T, Int64} = LinRange(0,Ly,Ny)

    # global parameters
    mass::T = π * 1.0   # total mass
    n::T = 0.0          # birth rate
    σ::T = 0.05          # noise
    β::T = 0.6          # production function 

    # local technological progress Al
    hₚ::T = 0.3           # bandwith convolution
    hₛ::T = 0.2           # bandwith sharp mollifier
    WₕᴾM::Matrix{T} = make_WDiscrete(Δx,hₚ)            
    WₕᴾS::Matrix{T} = make_WDiscrete(Δx,hₛ)            
    GM::Matrix{T} = make_smoothedBump(x,y,Δx,borderLength+2*hₛ,ones(1,1),0,ϵTol) 
    # GM::Matrix{T} = make_GCross(x,y,Δx,1.0,WₕᴾS,ϵTol)
    @assert hₚ < borderLength + hₛ

    # wages   
    ϵY::T = 0.01        # minimum threshold up to which l^β becomes linear 
    γw::T = 0.01        # speed 
    
    # endogenous amenities
    γEN::T = 0.06       # speed
    τ::T = 0.2          # share income to amenities
    φ::T = 0.5          # production function amenities
    γA::T = 0.2         # intensity congestion 

    # exogenous amenities
    γES::T = 10.0                                                        # speed 
    AES::Matrix{T} =  make_smoothedBump(x,y,Δx,borderLength,WₕᴾS,hₛ,ϵTol)  # spatial distribution
    ∂xAES::Matrix{T} = ∂x(AES,Nx,Ny,Δx)                                 # precompute ∂x
    ∂yAES::Matrix{T} = ∂y(AES,Nx,Ny,Δx)                                 # precompute ∂y

    # initial condition
    # u₀::Matrix{T} = make_u₀Flat(x,y,Δx,borderLength+3.25*hₛ,WₕᴾM,hₚ,mass,ϵTol)
    u₀::Matrix{T} = make_u₀Flat(x,y,Δx,[Lx/3+0.5,Ly/2],[1.6,2.],WₕᴾM,hₚ,mass,ϵTol)
    # u₀::Matrix{T} = make_u₀Circle(x,y,Δx,[Lx/2,Ly/2],1.0,WₕᴾM,hₚ,borderLength,mass)
    # u₀::Matrix{T} = make_u₀Gauss(x,y,[Lx/2,Ly/2],0.9)
    # u₀::Matrix{T} = make_u₀Cross(x,y,Δx,1.0,WₕᴾS,ϵTol,1.0)
    # u₀::Matrix{T} = make_u₀GaussCross(x,y,Δx,1.0,WₕᴾS,ϵTol,1.0,[Lx/2,Ly/2],1.0,mass)

    # saving 
    folder_name::String = "Lx=6,Ly=6,GUnif,u0FakeViaEmilia" 
    show::Bool = true
    

end

function Wᴾ(x)
    """ smoothing kernel """
    return norm(x) <= 1 ? 1-norm(x) : 0.0
end

function Wᴾ(x,h)
    """ rescaled kernel """
    return 1/h*Wᴾ(x/h)
end

function fY(x,p) return fY(x,p.ϵY,p.β)  end

function fY(x,ϵY,β)
    if x > ϵY 
        return x^β + (-1/2*β^2 + 3/2*β - 1)*ϵY^β
    else 
        return 1/2*β*(β-1)*ϵY^(β-2)*x^2 + β*(2-β)*ϵY^(β-1)*x
    end
end

function fw(x,p) return fw(x,p.ϵY,p.β) end

function fw(x,ϵY,β)
    if x > ϵY 
        return x^(β-1)
    else 
        return (β*(β-1)*ϵY^(β-2)*x + β*(2-β)*ϵY^(β-1))/β
    end
end

function make_WDiscrete(Δx,bandwith)
    Npt = ceil(Int,bandwith/Δx)
    Wd = zeros(2Npt+1,2Npt+1)
    x = -Npt*Δx:Δx:Npt*Δx
    y = -Npt*Δx:Δx:Npt*Δx
    for i in eachindex(x), j in eachindex(y)
        Wd[i,j] = Wᴾ([x[i],y[j]],bandwith)
    end
    return normalize(Wd,1)
end

function make_smoothedBump(x,y,Δx,borderLength,Wd,bandwidth,ϵTol) 
    Nx = length(x)
    Ny = length(y)
    A = bump(borderLength+bandwidth,1.0,Nx,Ny,Δx)

    A = imfilter(A,Wd,Fill(0,A))
    A[A .< ϵTol] .= 0
    return A
end

function make_GCross(x,y,Δx,crossWidth,Wd,ϵTol)
    Nx = length(x)
    Ny = length(y)

    crossWidthNpt = floor(Int,crossWidth/Δx)
    crossNxLeft = round(Int,(Nx-crossWidthNpt)/2)
    crossNyLeft = round(Int,(Ny-crossWidthNpt)/2)
    G = ones(Nx,Ny)
    G[1:crossNxLeft,1:crossNyLeft] .= 0
    G[crossNxLeft+crossWidthNpt:Nx,1:crossNyLeft] .= 0
    G[1:crossNxLeft,crossNyLeft+crossWidthNpt:Ny] .= 0
    G[crossNxLeft+crossWidthNpt:Nx,crossNyLeft+crossWidthNpt:Ny] .= 0

    G = imfilter(G,Wd,Fill(0,G))
    G[G .< ϵTol] .= 0

    return G
end

function make_u₀Cross(x,y,Δx,crossWidth,Wd,ϵTol,borderLength)
    Nx = length(x)
    Ny = length(y)
    rNpts = ceil(Int,borderLength/Δx)

    crossWidthNpt = floor(Int,crossWidth/Δx)
    crossNxLeft = round(Int,(Nx-crossWidthNpt)/2)
    crossNyLeft = round(Int,(Ny-crossWidthNpt)/2)
    G = ones(Nx,Ny)
    G[1:crossNxLeft,1:crossNyLeft] .= 0
    G[crossNxLeft+crossWidthNpt:Nx,1:crossNyLeft] .= 0
    G[1:crossNxLeft,crossNyLeft+crossWidthNpt:Ny] .= 0
    G[crossNxLeft+crossWidthNpt:Nx,crossNyLeft+crossWidthNpt:Ny] .= 0
    G[1:rNpts,:] .= 0
    G[:,1:rNpts] .= 0
    G[(Nx-rNpts):Nx,:] .= 0
    G[:,(Ny-rNpts):Ny] .= 0

    G = imfilter(G,Wd,Fill(0,G))
    G[G .< ϵTol] .= 0

    return G
end

function make_u₀Flat(x,y,Δx,borderLength,Wd,bandwidth,mass,ϵTol)
    Nx = length(x)
    Ny = length(y)
    u₀ = bump(borderLength+bandwidth,1,Nx,Ny,Δx)

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀ = u₀ / (sum(u₀)*Δx^2) * mass
    u₀[u₀ .< ϵTol] .= 0

    return u₀
end

function make_u₀Flat(x,y,Δx,center,r,Wd,bandwidth,mass,ϵTol)
    Nx = length(x)
    Ny = length(y)
    u₀ = bump(center,r,mass,Nx,Ny,Δx)

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀ = u₀ / (sum(u₀)*Δx^2) * mass
    u₀[u₀ .< ϵTol] .= 0

    return u₀
end

function make_u₀Gauss(x,y,center,sd)
    X = MvNormal(center,sd^2*I(2))
    Nx = length(x)
    Ny = length(y)
    u₀ = reshape(pdf(X,[[xi,yi] for xi in x for yi in y]),Nx,Ny)
    return u₀
end

function make_u₀GaussCross(x,y,Δx,crossWidth,Wd,ϵTol,borderLength,center,sd,mass)
    X = MvNormal(center,sd^2*I(2))
    Nx = length(x)
    Ny = length(y)
    u₀ = reshape(pdf(X,[[xi,yi] for xi in x for yi in y]),Nx,Ny)

    rNpts = ceil(Int,borderLength/Δx)

    crossWidthNpt = floor(Int,crossWidth/Δx)
    crossNxLeft = round(Int,(Nx-crossWidthNpt)/2)
    crossNyLeft = round(Int,(Ny-crossWidthNpt)/2)
    u₀[1:crossNxLeft,1:crossNyLeft] .= 0
    u₀[crossNxLeft+crossWidthNpt:Nx,1:crossNyLeft] .= 0
    u₀[1:crossNxLeft,crossNyLeft+crossWidthNpt:Ny] .= 0
    u₀[crossNxLeft+crossWidthNpt:Nx,crossNyLeft+crossWidthNpt:Ny] .= 0
    u₀[1:rNpts,:] .= 0
    u₀[:,1:rNpts] .= 0
    u₀[(Nx-rNpts):Nx,:] .= 0
    u₀[:,(Ny-rNpts):Ny] .= 0

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀[u₀ .< ϵTol] .= 0
    u₀ = u₀ * mass

    return u₀

end

function make_u₀Circle(x,y,Δx,center,radius,Wd,ϵTol,borderLength,mass)
    Nx = length(x)
    Ny = length(y)
    u₀ = zeros(Nx,Ny)

    rNpts = ceil(Int,borderLength/Δx)
    for i=1:Nx, j=1:Ny
        if (x[i]-center[1])^2 + (y[j]-center[2])^2 ≤ radius^2
            u₀[i,j] = 1
        end
    end

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀[u₀ .< ϵTol] .= 0
    u₀ = u₀ * mass

    return u₀
end

function bump(borderLength::T,height,Nx,Ny,Δx) where T
    cube = zeros(Nx,Ny)
    rNpts = ceil(Int,borderLength/Δx)
    cube[rNpts:(Nx-rNpts+1),rNpts:(Ny-rNpts+1)] .= height

    return cube
end

function bump(center::Vector{T},r::Vector{T},height,Nx,Ny,Δx) where T
    rx, ry = r
    rNptsX = ceil(Int,rx/Δx)
    rNptsY = ceil(Int,ry/Δx)
    cᵢ,cⱼ = round.(Int,center/Δx)
    
    cube = zeros(Nx,Ny)
    cube[cᵢ-rNptsX:cᵢ+rNptsX+1,cⱼ-rNptsY:cⱼ+rNptsY+1] .= height
    
    return cube
end

