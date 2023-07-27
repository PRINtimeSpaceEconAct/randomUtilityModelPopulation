
@with_kw struct RUM{T}

    # domain
    Lx::T = 8.0
    Ly::T = 4.0
    T_end::T = 1.0
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
    # GM::Matrix{T} = make_smoothedBump(x,y,Δx,borderLength+2*hₛ,ones(1,1),0,ϵTol) 
    GM::Matrix{T} = make_GCross(x,y,Δx,1.0,WₕᴾS,ϵTol)
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
    u₀::Matrix{T} = make_u₀(x,y,Δx,borderLength+3.25*hₛ,WₕᴾM,hₚ,mass,ϵTol)

    # saving 
    folder_name::String = "Lx=8,Ly=4,GCross" 
    show::Bool = false
    

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

function make_u₀(x,y,Δx,borderLength,Wd,bandwidth,mass,ϵTol)
    Nx = length(x)
    Ny = length(y)
    u₀ = bump(borderLength+bandwidth,1,Nx,Ny,Δx)

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀ = u₀ / (sum(u₀)*Δx^2) * mass
    u₀[u₀ .< ϵTol] .= 0

    return u₀
end

function make_u₀(x,y,Δx,center::Vector{T},r::Vector{T},Wd,bandwidth,mass,ϵTol) where T
    Nx = length(x)
    Ny = length(y)
    u₀ = bump(center,r-bandwidth,1,Nx,Ny,Δx)

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀ = u₀ / (sum(u₀)*Δx^2) * mass
    u₀[u₀ .< ϵTol] .= 0

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

