
@with_kw struct RUM{T}

    # domain
    Lx::T = 2.0
    Ly::T = 1.0
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
    n::T = 0.0          # birth rate
    σ::T = 0.1          # noise
    β::T = 0.4          # production function 

    # local technological progress Al
    hₚ::T = 0.15           # bandwith convolution
    hₛ::T = 0.1            # bandwith sharp mollifier
    WₕᴾM::Matrix{T} = make_WDiscrete(Δx,hₚ)
    WₕᴾS::Matrix{T} = make_WDiscrete(Δx,hₛ)            
    GM::Matrix{T} = make_smoothedBump(x,y,Δx,borderLength,WₕᴾS,ϵTol) 
    @assert hₚ < borderLength + hₛ

    # wages   
    ϵY::T = 0.01        # minimum threshold up to which l^β becomes linear 
    γw::T = 0.01        # speed 
    # γw::T = 0.0        # speed 

    
    # endogenous amenities
    γEN::T = 0.06       # speed
    # γEN::T = 0.00       # speed
    τ::T = 0.2          # share income to amenities
    φ::T = 0.5          # production function amenities
    γA::T = 0.2         # intensity congestion 

    # exogenous amenities
    γES::T = 10.0                                                        # speed 
    AES::Matrix{T} =  make_smoothedBump(x,y,Δx,borderLength,WₕᴾS,ϵTol)  # spatial distribution
    ∂xAES::Matrix{T} = ∂x(AES,Nx,Ny,Δx)                                 # precompute ∂x
    ∂yAES::Matrix{T} = ∂y(AES,Nx,Ny,Δx)                                 # precompute ∂y

    # initial condition
    u₀::Matrix{T} = make_u₀(x,y,Δx,borderLength+hₛ,WₕᴾM,ϵTol)

    # saving 
    folder_name::String = "pics" 
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

function make_smoothedBump(x,y,Δx,borderLength,Wd,ϵTol) 
    Nx = length(x)
    Ny = length(y)
    borderNpts = ceil(Int,borderLength/Δx)
    
    A = ones(Nx,Ny)
    A[1:borderNpts,:] .= 0.0
    A[:,1:borderNpts] .= 0.0
    A[end-borderNpts+1:end,:] .= 0.0
    A[:,end-borderNpts+1:end] .= 0.0

    A = imfilter(A,Wd,Fill(0,A))
    A[A .< ϵTol] .= 0
    return A
end

function make_u₀(x,y,Δx,borderLength,Wd,ϵTol)
    Nx = length(x)
    Ny = length(y)
    u₀ = bump(borderLength,1,Nx,Ny,Δx)

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀ = u₀ / (sum(u₀)*Δx^2)
    u₀[u₀ .< ϵTol] .= 0

    return u₀
end

function make_u₀(x,y,Δx,center::Vector{T},r::Vector{T},Wd,ϵTol) where T
    Nx = length(x)
    Ny = length(y)
    u₀ = bump(center,r,1,Nx,Ny,Δx)

    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀ = u₀ / (sum(u₀)*Δx^2)
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

