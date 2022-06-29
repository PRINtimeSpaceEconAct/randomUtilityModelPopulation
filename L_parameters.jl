
@with_kw struct RUM{T}

    # domain
    Lx::T = 1.0
    Ly::T = 1.0
    T_end::T = 1.0
    borderLength::T = 0.2

    # numerical
    Δx::T = 1e-2
    Nx::Int = Int(Lx/Δx)
    Ny::Int = Int(Ly/Δx)
    x::LinRange{T, Int64} = LinRange(0,Lx,Nx)
    y::LinRange{T, Int64} = LinRange(0,Ly,Ny)


    # global parameters
    n::T = 0.0          # birth rate
    σ::T = 0.01          # noise
    β::T = 0.5          # production function 

    # local technological progress Al
    hₚ::T = 0.15         # bandwith 
    WₕᴾM::Matrix{T} = make_WDiscrete(Δx,hₚ)
    GM::Matrix{T} = G(Nx,Ny,Δx,borderLength,WₕᴾM) 
    @assert hₚ < borderLength

    # wages   
    γw::T = 0.01        # speed 
    
    # endogenous amenities
    γEN::T = 0.06       # speed
    τ::T = 0.2          # share income to amenities
    φ::T = 0.5          # production function amenities
    γA::T = 0.1         # intensity congestion 

    # exogenous amenities
    γES::T = 1.0                                    # speed 
    AES::Matrix{T} = make_Vflat(Nx,Ny,Δx,hₚ/2,WₕᴾM)  # spatial distribution
    ∂xAES::Matrix{T} = ∂x(AES,Nx,Ny,Δx)             # precompute ∂x
    ∂yAES::Matrix{T} = ∂y(AES,Nx,Ny,Δx)             # precompute ∂y


    # initial condition
    u₀::Matrix{T} = make_u₀(Nx,Ny,Δx,Lx,Ly,WₕᴾM)


end

function make_u₀(Nx,Ny,Δx,Lx,Ly,Wd)
    # u₀ = zeros(Nx,Ny)
    # u₀[Int(Nx/2):Int(Nx/2)+10,Int(Nx/2):Int(Nx/2)+10] .= 1
    
    # u₀ = bump([L/4,L/4],0.3,0.5,Nx,Ny,Δx) +   
    #     bump([3/4*L,L/4],0.3,0.5,Nx,Ny,Δx) +
    #     bump([L/2,3/4*L],0.3,1.0,Nx,Ny,Δx)   
    # u₀ = bump([Lx/2,Ly/2],[3.0,0.4],1.0,Nx,Ny,Δx) 

    u₀ = bump([Lx/2,Ly/2],[Lx/8,Ly/8],1,Nx,Ny,Δx)
    u₀ = imfilter(u₀,Wd,Fill(0,u₀))
    u₀ = u₀ / (sum(u₀)*Δx^2)
    return u₀
end

function bump(center,r,height,Nx,Ny,Δx)
    rx, ry = r
    rNptsX = ceil(Int,rx/Δx)
    rNptsY = ceil(Int,ry/Δx)
    cᵢ,cⱼ = round.(Int,center/Δx)
    
    cube = zeros(Nx,Ny)
    #symmetric
    cube[cᵢ-rNptsX:cᵢ+rNptsX+1,cⱼ-rNptsY:cⱼ+rNptsY+1] .= height
    #asymmetric
    # cube[cᵢ-rNptsX:cᵢ+rNptsX,cⱼ-rNptsY:cⱼ+rNptsY+1] .= height
    
    return cube
end

function make_Vflat(Nx,Ny,Δx,borderLength,Wd)
    V = ones(Nx,Ny)
    borderNpts = ceil(Int,borderLength/Δx)
    V[1:borderNpts,:] .= 0.0
    V[:,1:borderNpts] .= 0.0
    V[end-borderNpts+1:end,:] .= 0.0
    V[:,end-borderNpts+1:end] .= 0.0

    V = imfilter(V,Wd,Fill(0,V))
    return V
end

function G(Nx,Ny,Δx,borderLength,Wd) 
    GM = ones(Nx,Ny)
    borderNpts = ceil(Int,borderLength/Δx) + 10
    GM[1:borderNpts,:] .= 0.0
    GM[:,1:borderNpts] .= 0.0
    GM[end-borderNpts+1:end,:] .= 0.0
    GM[:,end-borderNpts+1:end] .= 0.0
    GM = imfilter(GM,Wd,Fill(0,GM))
    return GM
end

function Wᴾ(x)
    """ smoothing kernel """
    return norm(x) <= 1 ? 1-norm(x) : 0.0
end

function Wᴾ(x,h)
    """ rescaled kernel """
    return 1/h*Wᴾ(x/h)
end

function make_WDiscrete(Δx,bandwith)
    Npt = ceil(Int,bandwith/Δx)
    Wd = zeros(2Npt+1,2Npt+1)
    x = -Npt*Δx:Δx:Npt*Δx
    y = -Npt*Δx:Δx:Npt*Δx
    for i in 1:length(x), j in 1:length(y)
        Wd[i,j] = Wᴾ([x[i],y[j]],bandwith)
    end
    return normalize(Wd,1)
end

