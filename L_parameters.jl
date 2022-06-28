
@with_kw struct RUM{T}
    # aggregation diffusion
    γa::T = - 0.06      # aggregation (< 0)
    γd::T = 0.03       # diffusion
    # γd::T = 0.025       # diffusion
    h::T = 0.45          # bandwith  
    # h::T = 0.3          # bandwith  
    # γ₁::T = 0.01        # growth linear
    # γ₂::T = -0.02       # growth quadratic
    γ₁::T = 0.0         # growth linear
    γ₂::T = 0.0         # growth quadratic
    γᵥ::T = 1.0         # potential strength
    
    # reallocation gain
    V_GRD::T = 0.0      # reallocation gain physical capital
    V_GRA::T = 0.0      # reallocation gain labor 

    # domain
    Lx::T = 8.0
    Ly::T = 2.0
    T_end::T = 10000.0
    borderLength::T = 0.5
    @assert h < borderLength

    # numerical
    Δx::T = 1e-2
    Nx::Int = Int(Lx/Δx)
    Ny::Int = Int(Ly/Δx)
    x::LinRange{T, Int64} = LinRange(0,Lx,Nx)
    y::LinRange{T, Int64} = LinRange(0,Ly,Ny)
    Wd::Matrix{T} = make_WDiscrete(Δx,h)

    # confining potential around border
    V::Matrix{T} = make_Vflat(Nx,Ny,Δx,h/2,Wd)
    ∂xV::Matrix{T} = ∂x(V,Nx,Ny,Δx)
    ∂yV::Matrix{T} = ∂y(V,Nx,Ny,Δx)

    # initial condition
    u₀::Matrix{T} = make_u₀(Nx,Ny,Δx,Lx,Ly)

    # human capital (optional)
    β::Float64 = 1.0
    γ::Float64 = 1.0
    δ::Float64 = 1.0
    ϕ::Float64 = 1.0
    ψ::Float64 = 1.0

end

function make_u₀(Nx,Ny,Δx,Lx,Ly)
    # u₀ = zeros(Nx,Ny)
    # u₀[Int(Nx/2):Int(Nx/2)+10,Int(Nx/2):Int(Nx/2)+10] .= 1
    
    # u₀ = bump([L/4,L/4],0.3,0.5,Nx,Ny,Δx) +   
    #     bump([3/4*L,L/4],0.3,0.5,Nx,Ny,Δx) +
    #     bump([L/2,3/4*L],0.3,1.0,Nx,Ny,Δx)   
    u₀ = bump([Lx/2,Ly/2],[3.0,0.4],1.0,Nx,Ny,Δx) 
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
    V = -ones(Nx,Ny)
    borderNpts = ceil(Int,borderLength/Δx)
    V[1:borderNpts,:] .= 0.0
    V[:,1:borderNpts] .= 0.0
    V[end-borderNpts+1:end,:] .= 0.0
    V[:,end-borderNpts+1:end] .= 0.0

    V = imfilter(V,Wd,Fill(0,V))
    return V
end

function W(x)
    """ smoothing kernel """
    return norm(x) <= 1 ? 1-norm(x) : 0.0
end

function W(x,h)
    """ rescaled kernel """
    return 1/h*W(x/h)
end

function make_WDiscrete(Δx,bandwith)
    Npt = ceil(Int,bandwith/Δx)
    Wd = zeros(2Npt+1,2Npt+1)
    x = -Npt*Δx:Δx:Npt*Δx
    y = -Npt*Δx:Δx:Npt*Δx
    for i in 1:length(x), j in 1:length(y)
        Wd[i,j] = W([x[i],y[j]],bandwith)
    end
    return normalize(Wd,1)
end

