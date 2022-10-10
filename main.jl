# Started on April 6th 2022
# Cristiano Ricci - cristiano.ricci6@gmail.com

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

include("L_LoadAll.jl")

global T_plot = 0.1
function main(p)
    
    T_span = (0.0,p.T_end);
    prob = ODEProblem(df!, p.u₀, T_span, p);
    t_save = 0.0:0.1:p.T_end


    function affect!(integrator)
        # plotAll(integrator.sol,p,t=integrator.t)
        @show integrator.t 
    end
    cb = PeriodicCallback(affect!,T_plot; initial_affect=true)

    # fig = figure(figsize=(15,10/(p.Nx/p.Ny)))
    sol = solve(prob,Tsit5(),
        callback=cb)
        # saveat=t_save,
        # progress=true,progress_steps = 1,
        # isoutofdomain = (u,p,t) -> any(x->x<0, u),

    fig.close()
    return sol,p
end

