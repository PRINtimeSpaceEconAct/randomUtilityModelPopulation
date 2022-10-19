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
        if integrator.p.show 
            plotAll(integrator.sol,p,t=integrator.t, saveFig = true)
        end
        @show integrator.t 
    end
    cb = PeriodicCallback(affect!,T_plot; initial_affect=true)

    if p.show 
        fig = figure(figsize=(15,max(10/(p.Nx/p.Ny),3.0)))
    end
    sol = solve(prob,Tsit5(),
        callback=cb)
        # saveat=t_save,
        # progress=true,progress_steps = 1,
        # isoutofdomain = (u,p,t) -> any(x->x<0, u),

    if !isdir(p.folder_name)
        mkdir(p.folder_name)
    end
    print("saving...\n")
    @save "$(p.folder_name)/solp.jld2" sol p
    print("done\n")

    return sol,p

end

