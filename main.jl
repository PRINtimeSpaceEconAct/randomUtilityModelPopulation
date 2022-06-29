# Started on April 6th 2022
# Cristiano Ricci - cristiano.ricci6@gmail.com

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

include("L_LoadAll.jl")


function main(p)
    
    T_span = (0.0,p.T_end);
    prob = ODEProblem(df!, p.u₀, T_span, p);
    sol = solve(prob,ROCK2(),progress=true,progress_steps = 1);
    return sol,p
end

