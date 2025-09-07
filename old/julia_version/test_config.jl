using Pkg; Pkg.activate(@__DIR__); include("RBFNS.jl"); using .RBFNS; cfg = default_config(); println("Reynolds number: $(cfg.reynolds_number)")
