using DifferentialEquations, DiffEqBiological
using Gadfly, RDatasets
using Latexify

# iris = dataset("datasets", "iris")

# Model1 without DM cycle 
Oct4_auto_active = @reaction_network begin
  (a₀,d₀),       O + O ↔ O₂         # OCT4 Dimerization
  (r₁,r₀ ),      O₂ + D ↔ Dₒ        # OCT4 Auto-activation
   α₀,           D → m + D          # OCT4 Transcription (leaky)
   α₁,           Dₒ → m + Dₒ        # OCT4 Transcription (activ.)
   k,            m → m + O          # OCT4 Translation
  (β,β),         (O₂, O)→∅          # OCT4 Dilution
  (δₘ,δₒ),       (m,O) → ∅          # mRNA/Protein Degradation
    end a₀ d₀ r₁ r₀ α₀ α₁ k β
p = (0.00166,0.0001,0.1, 0.2,1,0.1)
tspan = (0., 100.)
u0 = [301., 100., 0., 0.,0.,0.]  # S = 301, E = 100, SE = 0, P = 0

# solve ODEs
oprob = ODEProblem(Oct4_auto_active, u0, tspan, p)
osol  = solve(oprob, Tsit5())
using IterableTables, DataFrames
df = DataFrame(osol')
# plot(DataFrame(osol'), x=:x1, Geom.point, Geom.line)

