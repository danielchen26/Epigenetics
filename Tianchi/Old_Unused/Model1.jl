using DifferentialEquations, DiffEqBiological
# Model1 with mRNA levels
Oct4_auto_active = @reaction_network begin
  # Oct4 cycle
  (a₀,d₀),       O + O ↔ O₂          # Oct4 Dimerization
  (r₁,r₀ ),      O₂ + D ↔ Dₒ         # Oct4 Auto-activation
   α₀,           D → m + D          # Oct4 Transcription (leaky)
   α₁,           Dₒ → m + Dₒ         # Oct4 Transcription (activ.)
   k,            m → m + O          # Oct4 Translation
  (β,β,β),      (O₂, O, m)→∅       # Oct4/mRNA Dilution
  (δₘ,δₒ),       (m,O) → ∅          # mRNA/Protein Degradation

 # Oct4 de-Methylation cycle
  (a₁,d₁),       D + DNMT ↔ C₁      # De-novo Methylation
  k₁,            C₁ → D̄ + DNMT      # De-novo Methylation
  # adk,           D → D̄           # De-novo Methylation


  (a₂,d₂),        D̄ +TET ↔ C₂        # Active de-Methylation
  k₂,             C₂ → Dₕ + TET      # Active de-Methylation
  β,              Dₕ →  D            # Passive De-Methylation

# TET protein
  γD,            ∅ → DNMT          # DNMT Constitutive Transcription
  ϕ₀,            Dₜ → Dₜ + TET       # TET Transcription/Translation (leaky)
  ϕ₁,            D♇ → D♇ + TET      # TET Transcription/Translation (activ.)
  (r₂,rᵦ),        Dₜ + O₂ ↔ D♇       # OCT4 Activating TET
  u,             ∅ → m             # Ectopic Overexpression


 # Nanog cycle
 (a₀,d₀),        N + N ↔ N₂          # Nanog Dimerization
 (r₁,r₀ ),       N₂ + Dₙ ↔ Dᴺ        # Nanog Auto-activation
  α₀,            Dₙ → mₙ + Dₙ        # Nanog Transcription (leaky)
  α₁,            Dᴺ → mₙ + Dᴺ        # Nanog Transcription (activ.)
  k,             mₙ → mₙ + N         # Nanog Translation
 (β,β,β),        (N₂, N, mₙ)→∅       # Nanog/mRNA Dilution
 (δₘ,δₒ),        (mₙ,N) → ∅          # mRNA/Protein Degradation

 # Oct4 Nanog interaction
 (ζ₁, ζ₀),       O₂ + Dₙ ↔ Dᴺ        # Oct4 dimmer activate Nanog promoter
# Nanog2 regulate Oct4?

 # Nanog + TET(1) activate OCT4
(η₁, η₀),         N + TET ↔ NT       # Nanog + TET form complex in preparation for Oct4 demethylation
 k₃,             NT + D̄ → D + NT    # Nanog-TET complex de-Methylate promoter of the Oct4

 # TET1 binding to its own promoter
 # Nanog + TET1 induce not only Oct4 but TET1
end a₀ d₀ r₁ r₀ α₀ α₁ k β δₘ δₒ a₁ d₁ k₁ a₂ d₂ k₂ γD ϕ₀ ϕ₁ r₂ rᵦ u ζ₁ ζ₀ η₁ η₀ k₃


# #        a₀    d₀     r₁     r₀     α₀   α₁    k     β       δₘ     δₒ    a₁    d₁   k₁    a₂    d₂   k₂    γD    ϕ₀    ϕ₁   r₂     rᵦ  u   ζ₁   ζ₀     η₁   η₀    k₃
# p = (100000, 1000, 100000,  1000,  0,  0.2,  0.5,  0.1,  0.069,  0.09,  100, 100, 100, 100, 100,  100, 0.01, 0.01, 0.5, 1e5, 1e3, 100 , 100,  100, 100, 100, 100 )
# tspan = (0., 1e7)
# u0 = rand(Float64,19)
#
#
# # solve ODEs
# oprob = ODEProblem(Oct4_auto_active, u0, tspan, p)
# osol  = solve(oprob, Tsit5(), maxiters = 1e6)
# plot(osol, vars = getindex(Oct4_auto_active.syms, [1,9,14]), legend = true)
#
# # Plot O,D,TET,N,NT
# f1 = plot(osol, vars = getindex(Oct4_auto_active.syms, [1,3,9,14,19]))
# f2 = plot(osol, vars = getindex(Oct4_auto_active.syms, [1,14]), xlabel = "time");
# l = @layout [ a; b]
# plot(f1,f2,layout =l,legend = true)
