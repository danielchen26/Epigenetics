using DifferentialEquations, DiffEqBiological
using ParameterizedFunctions
# model 2

O_N_auto = @reaction_network begin
  # Oct4 cycle
  (a₀,d₀),       O + O ↔ O₂          # Oct4 Dimerization
  (r₁,r₀ ),      O₂ + D ↔ Dₒ         # Oct4 Auto-activation
   # α₀*k,         D → D + O          # Oct4 translation (leaky)
   α₁*k,         Dₒ → Dₒ + O          # Oct4 translation (activ.)
  (β,β),        (O₂, O)→∅           # Oct4 Dilution
    δₒ,           O → ∅             # Protein Degradation

 # Oct4 de-Methylation cycle
  # (a₁,d₁),       D + DNMT ↔ C₁      # De-novo Methylation
  # k₁,            C₁ → D̄ + DNMT      # De-novo Methylation
  aₐ,              D → D̄           # De-novo Methylation
  (a₂,d₂),        D̄ +TET ↔ C₂        # Active de-Methylation
  k₂,             C₂ → Dₕ + TET      # Active de-Methylation
  β,              Dₕ →  D            # Passive De-Methylation


# TET protein
  # γD,            ∅ → DNMT          # DNMT Constitutive Transcription
  # ϕ₀,            Dₜ → Dₜ + TET       # TET Transcription/Translation (leaky) ------ set 0
 (r₂,rᵦ),        Dₜ + O₂ ↔ Dᵗ       # OCT4 Activating TET
  ϕ₁,            Dᵗ → Dᵗ + TET      # TET Transcription/Translation (activ.)
  (β, β),         (C₂,TET) → ∅            # TET Dilution
  δₜ,                TET → ∅              # TET Degradation

 # Nanog cycle
 (a₀,d₀),        N + N ↔ N₂          # Nanog Dimerization
 (r₁,r₀ ),       N₂ + Dₙ ↔ Dᴺ        # Nanog Auto-activation
  # α₀*k,         Dₙ → Dₙ + N          # Nanog translation (leaky)
  α₁*k,         Dᴺ → Dᴺ + N         # Nanog translation (activ.)
 (β,β),        (N₂, N)→∅            # Nanog/mRNA Dilution
  δₒ,            N → ∅               # mRNA/Protein Degradation

 # Oct4 Nanog interaction
 (ζ₁, ζ₀),       O₂ + Dₙ ↔ Dᴺ        # Oct4 dimmer activate Nanog promoter

 # Nanog + TET(1) activate OCT4
(η₁, η₀),         N + TET ↔ NT       # Nanog + TET form complex in preparation for Oct4 demethylation
 k₃,             NT + D̄ → D + NT    # Nanog-TET complex de-Methylate promoter of the Oct4
 # Nanog + TET1 induce not only Oct4 but TET1
end a₀ d₀ r₁ r₀ α₀ α₁ k β δₒ aₐ a₂ d₂ k₂ ϕ₀ ϕ₁ r₂ rᵦ δₜ ζ₁ ζ₀ η₁ η₀ k₃




#  ODE version of CRN

O_N_auto_ode = @ode_def begin
    dO = -O*β - O*δₒ - O^2*a₀ + 2*O₂*d₀ + k*Dₒ*α₁;
    dO₂ = Dᴺ*ζ₀ + (1/2)*O^2*a₀ - O₂*d₀ + rᵦ*Dᵗ + r₀*Dₒ - β*O₂ - D*r₁*O₂ - O₂*Dₙ*ζ₁ - r₂*O₂*Dₜ;
    dD = -D*aₐ + r₀*Dₒ + β*Dₕ - D*r₁*O₂ + NT*D̄*k₃;
    dDₒ = -r₀*Dₒ + D*r₁*O₂;
    dD̄ = D*aₐ + d₂*C₂ - NT*D̄*k₃ - TET*D̄*a₂;
    dTET = Dᵗ*ϕ₁ + NT*η₀ - TET*δₜ + d₂*C₂ + k₂*C₂ - β*TET - N*TET*η₁ - TET*D̄*a₂;
    dC₂ = -d₂*C₂ - k₂*C₂ - β*C₂ + TET*D̄*a₂;
    dDₕ = k₂*C₂ - β*Dₕ;
    dDₜ = rᵦ*Dᵗ - r₂*O₂*Dₜ;
    dDᵗ = -rᵦ*Dᵗ + r₂*O₂*Dₜ;
    dN = -N*β - N*δₒ - N^2*a₀ + NT*η₀ + 2*N₂*d₀ - N*TET*η₁ + k*Dᴺ*α₁;
    dN₂ = (1/2)*N^2*a₀ - N₂*d₀ + r₀*Dᴺ - β*N₂ - r₁*N₂*Dₙ;
    dDₙ =  Dᴺ*ζ₀ + r₀*Dᴺ - O₂*Dₙ*ζ₁ - r₁*N₂*Dₙ;
    dDᴺ = -Dᴺ*ζ₀ - r₀*Dᴺ + O₂*Dₙ*ζ₁ + r₁*N₂*Dₙ;
    dNT = -NT*η₀ + N*TET*η₁
end a₀ d₀ r₁ r₀ α₀ α₁ k  β  δₒ aₐ a₂ d₂ k₂ ϕ₀ ϕ₁ r₂ rᵦ δₜ ζ₁ ζ₀ η₁ η₀ k₃














  # :O   => 1
  # :O₂  => 2
  # :D   => 3
  # :Dₒ  => 4
  # :D̄  => 5
  # :TET => 6
  # :C₂  => 7
  # :Dₕ  => 8
  # :Dₜ  => 9
  # :Dᵗ  => 10
  # :N   => 11
  # :N₂  => 12
  # :Dₙ  => 13
  # :Dᴺ  => 14
  # :NT  => 15
