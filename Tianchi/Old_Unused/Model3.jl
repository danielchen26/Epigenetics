using DifferentialEquations, DiffEqBiological
using ParameterizedFunctions
# model 3

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
  aₐ,               D → D̄           # De-novo Methylation
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
 # k₃,             NT + D̄ → D + NT    # Nanog-TET complex de-Methylate promoter of the Oct4
 k₃,             NT + D̄ → Dₕ + NT    # Nanog-TET complex de-Methylate promoter of the Oct4
 # Nanog + TET1 induce not only Oct4 but TET1
 k₄,             NT + Dₜ → Dᵗ + NT
end a₀ d₀ r₁ r₀ α₀ α₁ k β δₒ aₐ a₂ d₂ k₂ ϕ₀ ϕ₁ r₂ rᵦ δₜ ζ₁ ζ₀ η₁ η₀ k₃ k₄
