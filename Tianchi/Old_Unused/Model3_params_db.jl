# Dataframe for models parameters
using DataFrames
params = DataFrame(names=["a₀", "d₀", "r₁", "r₀", "α₀","α₁", "k",  "β", "δₒ",  "aₐ",  "a₂","d₂","k₂","ϕ₀","ϕ₁", "r₂","rᵦ","δₜ", "ζ₁", "ζ₀","η₁", "η₀" ,"k₃", "k₄"],
                values=[100000, 1000, 100000,  1000,  0,  0.2,  0.5,  0.1,   0.09,  100, 100, 100,  100, 0,   0.5, 1e5,  1e3,  0.0693,  100,  100, 100, 100, 100, 100])


param_dict =
DataFrame(
Dict(
"a₀"  =>     100000.0,  # OCT4 Dimerization on
"d₀"  =>     1000.0,    # OCT4 Dimerization off
"r₁"  =>     100000.0,  # Oct4 Auto-activation on
"r₀"  =>     1000.0,    # Oct4 Auto-activation off
"α₀"  =>     0.0,       # Oct4 mRNA (leaky)
"α₁"  =>     0.2,       # Oct4 mRNA on (activ.)
"k "  =>    0.5,        # Oct4 protein on (activ.)
"β "  =>    0.1,        # Dilution rate
"δₒ"  =>     0.09,       # OCT4 degradation
"aₐ"  =>     100.0,      # De-novo Methylation
# "d₁"  =>     100.0,
# "k₁"  =>     100.0,
"a₂"  =>     100.0,     # Active de-Methylation on
"d₂"  =>     100.0,     # Active de-Methylation off
"k₂"  =>     100.0,     # Active de-Methylation oxy
"ϕ₀"  =>     0.0,       # TET Transcription/Translation (leaky)
"ϕ₁"  =>     0.5,       #TET Transcription/Translation (activ.)
"r₂"  =>     100000.0,  # OCT4 Activating TET on
"rᵦ"  =>     1000.0,    # TET Transcription/Translation (activ.)
"δₜ"  =>     0.0693,    # TET Degradation
"ζ₁"  =>     100.0,    # Oct4 dimmer activate Nanog on
"ζ₀"  =>     100.0,    # Oct4 dimmer activate Nanog off
"η₁"  =>     100.0,    # NT forming on
"η₀"  =>     100.0,    # NT forming off
"k₃"  =>     100.0,    # NT assist Methylate Oct4 promoter
"k₄"  =>     100.0    # NT assist activate TET promoter
))
