# Rate constant and DF sampling ranges, from within physiological ranges of biological systems
const KF_RANGE = (1e-3, 1e2) #forward rate constants, units of 1/µM/s
const KR_RANGE = (1e-3, 1e3) #reverse rate constants, units of 1/s
const KCAT_RANGE = (1e-3, 1e3) #catalytic rate constants, units of 1/s
const DF_RANGE = (1.0, 1e4) #dimensionality factor, defined as V/Ah, where V is the volume of the reaction system and A is the surface area of the membrane, and h is a thermodynamic length scale (typically 10nm) specific to each binding pair (we assume h = 10nm for all binding pairs).

# Species sampling ranges, in units of µM
const L_RANGE = (1e-1, 1e2) # PIP/PIP2
const K_RANGE = (1e-1, 1e2) # PIP5K
const P_RANGE = (1e-1, 1e2) # Synaptojanin
const A_RANGE = (1e-1, 1e2) # AP2

# Bounds dictionary for input parameters and initial species concentrations
const BOUNDS = dictionary([:kfᴸᴬ => KF_RANGE, :krᴸᴬ => KR_RANGE, :kfᴸᴷ => KF_RANGE, :krᴸᴷ => KR_RANGE, :kcatᴸᴷ => KCAT_RANGE, :kfᴸᴾ => KF_RANGE, :krᴸᴾ => KR_RANGE, :kcatᴸᴾ => KCAT_RANGE, :kfᴬᴷ => KF_RANGE, :krᴬᴷ => KR_RANGE, :kfᴬᴾ => KF_RANGE, :krᴬᴾ => KR_RANGE, :DF => DF_RANGE, :L => L_RANGE, :K => K_RANGE, :P => P_RANGE, :A => A_RANGE])

# Default values for input parameters and initial species concentrations
const DEFAULTS = dictionary([:kfᴸᴬ => 1.8438562425888723, :krᴸᴬ => 14.208807856126104, :kfᴸᴷ => 0.0015636561694145224, :krᴸᴷ => 79.75320912981275, :kcatᴸᴷ => 49.35726678915171, :kfᴸᴾ => 24.785882197512002, :krᴸᴾ => 122.13694437630977, :kcatᴸᴾ => 39.8567432994366, :kfᴬᴷ => 8.915650787390948, :krᴬᴷ => 0.03319110625610763, :kfᴬᴾ => 0.024471394984824354, :krᴬᴾ => 0.03517891784666366, :DF => 2582.7847577988523, :L => 43.15801010657622, :K => 21.026988187748856, :P => 5.484455252993349, :A => 9.83564656410484])


#* Base reaction network, AKA the "restricted" 12-variable model that doesn't allow A and enzyme binding off the membrane. 
# This is the model used for the mechanistic analysis, as the additional peripheral reactions don't seem to contribute any new oscillatory behaviors or features, just makes this model more sensitive. 
const base_rn = @reaction_network base_rn begin

    @parameters begin
        kfᴸᴬ::Float64 = $(DEFAULTS[:kfᴸᴬ]), [bounds=$(BOUNDS[:kfᴸᴬ])] 
        krᴸᴬ::Float64 = $(DEFAULTS[:krᴸᴬ]), [bounds=$(BOUNDS[:krᴸᴬ])] 
        kfᴸᴷ::Float64 = $(DEFAULTS[:kfᴸᴷ]), [bounds=$(BOUNDS[:kfᴸᴷ])] 
        krᴸᴷ::Float64 = $(DEFAULTS[:krᴸᴷ]), [bounds=$(BOUNDS[:krᴸᴷ])] 
        kcatᴸᴷ::Float64 = $(DEFAULTS[:kcatᴸᴷ]), [bounds=$(BOUNDS[:kcatᴸᴷ])]
        kfᴸᴾ::Float64 = $(DEFAULTS[:kfᴸᴾ]), [bounds=$(BOUNDS[:kfᴸᴾ])]
        krᴸᴾ::Float64 = $(DEFAULTS[:krᴸᴾ]), [bounds=$(BOUNDS[:krᴸᴾ])] 
        kcatᴸᴾ::Float64 = $(DEFAULTS[:kcatᴸᴾ]), [bounds=$(BOUNDS[:kcatᴸᴾ])] 
        kfᴬᴷ::Float64 = $(DEFAULTS[:kfᴬᴷ]), [bounds=$(BOUNDS[:kfᴬᴷ])] 
        krᴬᴷ::Float64 = $(DEFAULTS[:krᴬᴷ]), [bounds=$(BOUNDS[:krᴬᴷ])] 
        kfᴬᴾ::Float64 = $(DEFAULTS[:kfᴬᴾ]), [bounds=$(BOUNDS[:kfᴬᴾ])] 
        krᴬᴾ::Float64 = $(DEFAULTS[:krᴬᴾ]), [bounds=$(BOUNDS[:krᴬᴾ])]
        DF::Float64 = $(DEFAULTS[:DF]), [description="Unitless dimensionality factor", bounds=$(BOUNDS[:DF])]
    end 

    @species L(t) = $(DEFAULTS[:L]) [description="PIP", bounds=$(BOUNDS[:L])] K(t) = $(DEFAULTS[:K]) [description="PIP5K", bounds=$(BOUNDS[:K])] P(t) = $(DEFAULTS[:P]) [description="Synaptojanin", bounds=$(BOUNDS[:P])] A(t) = $(DEFAULTS[:A]) [description="AP2", bounds=$(BOUNDS[:A])] Lp(t) = 0.0 [description="PIP2", tunable=false] LpA(t) = 0.0 [description="PIP2-AP2", tunable=false] LK(t) = 0.0 [description="PIP-Kinase", tunable=false] LpP(t) = 0.0 [description="PIP2-Phosphatase", tunable=false] LpAK(t) = 0.0 [description="PIP2-AP2-Kinase", tunable=false] LpAP(t) = 0.0 [description="PIP2-AP2-Phosphatase", tunable=false] LpAKL(t) = 0.0 [description="PIP2-AP2-Kinase-PIP", tunable=false] LpAPLp(t) = 0.0 [description="PIP2-AP2-Phosphatase-PIP2", tunable=false] 

    @observables Amem ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp)
    

    #* ALIASES: L = PIP, Lp = PIP2, K = Kinase, P = Phosphatase, A = AP2 
    #* reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
    #* DF is the dimensionless dimensionality factor that scales the forward rate constants when in 2D, and is proportional to the membrane surface to volume ratio of the reaction system. 

    #$ Reaction 1 (Kdᴸᴬ)
    (kfᴸᴬ,krᴸᴬ), Lp + A <--> LpA # Lp binding to AP2 adaptor, providing anchor for phosphatase to phosphorylate kinase to localize to the membrane.

    #< Enzyme localization to membrane. Only instance of competition between kinases and phosphatases for LpA anchor points.
    #$ Reaction 2 (Kdᴬᴷ)
    (kfᴬᴷ,krᴬᴷ), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    #$ Reaction 3 (Kdᴬᴾ)
    (kfᴬᴾ,krᴬᴾ), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 

    #< Enzymatic reactions controlling the L to Lp ratio of the lipid membrane, but operating from 3D cytosol.
    #$ Reaction 4a (Kmᴸᴷ)
    (kfᴸᴷ,krᴸᴷ), L + K <--> LK # L binding to kinase in 3D cytosol.
    kcatᴸᴷ, LK --> Lp + K # L phosphorylation by kinase into Lp

    #$ Reaction 5a (Kmᴸᴾ)
    (kfᴸᴾ,krᴸᴾ), Lp + P <--> LpP # Lp binding to phosphatase in 3D cytosol.
    kcatᴸᴾ, LpP --> L + P # L dephosphorylation by phosphatase

    #< Membrane-localized versions of enzymatic reactions, scaled by DF!
    #$ Reaction 4b (Kmᴸᴷ / DF)
    (kfᴸᴷ*DF,krᴸᴷ), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by DF. DF affects enzyme-substrate binding only, and so multiplied by the forward rate constant (thus is on the denominator of the effective Kmᴸᴷ). POSITIVE FEEDBACK
    kcatᴸᴷ, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 

    #$ Reaction 5b (Kmᴸᴾ / DF)
    (kfᴸᴾ*DF,krᴸᴾ), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by DF. AMPLIFIED NEGATIVE FEEDBACK
    kcatᴸᴾ, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality
end  

#* Peripheral reactions that allow A and enzymes to bind in the cytosol, extending the base reaction network with this generates the full 16-variable "unrestricted" model
const peripheral_rn = @reaction_network peripheral_rn begin
    @parameters begin
        kfᴸᴬ::Float64 = $(DEFAULTS[:kfᴸᴬ]), [bounds=$(BOUNDS[:kfᴸᴬ])] 
        krᴸᴬ::Float64 = $(DEFAULTS[:krᴸᴬ]), [bounds=$(BOUNDS[:krᴸᴬ])] 
        kfᴸᴷ::Float64 = $(DEFAULTS[:kfᴸᴷ]), [bounds=$(BOUNDS[:kfᴸᴷ])] 
        krᴸᴷ::Float64 = $(DEFAULTS[:krᴸᴷ]), [bounds=$(BOUNDS[:krᴸᴷ])] 
        kcatᴸᴷ::Float64 = $(DEFAULTS[:kcatᴸᴷ]), [bounds=$(BOUNDS[:kcatᴸᴷ])]
        kfᴸᴾ::Float64 = $(DEFAULTS[:kfᴸᴾ]), [bounds=$(BOUNDS[:kfᴸᴾ])]
        krᴸᴾ::Float64 = $(DEFAULTS[:krᴸᴾ]), [bounds=$(BOUNDS[:krᴸᴾ])] 
        kcatᴸᴾ::Float64 = $(DEFAULTS[:kcatᴸᴾ]), [bounds=$(BOUNDS[:kcatᴸᴾ])] 
        kfᴬᴷ::Float64 = $(DEFAULTS[:kfᴬᴷ]), [bounds=$(BOUNDS[:kfᴬᴷ])] 
        krᴬᴷ::Float64 = $(DEFAULTS[:krᴬᴷ]), [bounds=$(BOUNDS[:krᴬᴷ])] 
        kfᴬᴾ::Float64 = $(DEFAULTS[:kfᴬᴾ]), [bounds=$(BOUNDS[:kfᴬᴾ])] 
        krᴬᴾ::Float64 = $(DEFAULTS[:krᴬᴾ]), [bounds=$(BOUNDS[:krᴬᴾ])]
        DF::Float64 = $(DEFAULTS[:DF]), [description="Unitless dimensionality factor", bounds=$(BOUNDS[:DF])]
    end 

    @species L(t) = $(DEFAULTS[:L]) [description="PIP", bounds=$(BOUNDS[:L])] K(t) = $(DEFAULTS[:K]) [description="PIP5K", bounds=$(BOUNDS[:K])] P(t) = $(DEFAULTS[:P]) [description="Synaptojanin", bounds=$(BOUNDS[:P])] A(t) = $(DEFAULTS[:A]) [description="AP2", bounds=$(BOUNDS[:A])] Lp(t) = 0.0 [description="PIP2", tunable=false] LpA(t) = 0.0 [description="PIP2-AP2", tunable=false] LK(t) = 0.0 [description="PIP-Kinase", tunable=false] LpP(t) = 0.0 [description="PIP2-Phosphatase", tunable=false] LpAK(t) = 0.0 [description="PIP2-AP2-Kinase", tunable=false] LpAP(t) = 0.0 [description="PIP2-AP2-Phosphatase", tunable=false] LpAKL(t) = 0.0 [description="PIP2-AP2-Kinase-PIP", tunable=false] LpAPLp(t) = 0.0 [description="PIP2-AP2-Phosphatase-PIP2", tunable=false] AK(t) = 0.0 [description = "AP2-Kinase", tunable = false] AP(t) = 0.0 [description = "AP2-Phosphatase", tunable = false] AKL(t) = 0.0 [description = "AP2-Kinase-PIP", tunable = false] APLp(t) = 0.0 [description = "AP2-Phosphatase-PIP2", tunable = false] 

    #- Amem is the fraction of AP2 that is membrane-bound, AKA complexed in a 2D membrane-bound complex. The denominator simply represents the total amount of AP2 in the mass conserved system.
    @observables begin
        #* Old definition of Amem that doesn't include AKL and APLp, but is much better behaved during optimization (smoother fitness landscape and less instability). This is what is scored via the FFT based fitness function
        Amem_old ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp)
        #* New definition includes AKL and APLp, which should be included in order to match how any real experiemnt would measure it, even though A is only indirectly localized to the membrane in those species. Period and amplitude are measured from this definition, and all downstream analysis assumes this definition.
        Amem ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp + AKL + APLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp)
    end

    (kfᴸᴬ,krᴸᴬ), Lp + AK <--> LpAK
    (kfᴸᴬ*DF,krᴸᴬ), Lp + AKL <--> LpAKL

    (kfᴸᴬ,krᴸᴬ), Lp + AP <--> LpAP
    (kfᴸᴬ*DF,krᴸᴬ), Lp + APLp <--> LpAPLp

    (kfᴬᴷ,krᴬᴷ), A + K <--> AK
    (kfᴬᴾ,krᴬᴾ), A + P <--> AP

    (kfᴬᴷ,krᴬᴷ), A + LK <--> AKL
    (kfᴬᴾ,krᴬᴾ), A + LpP <--> APLp

    (kfᴬᴷ*DF,krᴬᴷ), LpA + LK <--> LpAKL
    (kfᴬᴾ*DF,krᴬᴾ), LpA + LpP <--> LpAPLp

    (kfᴸᴷ,krᴸᴷ), AK + L <--> AKL #binding of kinase to lipid
    kcatᴸᴷ, AKL --> Lp + AK #phosphorylation of lipid

    (kfᴸᴾ,krᴸᴾ), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcatᴸᴾ, APLp --> L + AP #dephosphorylation of lipid
end

"Full model is the peripheral reaction network extended by the base reaction network"
# Create a composed reaction network instead of extending completed systems
const fullrn = @reaction_network fullrn begin
    @parameters begin
        kfᴸᴬ::Float64 = $(DEFAULTS[:kfᴸᴬ]), [bounds=$(BOUNDS[:kfᴸᴬ])] 
        krᴸᴬ::Float64 = $(DEFAULTS[:krᴸᴬ]), [bounds=$(BOUNDS[:krᴸᴬ])] 
        kfᴸᴷ::Float64 = $(DEFAULTS[:kfᴸᴷ]), [bounds=$(BOUNDS[:kfᴸᴷ])] 
        krᴸᴷ::Float64 = $(DEFAULTS[:krᴸᴷ]), [bounds=$(BOUNDS[:krᴸᴷ])] 
        kcatᴸᴷ::Float64 = $(DEFAULTS[:kcatᴸᴷ]), [bounds=$(BOUNDS[:kcatᴸᴷ])]
        kfᴸᴾ::Float64 = $(DEFAULTS[:kfᴸᴾ]), [bounds=$(BOUNDS[:kfᴸᴾ])]
        krᴸᴾ::Float64 = $(DEFAULTS[:krᴸᴾ]), [bounds=$(BOUNDS[:krᴸᴾ])] 
        kcatᴸᴾ::Float64 = $(DEFAULTS[:kcatᴸᴾ]), [bounds=$(BOUNDS[:kcatᴸᴾ])] 
        kfᴬᴷ::Float64 = $(DEFAULTS[:kfᴬᴷ]), [bounds=$(BOUNDS[:kfᴬᴷ])] 
        krᴬᴷ::Float64 = $(DEFAULTS[:krᴬᴷ]), [bounds=$(BOUNDS[:krᴬᴷ])] 
        kfᴬᴾ::Float64 = $(DEFAULTS[:kfᴬᴾ]), [bounds=$(BOUNDS[:kfᴬᴾ])] 
        krᴬᴾ::Float64 = $(DEFAULTS[:krᴬᴾ]), [bounds=$(BOUNDS[:krᴬᴾ])]
        DF::Float64 = $(DEFAULTS[:DF]), [description="Unitless dimensionality factor", bounds=$(BOUNDS[:DF])]
    end 

    @species L(t) = $(DEFAULTS[:L]) [description="PIP", bounds=$(BOUNDS[:L])] K(t) = $(DEFAULTS[:K]) [description="PIP5K", bounds=$(BOUNDS[:K])] P(t) = $(DEFAULTS[:P]) [description="Synaptojanin", bounds=$(BOUNDS[:P])] A(t) = $(DEFAULTS[:A]) [description="AP2", bounds=$(BOUNDS[:A])] Lp(t) = 0.0 [description="PIP2", tunable=false] LpA(t) = 0.0 [description="PIP2-AP2", tunable=false] LK(t) = 0.0 [description="PIP-Kinase", tunable=false] LpP(t) = 0.0 [description="PIP2-Phosphatase", tunable=false] LpAK(t) = 0.0 [description="PIP2-AP2-Kinase", tunable=false] LpAP(t) = 0.0 [description="PIP2-AP2-Phosphatase", tunable=false] LpAKL(t) = 0.0 [description="PIP2-AP2-Kinase-PIP", tunable=false] LpAPLp(t) = 0.0 [description="PIP2-AP2-Phosphatase-PIP2", tunable=false] AK(t) = 0.0 [description = "AP2-Kinase", tunable = false] AP(t) = 0.0 [description = "AP2-Phosphatase", tunable = false] AKL(t) = 0.0 [description = "AP2-Kinase-PIP", tunable = false] APLp(t) = 0.0 [description = "AP2-Phosphatase-PIP2", tunable = false] 

    @observables begin
        Amem_old ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp)
        Amem ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp + AKL + APLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp)
    end

    # Base reaction network reactions
    (kfᴸᴬ,krᴸᴬ), Lp + A <--> LpA
    (kfᴬᴷ,krᴬᴷ), LpA + K <--> LpAK
    (kfᴬᴾ,krᴬᴾ), LpA + P <--> LpAP
    (kfᴸᴷ,krᴸᴷ), L + K <--> LK
    kcatᴸᴷ, LK --> Lp + K
    (kfᴸᴾ,krᴸᴾ), Lp + P <--> LpP
    kcatᴸᴾ, LpP --> L + P
    (kfᴸᴷ*DF,krᴸᴷ), LpAK + L <--> LpAKL
    kcatᴸᴷ, LpAKL --> Lp + LpAK
    (kfᴸᴾ*DF,krᴸᴾ), Lp + LpAP <--> LpAPLp
    kcatᴸᴾ, LpAPLp --> L + LpAP

    # Peripheral reaction network reactions
    (kfᴸᴬ,krᴸᴬ), Lp + AK <--> LpAK
    (kfᴸᴬ*DF,krᴸᴬ), Lp + AKL <--> LpAKL
    (kfᴸᴬ,krᴸᴬ), Lp + AP <--> LpAP
    (kfᴸᴬ*DF,krᴸᴬ), Lp + APLp <--> LpAPLp
    (kfᴬᴷ,krᴬᴷ), A + K <--> AK
    (kfᴬᴾ,krᴬᴾ), A + P <--> AP
    (kfᴬᴷ,krᴬᴷ), A + LK <--> AKL
    (kfᴬᴾ,krᴬᴾ), A + LpP <--> APLp
    (kfᴬᴷ*DF,krᴬᴷ), LpA + LK <--> LpAKL
    (kfᴬᴾ*DF,krᴬᴾ), LpA + LpP <--> LpAPLp
    (kfᴸᴷ,krᴸᴷ), AK + L <--> AKL
    kcatᴸᴷ, AKL --> Lp + AK
    (kfᴸᴾ,krᴸᴾ), AP + Lp <--> APLp
    kcatᴸᴾ, APLp --> L + AP
end

"""
    make_odeprob(;rn = fullrn, tspan = (0.0, 2000.0))

Create an ODEProblem from the full reaction network. Convenience function for solving the default ODEs.
"""
function make_odeprob(;rn = fullrn, tspan = (0.0, 2000.0))
    odeprob = ODEProblem(rn, Float64[], tspan; jac = true, saveat = 0.1, abstol = 1e-8, reltol = 1e-6)
    return odeprob
end













