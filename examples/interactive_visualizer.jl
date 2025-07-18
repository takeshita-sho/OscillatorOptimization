using DrWatson
@quickactivate "GeometricallyTunableOscillator"
include("/home/jfisch27/Desktop/ThesisStuff/GeometricallyTunableOscillator/src/OscTools/OscTools.jl")
using .OscTools

begin 
    using Catalyst 
    using OrdinaryDiffEq
    using ModelingToolkit
    using DataFrames
    using CSV

    using GLMakie
    using SymbolicIndexingInterface
    using StatsBase
    using FFTW
    using Peaks
    using ADTypes
end

sweepDF_dir = "sweepDF_12-16-24_12853a9657a00ad0c61672af4b66700fecb4acfa"
sweepDF_path = joinpath("cluster", "unrestricted", sweepDF_dir)
sweepDF_datapath = datadir(sweepDF_path)

# filtered_oscillatory_df = CSV.read(joinpath(sweepDF_datapath, "filtered_combined_sweepDF.csv"), DataFrame) 
filtered_oscillatory_df = CSV.read(joinpath(sweepDF_datapath, "filtered_combined_sweepDF-slow_lipid_binding.csv"), DataFrame) 
sort!(filtered_oscillatory_df, [:DF, :Amplitude]; rev = [false, true])

filtered_oscillatory_df = CSV.read(joinpath(sweepDF_datapath, "nerdss_selected_sweepDF.csv"), DataFrame)

# transform!(filtered_oscillatory_df, [:Kmᴸᴷ, :Kmᴸᴾ] => ((x, y) -> max.(x./y, y./x)) => :asymmetry)

# sort!(filtered_oscillatory_df, [:DF, :asymmetry])


function scientific(x::Float64)
    e = floor(Int, log10(abs(x)))
    c = x / 10.0^e
    return "$(round(c, digits=2))e$e"
end

function plot_interactive(df, rn; tend = 1000.0, dt = 0.1)
    odeprob = ODEProblem(rn, Float64[], (0.0, tend); jac = true, reltol = 1e-6, abstol = 1e-8, saveat = dt)
    osys = odeprob.f.sys
    chunk_size = length(odeprob.u0)
    tspan = 0.0:dt:tend

    autodiff = AutoForwardDiff(chunksize=chunk_size)


    parameter_syms = get_parameter_symbols(rn)
    species_syms = get_species_symbols(rn)
    bounds_dict = get_tunable_bounds_dictionary([parameter_syms; species_syms], rn)


    #- Get vector of indices for all the rate constants in the parameter vector needed to compute Kd and Km values 
    parameter_strings = string.(parameter_syms)
    kfᴸᴬ_idx, krᴸᴬ_idx = findall(x -> occursin("ᴸᴬ", x), parameter_strings)
    kfᴬᴷ_idx, krᴬᴷ_idx = findall(x -> occursin("ᴬᴷ", x), parameter_strings)
    kfᴬᴾ_idx, krᴬᴾ_idx = findall(x -> occursin("ᴬᴾ", x), parameter_strings)
    kfᴸᴷ_idx, krᴸᴷ_idx, kcatᴸᴷ_idx = findall(x -> occursin("ᴸᴷ", x), parameter_strings)
    kfᴸᴾ_idx, krᴸᴾ_idx, kcatᴸᴾ_idx = findall(x -> occursin("ᴸᴾ", x), parameter_strings)
    DF_idx = findfirst(x -> occursin("DF", x), parameter_strings)

    #* Dictionary mapping indices to corresponding color
    color_dict = Dict(kfᴸᴬ_idx => :red, krᴸᴬ_idx => :red,
                      kfᴬᴷ_idx => :green, krᴬᴷ_idx => :green,
                      kfᴬᴾ_idx => :purple, krᴬᴾ_idx => :purple,
                      kfᴸᴷ_idx => :orange, krᴸᴷ_idx => :orange, kcatᴸᴷ_idx => :orange,
                      kfᴸᴾ_idx => :blue, krᴸᴾ_idx => :blue, kcatᴸᴾ_idx => :blue,
                      DF_idx => :gray)


    println("Plotting...")
    fig = Figure(;size = (1800, 1400));


    #< GRIDS ##
    plot_grid = fig[1, 1:2] = GridLayout()
    # right_plot_grid = plot_grid[1, 2] = GridLayout()
    slider_grid = fig[2, :] = GridLayout()
    parameter_slider_grid = slider_grid[1, 1] = GridLayout()
    species_slider_grid = slider_grid[1, 2] = GridLayout()
    data_slider_grid = fig[3, :] = GridLayout()
    solver_grid = species_slider_grid[end+1, 1:3] = GridLayout()


    #< AXES ##
    sol_ax = Axis(plot_grid[1, 1], title = "Time Series",
                xlabelsize = 20, xlabel = "Time (s)",
                ylabel = "Normalized Concentration %", ylabelsize = 20,
                limits = (odeprob.tspan, (0.0, 1.0)))

    #< Data Row slider
    #* Make row slider to cycle through the dataframe
    n_rows = nrow(df)
    row_slider = Slider(data_slider_grid[1, 1:2], range = 1:1:n_rows, startvalue = 1, horizontal = true, color_active = :red, color_active_dimmed = :pink, valign = :center, tellheight = true)

    row_slider_label = lift(row_slider.value) do rownumber
        return "Row $(rownumber)/$(n_rows)"
    end 

    #* Label for row slider
    Label(data_slider_grid[2, 1], row_slider_label, fontsize = 20, color = :black, valign = :top, tellheight = true, tellwidth = false)


    buttongrid = data_slider_grid[2, 2] = GridLayout()
    #- Reset button that resets the row to what the slider is at
    reset_button = Button(buttongrid[1, 1], label = "Reset", labelcolor = :red)

    on(reset_button.clicks) do n 
        notify(row_slider.value)
    end


    #* Button to advance the row slider by 1
    advance_button = Button(buttongrid[1, -1], label = "Advance", labelcolor = :green)

    on(advance_button.clicks) do n 
        row_slider.value[] += 1
    end

    #* Button to go back a row
    back_button = Button(buttongrid[1, -2], label = "Back", labelcolor = :blue)

    on(back_button.clicks) do n 
        row_slider.value[] -= 1
    end

    # Mapping of each DF value to the start and end indices of the corresponding group of rows in the dataframe
    DF_group_indices = Dict{Float64, UnitRange{Int}}()
    DF_values = unique(df.DF)
    for (i, DF_value) in enumerate(DF_values)
        idx_range = searchsorted(df.DF, DF_value)
        DF_group_indices[DF_value] = idx_range
    end

    #* Button to advance to the next DF value 
    next_df_button = Button(buttongrid[2, -1], label = "Next DF", labelcolor = :blue)

    on(next_df_button.clicks) do n 
        current_DF = df.DF[row_slider.value[]]
        if current_DF == DF_values[end]
            return
        else
            next_DF = DF_values[findfirst(x -> x == current_DF, DF_values) + 1]
            row_slider.value[] = DF_group_indices[next_DF][1]
        end
    end

    #* Button to go back to the previous DF value 
    prev_df_button = Button(buttongrid[2, -2], label = "Prev DF", labelcolor = :blue)

    on(prev_df_button.clicks) do n 
        current_DF = df.DF[row_slider.value[]]
        if current_DF == DF_values[1]
            return
        else
            prev_DF = DF_values[findfirst(x -> x == current_DF, DF_values) - 1]
            row_slider.value[] = DF_group_indices[prev_DF][1]
        end
    end

    #- Parameter value slider grid
    #* Make slider grid with slider for each parameter in ind 
    parameter_sliders = SliderGrid(parameter_slider_grid[1, 1],
                                ((label = rich(string(label), fontsize = 20), range = logrange(bounds_dict[label]..., 500000), startvalue = odeprob.ps[label]) for label in parameter_syms)...; valign = :top, tellheight = true, tellwidth = true)

    #* Color code the sliders based on corresponding reaction 
    for (idx, slider) in enumerate(parameter_sliders.sliders)
        slider.color_active_dimmed = (color_dict[idx], 0.5)
        slider.color_active = color_dict[idx]
    end


    #- Initial conditions slider grid
    species_sliders = SliderGrid(species_slider_grid[1, 1:3],
                                ((label = rich(string(label), fontsize = 20), range = logrange(bounds_dict[label]..., 500000), startvalue = odeprob[label]) for label in species_syms)...; valign = :top, tellheight = false, tellwidth = true)


    Label(parameter_slider_grid[1, :, Top()], "Parameters", fontsize = 25, tellwidth = false)
    Label(species_slider_grid[1, :, Top()], "Initial Conditions", fontsize = 25, tellwidth = false, tellheight = true, valign = :top)


    observed_p = Observable{typeof(copy(odeprob.ps[parameter_syms]))}(copy(odeprob.ps[parameter_syms]))
    observed_u0 = Observable{typeof(copy(odeprob[species_syms]))}(copy(odeprob[species_syms]))


    # Just collecting slider observables into vectors
    parameter_slider_observables = [s.value for s in parameter_sliders.sliders]
    species_slider_observables = [s.value for s in species_sliders.sliders]


    # Needed to prevent circular loop between updates to the sliders and the row slider
    is_programmatic_update = Ref{Bool}(false)

    # On row slider update, update the parameter and species sliders
    on(row_slider.value; update = true, priority = 1) do rownumber
        # is_programmatic_update[] = true
        dfrow = df[rownumber, :]
        p = collect(dfrow[parameter_syms])
        u0 = collect(dfrow[species_syms])

        observed_p[] .= p
        observed_u0[] .= u0

        notify(observed_p)
    end

    Makie.on_latest(row_slider.value; update = true, spawn = false) do _
        is_programmatic_update[] = true
        set_close_to!.(parameter_sliders.sliders, to_value(observed_p))
        set_close_to!.(species_sliders.sliders, to_value(observed_u0))
        is_programmatic_update[] = false
    end

    onany(parameter_slider_observables...) do slvalues...
        if !is_programmatic_update[]
            # println("Parameter slider updated")
            observed_p[] .= slvalues
            notify(observed_p)
        end
    end
    
    onany(species_slider_observables...) do slvalues...
        if !is_programmatic_update[]
            # println("Species slider updated")
            observed_u0[] .= slvalues
            notify(observed_u0)
        end
    end

    # Create a slider for reltol (using log scale since it's typically varied by orders of magnitude)
    reltol_slider = Slider(solver_grid[1, 1:3],  # Make it span all columns like the species sliders
                        range = 10.0 .^ range(-12, -1, length=1000),
                        startvalue = 1e-6,
                        horizontal = true)

    # Add a label to show current reltol value
    reltol_label = lift(reltol_slider.value) do val
        return "RelTol: $(scientific(val))"
    end
    Label(solver_grid[2, 1:3], reltol_label, fontsize = 20)  # Make label span columns too

    #< Solve and get variable solutions for plotting 
    set_parameters! = setp(odeprob, parameter_syms)
    set_species! = setu(odeprob, species_syms)
    sol = lift(observed_p, observed_u0, reltol_slider.value) do observed_p, observed_u0, reltol
        set_parameters!(odeprob, observed_p)
        set_species!(odeprob, observed_u0)
        # Remake problem with new reltol
        new_prob = remake(odeprob; reltol=reltol)
        return solve(new_prob, Rodas5P(autodiff = autodiff))
    end

    # For peak finding and period/amplitude calculation 
    get_Amem = getu(osys, :Amem)


    #< Timeseries
    Amem = lift(sol) do sol
        return get_Amem(sol)
    end

    amem_peaks_and_troughs = lift(Amem) do amem
        return find_amem_peaks_no_simd(amem)
    end

    amem_period = lift(amem_peaks_and_troughs) do amem_peaks_and_troughs
        amem_peaks_indices = amem_peaks_and_troughs[1]
        if length(amem_peaks_indices) < 2
            return 0.0
        else
            return compute_period(tspan, amem_peaks_indices)
        end
    end

    amem_amplitude = lift(amem_peaks_and_troughs) do amem_peaks_and_troughs
        amem_peaks_indices, amem_troughs_indices = amem_peaks_and_troughs
        amem_peaks_heights = Amem.val[amem_peaks_indices]
        amem_troughs_heights = Amem.val[amem_troughs_indices]
        if length(amem_peaks_indices) < 2
            return 0.0
        else
            return OscTools.compute_amplitude(amem_peaks_heights, amem_troughs_heights)
        end
    end


    #< Period and amplitude labeling 
    Amem_period_label = lift(amem_period) do period
        period_string = period == 0.0 ? "No Periodic" : "Period: " * string(round(Int, period)) * " (s)"
        return period_string
    end


    Amem_amplitude_label = lift(amem_amplitude) do amplitude
        amplitude_string = amplitude == 0.0 ? "No Amplitude" : "Amplitude: " * string(round(amplitude, digits = 3))
        return amplitude_string
    end

    labels_grid = species_slider_grid[end+1, 1] = GridLayout()

    Label(labels_grid[end+1, 2], Amem_period_label, fontsize = 25, tellwidth = false, justification = :left)
    Label(labels_grid[end+1, 2], Amem_amplitude_label, fontsize = 25, tellwidth = false, justification = :left)


    #< Time series plot 
    lines!(sol_ax, tspan, Amem, color = :black, linewidth = 3, label = "Amem")
    axislegend(sol_ax)
    amem_peak_data = lift(amem_peaks_and_troughs) do amem_peaks_and_troughs
        amem_peaks, _ = amem_peaks_and_troughs
        amem_peaks_heights = Amem.val[amem_peaks]
        if isempty(amem_peaks)
            return Point2f[]
        else
            return [Point2f(tspan[i], h) for (i, h) in zip(amem_peaks, amem_peaks_heights)]
        end
    end

    amem_trough_data = lift(amem_peaks_and_troughs) do amem_peaks_and_troughs
        _, amem_troughs = amem_peaks_and_troughs
        amem_troughs_heights = Amem.val[amem_troughs]
        if isempty(amem_troughs)
            return Point2f[]
        else
            return [Point2f(tspan[i], h) for (i, h) in zip(amem_troughs, amem_troughs_heights)]
        end
    end

    scatter!(sol_ax, amem_peak_data, color = :red, markersize = 16)
    scatter!(sol_ax, amem_trough_data, color = :blue, markersize = 16)

    DataInspector()
    fig
end

plot_interactive(filtered_oscillatory_df, fullrn; tend = 2000.0, dt = 0.1)



