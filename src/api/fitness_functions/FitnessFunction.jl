# Main fitness function method used during optimization. 
# Takes the SavedValues object containing the values saved during `solve_odes!`, and calculates the fitness, period, and amplitude. Populates the zeros vector `phenotype`. 
# Any zeros after the first index representing fitness are returned as is if not filled and indicate NA. 
function calculate_fitness!(phenotype, saved_array, optsys)
    # Use original Amem for FFT
    Amem_original = @view saved_array[1, :]

    fftData = getFrequencies(Amem_original, optsys.rfft_plan)

    #* get the indexes of the peaks in the fft
    # fft_peakindexes = findmaxpeaks(fftData; height = 0.01) 
    fft_peakindexes = find_fft_peaks(fftData)

    #* if there is no signal in the frequency domain, return 0.0s
    # length(fft_peakindexes) < 2 && return phenotype
    if length(fft_peakindexes) < 2
        return phenotype
    end

    #* get the average standard deviation of the window around each peak in frequency domain
    standard_deviation = getSTD(fft_peakindexes, fftData) #/ length(fft_peakindexes)

    #* get the average summed differences between all the peaks in frequency domain
    fft_peakvals = @view fftData[fft_peakindexes]
    sum_diff = getDif(fft_peakvals) 

    #* fitness is the sum of the average standard deviation and the average summed differences between all the peaks in frequency domain
    phenotype[1] = standard_deviation + sum_diff

    #$ AMEM 
    Amem_full = @view saved_array[2, :]

    # indx_max, indx_min = findextrema(Amem_full)
    indx_max, indx_min = find_amem_peaks(Amem_full)

    #- get the values of the peaks and troughs in the time domain
    vals_max = @view Amem_full[indx_max]
    vals_min = @view Amem_full[indx_min]


    #* if there are not enough peaks in the time domain or if the solution is steady state, return 0.0s
    # length(indx_max) < 2 || !check_oscillatory(Amem_full) || !check_oscillation_regularity(Amem_full) && return phenotype
    if length(indx_max) < 2 || length(indx_min) < 2 || !check_oscillatory(Amem_full) || !check_oscillation_regularity(Amem_full)
        return phenotype
    end

    period, amplitude = getPerAmp(optsys.t, indx_max, vals_max, vals_min)
    
    phenotype[2:end] .= period, amplitude
    return phenotype
end

"""
    find_fft_peaks(fftData)

Find the indexes of the peaks in the fftData.
"""
function find_fft_peaks(fftData)
    fft_peakindexes = simplemaxima(fftData)
    peakheights!(fft_peakindexes, fftData[fft_peakindexes]; min = 1e-2)
    # fft_peaks = findmaxima(fftData) |> peakheights!(;min = 1e-2)
    return fft_peakindexes
end

"""
    find_amem_peaks(amem)

Find the peaks and troughs in the amem. Returns a tuple arrays for the indices of the peaks and troughs.
"""
function find_amem_peaks(amem)
    amem_peaks = simplemaxima(amem) 
    amem_troughs = simpleminima(amem) 
    return amem_peaks, amem_troughs
end

"""
    find_amem_peaks_no_simd(amem)

Find the peaks and troughs in the amem. Returns a tuple arrays for the indices of the peaks and troughs. Used when plateaus are present in the amem, throwing errors with SIMD in `simplemaxima` and `simpleminima`.
"""
function find_amem_peaks_no_simd(amem)
    amem_peaks = findmaxima(amem).indices 
    amem_troughs = findminima(amem).indices 
    return amem_peaks, amem_troughs
end

"""
    check_oscillatory(amem)

Check if the solution is oscillatory by checking the standard deviation of the last 200 seconds (2000 points). If the standard deviation is greater than 1e-6, return true, otherwise return false. 
"""
function check_oscillatory(amem)
    last_200_points = @view amem[end-2000:end]
    return std(last_200_points) >= 1e-6 |> float
end


"""
    check_oscillation_regularity(amem::Vector{Float64}; variation_threshold=0.01)

Checks if the given amem vector represents a regular oscillation by analyzing the consistency
of peak heights and peak-to-trough distances. Returns true if it finds any consecutive peaks
that maintain consistent heights and timing within the variation threshold.

Args:
    amem: Vector of amplitude values over time
    variation_threshold: Maximum allowed fractional variation (default: 0.01 = 1%)
        between consecutive peak heights and peak-to-trough distances

Returns:
    Bool: true if regular oscillation is detected, false otherwise
"""
function check_oscillation_regularity(amem; variation_threshold=0.01)
    # Get indices of all peaks and troughs in the signal
    peaks_indices, troughs_indices = find_amem_peaks_no_simd(amem)
    n_peaks = length(peaks_indices)
    
    # Initialize with first peak and its closest trough
    current_trough_idx = 1
    prev_height = amem[peaks_indices[1]]
    prev_peak_trough_distance = abs(troughs_indices[current_trough_idx] - peaks_indices[1])
    
    # Examine each subsequent peak
    for i in 2:n_peaks
        # Find the trough that's closest to the current peak
        # This while loop handles cases where troughs might be irregularly spaced
        while current_trough_idx < length(troughs_indices) && 
              abs(troughs_indices[current_trough_idx + 1] - peaks_indices[i]) < 
              abs(troughs_indices[current_trough_idx] - peaks_indices[i])
            current_trough_idx += 1
        end
        
        # Get current peak's properties
        current_height = amem[peaks_indices[i]]
        current_peak_trough_distance = abs(troughs_indices[current_trough_idx] - peaks_indices[i])
        
        # Calculate relative variations in height and timing
        # Uses mean normalization to get fractional differences
        height_variation = abs(current_height - prev_height) / 
                         mean((current_height, prev_height))
        distance_variation = abs(current_peak_trough_distance - prev_peak_trough_distance) / 
                           mean((current_peak_trough_distance, prev_peak_trough_distance))
        
        # If both variations are below threshold, we've found a regular oscillation
        height_variation < variation_threshold && 
        distance_variation < variation_threshold && return true
        
        # Store current values for next iteration's comparison
        prev_height = current_height
        prev_peak_trough_distance = current_peak_trough_distance
    end
    
    # No regular oscillation found
    return false
end