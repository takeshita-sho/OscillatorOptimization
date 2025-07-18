#< FITNESS FUNCTION HELPER FUNCTIONS ##
"""Get summed differences of the peak values from the FFT of the solution"""
function getDif(peakvals)
    differences = diff(peakvals)
    sum(abs.(differences)) #/ length(peakvals)
end

"""Get summed standard deviation of peaks values from the FFT of the solution"""
function getSTD(fft_peakindxs::Vector{Int}, fftData, window::Int =1) 
    arrLen = length(fftData)

    sum_std = sum(std(@view fftData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs; init=0.0) #* sum rolling window of standard deviations

    return sum_std
end 
#> END OF FITNESS FUNCTION HELPER FUNCTIONS ##



#< FFT HELPER FUNCTIONS ##
"""
    getFrequencies(timeseries)
Return the real-valued FFT of a 1D ODE solution, will be half the length of the timeseries
"""
function getFrequencies(timeseries::AbstractVector{Float64}, rfft_plan::FFTW.rFFTWPlan)
    rfft_result = rfft_plan * timeseries
    norm_val = length(rfft_result) #* normalize by length of RFFT
    return abs.(rfft_result) ./ norm_val
end

function getFrequencies(timeseries::AbstractVector{Float64})
    rfft_result = rfft(timeseries)
    norm_val = length(rfft_result) 
    return abs.(rfft_result) ./ norm_val
end
#> END OF FFT HELPER FUNCTIONS ##



#<< PERIOD AND AMPLITUDE FUNCTIONS ##
function getPerAmp_Amem(sol; min_prominence::Float64 = 0.01)
    Amem = sol[:Amem]
    solt = sol.t

    indx_max, indx_min = findextrema(Amem; min_prominence=min_prominence) 
    min(length(indx_min), length(indx_max)) < 2 && return 0.0, 0.0
    vals_max = Amem[indx_max]
    vals_min = Amem[indx_min]

    #* get the period and amplitude of the solution
    return getPerAmp(solt, indx_max, vals_max, vals_min)
end



"""Calculates the period and amplitude of each individual in the population."""
function getPerAmp(solt, indx_max::Vector{Int}, vals_max::AbstractVector{Float64}, vals_min::AbstractVector{Float64}) 
    return compute_period(solt, indx_max), compute_amplitude(vals_max, vals_min)
end

# function compute_period(solt, indx_max::Vector{Int})
#     mean(solt[indx_max[i+1]] - solt[indx_max[i]] for i in 1:(length(indx_max)-1))
# end

# function compute_amplitude(vals_max::AbstractVector{Float64}, vals_min::AbstractVector{Float64})
#     mean(abs, (vals_max[i] - vals_min[i] for i in 1:min(length(vals_max), length(vals_min))))
# end

# Just returns the time difference between the last two peaks
function compute_period(solt, indx_max::Vector{Int})
    last_peak_time = solt[indx_max[end]]
    second_last_peak_time = solt[indx_max[end-1]]
    return last_peak_time - second_last_peak_time
end

# Just returns the height difference between the last peak and trough
function compute_amplitude(vals_max::AbstractVector{Float64}, vals_min::AbstractVector{Float64})
    last_peak_value = vals_max[end]
    last_trough_value = vals_min[end]
    return last_peak_value - last_trough_value
end
#> END OF PERIOD AND AMPLITUDE FUNCTIONS ##






















