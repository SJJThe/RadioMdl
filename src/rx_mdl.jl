
"""
    freq_range(freq_res::T,
               freq_center::T,
               nb_freq_bins::Int,
               bandwidth::T) where T

Yields the frequency range centered on freq_center and of length nb_freq_bins.

"""
function freq_range(freq_res::T,
    freq_center::T,
    nb_freq_bins::Int,
    bandwidth::T) where T

    rng = range(-bandwidth/2 + freq_res/2, bandwidth/2 - freq_res/2, length=nb_freq_bins)
    return freq_center .+ rng
end



"""
    power_to_temperature(power::T,
                         bandwidth::T) where T

Yields the temperature of a source given its power in watt, bandwidth in hertz.
"""
function power_to_temperature(power::T,
    bandwidth::T) where T

    return power / (k_boltz * bandwidth)
end



"""
    temperature_to_power(temp::T,
                         bandwidth::T) where T

Yields the power of a source given its temperature in kelvin, bandwidth in
hertz.

"""
function temperature_to_power(temp::T,
    bandwidth::T) where T

    return k_boltz * bandwidth * temp
end



"""
    adc_noise_temperature(Vfs::T,
                          nb_bits::Int,
                          nb_freq_bins::Int;
                          instru_imp::T = impedance) where T

Yields the noise temperature of the ADC for a given full-scale voltage 'Vfs',
number of bits 'nb_bits' and number of frequency bins 'nb_freq_bins'. The 
instrument impedance 'instru_imp' is set to 50 Ohm by default.

"""
function adc_noise_temperature(Vfs::T,#Full scale ADC voltage
    nb_bits::Int,
    nb_freq_bins::Int; # in FFT slice
    instru_imp::T = impedance) where T
    
    # noise power of the ADC
    P_adc = Vfs^2 / (12*instru_imp) * 2^(-2. *nb_bits)

    return power_to_temperature(P_adc, T(nb_freq_bins))
end



"""
    friis_noise_temp(stages::Tuple{T,T}...) where T<:AbstractFloat

Yields the total noise temperature of a receiver given the noise temperatures
and gains of each stage.

"""
function friis_noise_temp(stages::Tuple{T,T}...) where T<:AbstractFloat
    
    # Check that at least one stage is provided
    if isempty(stages)
        return zero(T)
    end

    T_total = zero(T)
    current_cumulative_gain = one(T)
    
    # Iterate through each stage (T_i, G_i)
    for i in eachindex(stages)
        T_i, G_i = stages[i]
        
        # The first stage adds its full noise temperature.
        # Subsequent stages' noise temperatures are divided by the gain 
        # of ALL preceding stages.
        if i == 1
            T_total += T_i
        else
            T_total += T_i / current_cumulative_gain
        end

        # Update the cumulative gain for the *next* stage calculation
        current_cumulative_gain *= G_i
    end

    return T_total
end

