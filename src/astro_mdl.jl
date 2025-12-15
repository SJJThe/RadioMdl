
"""
"""
freq_to_wave(freq::T) where {T} = speed_c / freq
wave_to_freq(wave::T) where {T} = speed_c / wave



"""
power in watt, bandwidth in hertz
"""
function power_to_temperature(power::T,
    bandwidth::T) where T

    return power / (k_boltz * bandwidth)
end



"""
"""
function temperature_to_power(temp::T,
    bandwidth::T) where T

    return k_boltz * bandwidth * temp
end



"""
    flux_to_temperature(flux::T,
                        effective_apperture::T) where T

estimates the temperature of a point-like source from its flux and the antenna effective
aperture. flux must be in Jansky

"""
function flux_to_temperature(flux::T,
    effective_apperture::T) where T

    return flux*1e-26 / (2*k_boltz) * effective_apperture
end



"""
in Jansky
"""
function temperature_to_flux(temp::T,
    effective_apperture::T) where T

    return (2*k_boltz) * temp / effective_apperture * 1e26
end



"""
"""
function get_geometric_effective_aperture(aperture_efficiency::T,
    diameter::T) where T

    @assert T(0) <= aperture_efficiency <= T(1)
    return aperture_efficiency * pi * (diameter/2)^2
end



"""
    estim_casA_flux(center_freq::T) where T

estimates the flux of Cas A, given a frequency. Based on Baars et al. 2014

"""
function estim_casA_flux(center_freq::T;
    year::Int = 2025) where T
    
    if 22e6 < center_freq < 300e6
        a = 5.625
        var_a = .021^2
        b = -.634
        var_b = .015^2
        c = -.023
        var_c = .001^2
    elseif 300e6 < center_freq < 31e9
        a = 5.880
        var_a = .025^2
        b = -0.792
        var_b = .007^2
        c = 0.
        var_c = 0.
    else
        @error "the model is not valid for this frequency"
    end
    
    # decay
    decay = 0.97 - 0.3*log10(center_freq*1e-9) # in %/year since 1980
    var_decay = .04^2 + .04^2*log10(center_freq*1e-9)^2
    
    # log flux
    log_S_Jy = a + b*log10(center_freq*1e-6) + c*log10(center_freq*1e-6)^2
    var_log_S_Jy = var_a + var_b*log10(center_freq*1e-6)^2 + 
                   var_c*log10(center_freq*1e-6)^4

    # flux in Jy
    S_Jy = 10^log_S_Jy * (1 - decay*(year-1980)/100)
    var_S_Jy = var_log_S_Jy / S_Jy^2 * (1 - decay*(year-1980)/100)^2 + 
               10^log_S_Jy * var_decay 
    
    return S_Jy, var_S_Jy
end



"""
"""
function estim_virgoA_flux(center_freq::T) where T
    return 10^(5.023 - 0.856*log10(center_freq*1e-6))
end



"""
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

