
"""
    estim_temp(flux::T,
               effective_apperture::T) where T

estimates the temperature of a point-like source from its flux and the antenna effective
aperture. flux must be in Jansky

"""
function estim_temp(flux::T,
    effective_apperture::T) where T

    return flux*1e-26 / (2*k_boltz) * effective_apperture
end

function estim_temp(flux::Real,
    obs::Observation)

    instru = get_instrument(obs)
    frequency = get_center_freq(instru)
    ant = get_antenna(instru)
    max_gain = get_boresight_gain(ant)
    A_eff_max = gain_to_effective_aperture(max_gain, frequency)

    return estim_temp(flux, A_eff_max)
end



"""
power in watt, bandwidth in hertz
"""
function power_to_temperature(power::T,
    bandwidth::T) where T

    return power / (k_boltz*bandwidth)
end



"""
"""
function temperature_to_power(temp::T,
    bandwidth::T) where T

    return k_boltz*bandwidth * temp
end



"""
in Jansky
"""
function temperature_to_flux(temp::T,
    effective_apperture::T) where T

    return (2*k_boltz) * temp / effective_apperture * 1e26
end



"""
    estim_casA_flux(center_freq::T) where T

estimates the flux of Cas A, given a frequency. Based on Baars et al. 1977

"""
function estim_casA_flux(center_freq::T) where T

    decay = 0.97 - 0.3*log10(center_freq*1e-9) # in %/year since 1980
    
    return 10^(5.745 - 0.770*log10(center_freq*1e-6))*(1 - decay*43/100) # in Jy
end



"""
"""
function estim_virgoA_flux(center_freq::T) where T
    return 10^(5.023 - 0.856*log10(center_freq*1e-6))
end


