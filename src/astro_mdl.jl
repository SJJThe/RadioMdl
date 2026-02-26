
"""
    freq_to_wave(freq::T) where {T}

Converts frequency to wavelength.

"""
freq_to_wave(freq::T) where {T} = speed_c / freq



"""
    wave_to_freq(wave::T) where {T}

Converts wavelength to frequency.

"""
wave_to_freq(wave::T) where {T} = speed_c / wave



"""
    flux_to_temperature(flux::T,
                        effective_apperture::T) where T

estimates the temperature of a point-like source from its flux and the antenna effective
aperture. flux must be in Jansky

---
    flux_to_temperature(flux::AbstractVector{T},
                        effective_apperture::T) where T

Yields a 'Vector' of temperatures from a vector of fluxes.

"""
function flux_to_temperature(flux::T,
    effective_apperture::T) where T

    return flux*1e-26 / (2*k_boltz) * effective_apperture
end

function flux_to_temperature(flux::AbstractVector{T},
    effective_apperture::T) where T

    return flux_to_temperature.(flux, effective_apperture)
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

estimates the flux of Cas A, given a frequency. Based on Baars et al. 2014

"""
function estim_casA_flux(center_freq::T;
    year::Int = Year(today()).value) where T
    
    if 22e6 < center_freq < 300e6
        a = 5.625
        var_a = .021^2
        b = -.634
        var_b = .015^2
        c = -.023
        var_c = .001^2
    elseif 300e6 < center_freq
        if center_freq > 31e9
            @warn "the model is not valid for frequencies above 31GHz"
        end
        a = 5.880
        var_a = .025^2
        b = -0.792
        var_b = .007^2
        c = 0.
        var_c = 0.
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
    ground_model(alpha_grid::AbstractVector{T},
                 beta_grid::AbstractVector{T},
                 T_ground::AbstractMatrix{T}) where T

Yields a 'SphereMap' structure representing a ground model with temperature
'T_ground'. This can include ground local sources, terrain, etc.

---
    ground_model(alpha_grid::AbstractVector{T},
                 beta_grid::AbstractVector{T},
                 T_ground::T) where T

Yields a 'SphereMap' structure representing a ground model with constant
temperature 'T_ground'. 

---
    ground_model(T_ground::T) where T

Yields a 'SphereMap' structure representing a ground model with constant
temperature 'T_ground' and default alpha and beta grids sampled at 1 degree
resolution.
"""
function ground_model(alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T},
    T_ground::AbstractMatrix{T}) where T

    @assert size(mask) == (length(beta_grid), length(alpha_grid))

    return SphereMap(alpha_grid, beta_grid, T_ground)
end

function ground_model(alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T},
    T_ground::T) where T

    alpha_grid = alpha_grid[alpha_grid .!= T(360.)]
    
    mat_gnd = T_ground .* ones(T, length(alpha_grid), length(beta_grid))

    bellow_horizon = beta_grid .>= 90.
    mat_gnd .*= bellow_horizon'
    
    return SphereMap(alpha_grid, beta_grid, mat_gnd)
end

function ground_model(T_ground::T) where T

    alpha_grid = collect(0:1.:360)
    beta_grid = collect(0:1.:180)

    return ground_model(alpha_grid, beta_grid, T_ground)
end



"""
    atmos_opacity_impact(temp::T,
                         zenith_opacity::T,
                         zenith_angle::T) where T

Yields the 'temp' temperature altered by the atmosphere opacity at
'zenith_angle'.

---
    atmos_opacity_impact(temp::AbstractArray{T},
                         zenith_opacity::T,
                         zenith_angle::T) where T

Yields the 'temp' temperature altered by the atmosphere opacity at
'zenith_angle' for each element of 'temp'. 'temp being a 'AbstractArray', it can
be a 'DimArray' if needed.

---
    atmos_opacity_impact(temp::AbstractArray{T},
                         zenith_opacity::AbstractArray{T},
                         zenith_angle::T) where T

Yields the 'temp' temperature altered by the atmosphere opacity at
'zenith_angle' for each element of 'temp' and associated 'zenith_opacity'.

"""
function atmos_opacity_impact(temp::T,
    zenith_opacity::T,
    zenith_angle::T) where T

    @assert 0 <= zenith_angle < 90 "zenith angle must be at or over 0 and below 90 \
            degrees"

    return temp .* exp(- zenith_opacity / cos(deg2rad(zenith_angle)))
end

function atmos_opacity_impact(temp::AbstractArray{T},
    zenith_opacity::T,
    zenith_angle::T) where T

    return atmos_opacity_impact.(temp, zenith_opacity, zenith_angle)
end

function atmos_opacity_impact(temp::AbstractArray{T},
    zenith_opacity::AbstractArray{T},
    zenith_angle::T) where T

    @assert size(temp) == size(zenith_opacity)
    
    return atmos_opacity_impact.(temp, zenith_opacity, zenith_angle)
end



"""
    atmosphere_model(alpha_grid::AbstractVector{T},
                     beta_grid::AbstractVector{T},
                     T_eff::T,
                     zenith_opacity::T) where T

Yields a 'SphereMap' structure representing an atmosphere model with constant
temperature 'T_eff' and opacity 'zenith_opacity'.

---
    atmosphere_model(T_eff::T,
                     zenith_opacity::T) where T

Yields a 'SphereMap' structure representing an atmosphere model with constant
temperature 'T_eff' and opacity 'zenith_opacity'. The default alpha and beta grids
are sampled at 1 degree resolution.

"""
function atmosphere_model(alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T},
    T_eff::T,
    zenith_opacity::T,
    T_bkg::T = zero(T)) where T

    atm_els = zeros(T, length(beta_grid))
    els_horizon = beta_grid .< 90
    atm_els[els_horizon] = atmos_opacity_impact.(T_bkg - T_eff, zenith_opacity, 
                                                 beta_grid[els_horizon]) .+ T_eff
    
    alpha_grid = alpha_grid[alpha_grid .!= T(360.)]
    # beta_grid = beta_grid[beta_grid .!= T(90.)]

    atm_map = reduce(vcat, [atm_els' for _ in 1:length(alpha_grid)])

    return SphereMap(alpha_grid, beta_grid, atm_map)
end

function atmosphere_model(T_eff::T,
    zenith_opacity::T,
    T_bkg::T = zero(T)) where T

    alpha_grid = [0., 180., 360.]
    beta_grid = collect(0.:1.:180.)

    return atmosphere_model(alpha_grid, beta_grid, T_eff, zenith_opacity, T_bkg)
end
