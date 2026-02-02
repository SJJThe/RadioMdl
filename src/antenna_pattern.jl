
"""
    radiated_power_to_gain(rad_pow::AbstractDataFrame,
                           eta_rad::Real = 1.0;
                           alpha_col::Symbol = :az,
                           beta_col::Symbol = :polar,
                           map_col::Symbol = :power) where T

Yields the gain pattern of an antenna, given a radiated power pattern.

"""

function radiated_power_to_gain(rad_pow::AbstractDataFrame,
    eta_rad::Real = 1.0;
    alpha_col::Symbol = :az,
    beta_col::Symbol = :polar,
    map_col::Symbol = :power)

    @assert 0. <= eta_rad <= 1.

    # map the radiated power for interpolation
    rad_pow_map, a, b = map_sphere_coords(rad_pow; alpha_col=alpha_col, 
                                          beta_col=beta_col, map_col=map_col)

    # integrate over the sphere
    rad_pow_avg = trapz((a, b), rad_pow_map .* sin.(b')) / (4π)

    # directivity
    dir = rad_pow[:,map_col] ./ rad_pow_avg

    # gain
    return eta_rad .* dir
end



"""
    gain_to_effective_aperture(gain::Real,
                               wavelength::Real)

Yields the effective aperture of an antenna given its gain and wavelength.

"""
function gain_to_effective_aperture(gain::Real,
    wavelength::Real)
    
    return gain * (wavelength^2/(4π))
end



"""
    estim_hpbws(diameter::T,
                wavelength::T) where T

Yields the half power beamwidth (in degrees) of an antenna given its diameter
and wavelength.

"""
function estim_hpbws(diameter::T,
    wavelength::T) where T

    return rad2deg(67.6 * (wavelength / diameter))
end



"""
    effective_aperture_to_gain(effective_aperture::Real,
                               wavelength::Real)

Yields the gain of an antenna given its effective aperture and wavelength.

"""
function get_geometric_effective_aperture(aperture_efficiency::T,
    diameter::T) where T

    @assert T(0) <= aperture_efficiency <= T(1)

    return aperture_efficiency * pi * (diameter / 2)^2
end



"""
    antenna_mdl_ITU(gain_max::T,
                    half_beamwidth::T,
                    az::AbstractVector{T},
                    polar::AbstractVector{T};
                    single_rfi::Bool = false) where T

Create ITU recommended gain profile.

"""
function antenna_mdl_ITU(gain_max::T,
    half_beamwidth::T,
    az::AbstractVector{T},
    polar::AbstractVector{T};
    single_rfi::Bool = false) where T

    # gain profile container
    gain_profile = zeros(length(az))

    # select different parts of the gain profile
    gain_max_dB = 10. *log10(gain_max)
    parts = [0, half_beamwidth*sqrt(17/3), 10^((49-gain_max_dB)/25), 48, 80, 120, 180]
    part1 = findall(i -> parts[1] <= i < parts[2], az)
    part2 = findall(i -> parts[2] <= i < parts[3], az)
    part3 = findall(i -> parts[3] <= i < parts[4], az)
    part4 = findall(i -> parts[4] <= i < parts[5], az)
    part5 = findall(i -> parts[5] <= i < parts[6], az)
    part6 = findall(i -> parts[6] <= i <= parts[7], az)

    # calculate gain profile
    gain_profile[part1] .= gain_max_dB .- 3*(az[part1]./half_beamwidth).^2
    gain_profile[part2] .= gain_max_dB - (single_rfi ? 17 : 20)
    gain_profile[part3] .= (single_rfi ? 32 : 29) .- 25 .*log10.(az[part3])
    gain_profile[part4] .= (single_rfi ? -10 : -13)
    gain_profile[part5] .= (single_rfi ? -5 : -8)
    gain_profile[part6] .= (single_rfi ? -10 : -13)
    
    # create gain dataframe
    gain_pat = DataFrame(az=zeros(length(az)*length(polar)), 
                         polar=zeros(length(az)*length(polar)), 
                         gains=zeros(length(az)*length(polar)))
    for b in eachindex(polar)
        gain_pat[((b-1)*length(az)+1):b*length(az), :az] .= az
        gain_pat[((b-1)*length(az)+1):b*length(az), :polar] .= polar[b]
        gain_pat[((b-1)*length(az)+1):b*length(az), :gains] .= 10 .^(gain_profile./10)
    end

    return gain_pat
end



"""
    antenna_mdl_cst(gain::T,
                     az::AbstractVector{T},
                     polar::AbstractVector{T}) where T

Create constant gain pattern for omni-directional antennas.

"""
function antenna_mdl_cst(gain::T,
    az::AbstractVector{T},
    polar::AbstractVector{T}) where T

    gain_pat = DataFrame(az=zeros(length(az)*length(polar)), 
                         polar=zeros(length(az)*length(polar)), 
                         gains=zeros(length(az)*length(polar)))
    for b in eachindex(polar)
        gain_pat[((b-1)*length(az)+1):b*length(az), :az] .= az
        gain_pat[((b-1)*length(az)+1):b*length(az), :polar] .= polar[b]
        gain_pat[((b-1)*length(az)+1):b*length(az), :gains] .= gain
    end

    return gain_pat
end