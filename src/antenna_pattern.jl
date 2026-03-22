
"""
    integration_weights(alpha_grid::::AbstractVector{T},
                        beta_grid::AbstractVector{T}) where T

Yeilds a 2D-matrix of solid angle weights for the integration over a sphere. The
weights are defined by 'alpha_grid' and 'beta_grid'.

"""
function integration_weights(alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T}) where T

    dalpha = alpha_grid[2] - alpha_grid[1]
    @assert all(isapprox.(diff(alpha_grid), dalpha; atol=10.0^floor(Int, log10(dalpha))))
    dbeta = beta_grid[2] - beta_grid[1]
    @assert all(isapprox.(diff(beta_grid), dbeta; atol=10.0^floor(Int, log10(dbeta))))

    # convert to radians
    alpha_rads = deg2rad.(alpha_grid)
    beta_rads = deg2rad.(beta_grid)
    
    # dimension weights
    w_alpha = similar(alpha_rads)
    for i in eachindex(alpha_rads)
        lo = i == 1 ? zero(T) : alpha_rads[i] - alpha_rads[i-1]
        hi = i == lastindex(alpha_rads) ? zero(T) : alpha_rads[i+1] - alpha_rads[i]
        w_alpha[i] = (lo + hi) / 2
    end
    w_beta = similar(beta_rads)
    for j in eachindex(beta_rads)
        lo = j == 1 ? zero(T) : beta_rads[j] - beta_rads[j-1]
        hi = j == lastindex(beta_rads) ? zero(T) : beta_rads[j+1] - beta_rads[j]
        w_beta[j] = (lo + hi) / 2
    end

    return w_alpha * w_beta'
end



"""
    radiated_power_to_gain!(rad_pow::AbstractDataFrame,
                           eta_rad::Real = 1.0;

# using DimensionalData                           alpha_col::Symbol = :caz,
                           beta_col::Symbol = :polar,
                           map_col::Symbol = :power) where T

Yields the gain pattern of an antenna, given a radiated power pattern.

"""
function radiated_power_to_gain!(rad_pow::AbstractDataFrame,
    eta_rad::Real = 1.0;
    alpha_col::Symbol = :caz,
    beta_col::Symbol = :polar,
    map_col::Symbol = :power)

    @assert 0. <= eta_rad <= 1.

    # map the radiated power for interpolation
    rad_pow_map, a, b = map_sphere_coords(rad_pow; alpha_col=alpha_col, 
                                          beta_col=beta_col, map_col=map_col)

    # integrate over the sphere
    rad_pow_avg = trapz((deg2rad.(a), deg2rad.(b)), rad_pow_map .* sind.(b')) / (4π)

    # directivity
    rad_pow[:,map_col] ./= rad_pow_avg

    # gain
    rad_pow[:,map_col] .*= eta_rad
    map_col == :gains ? nothing : rename!(rad_pow, map_col => :gains)

    return rad_pow
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
    antenna_mdl_ITU_SA_509_3(gain_max::T,
                             half_beamwidth::T,
                             caz::AbstractVector{T},
                             pol::AbstractVector{T};
                             single_rfi::Bool = false) where T

Create ITU recommended gain profile according to ITU-R SA.509-3 "Space research
earth station and radio astronomy reference antenna radiation pattern for use in
interference calculations, including coordination procedures, for frequencies
less than 30 GHz". 

"""
function antenna_mdl_ITU_SA_509_3(gain_max::T,
    half_beamwidth::T,
    caz::AbstractVector{T},
    pol::AbstractVector{T};
    single_rfi::Bool = false) where T

    # gain profile container
    gain_profile = zeros(length(pol))

    # select different parts of the gain profile
    gain_max_dB = 10. *log10(gain_max)
    parts = [0, half_beamwidth*sqrt(17/3), 10^((49-gain_max_dB)/25), 48, 80, 120, 180]
    part1 = findall(i -> parts[1] <= i < parts[2], pol)
    part2 = findall(i -> parts[2] <= i < parts[3], pol)
    part3 = findall(i -> parts[3] <= i < parts[4], pol)
    part4 = findall(i -> parts[4] <= i < parts[5], pol)
    part5 = findall(i -> parts[5] <= i < parts[6], pol)
    part6 = findall(i -> parts[6] <= i <= parts[7], pol)

    # calculate gain profile
    gain_profile[part1] .= gain_max_dB .- 3*(pol[part1]./half_beamwidth).^2
    gain_profile[part2] .= gain_max_dB - (single_rfi ? 17 : 20)
    gain_profile[part3] .= (single_rfi ? 32 : 29) .- 25 .*log10.(pol[part3])
    gain_profile[part4] .= (single_rfi ? -10 : -13)
    gain_profile[part5] .= (single_rfi ? -5 : -8)
    gain_profile[part6] .= (single_rfi ? -10 : -13)
    
    # create gain dataframe
    gain_pat = DataFrame(polar=zeros(length(pol)*length(caz)), 
                         caz=zeros(length(pol)*length(caz)), 
                         gains=zeros(length(pol)*length(caz)))
    for b in eachindex(caz)
        gain_pat[((b-1)*length(pol)+1):b*length(pol), :polar] .= pol
        gain_pat[((b-1)*length(pol)+1):b*length(pol), :caz] .= caz[b]
        gain_pat[((b-1)*length(pol)+1):b*length(pol), :gains] .= 10 .^(gain_profile./10)
    end

    return gain_pat
end



"""
    antenna_mdl_ITU_RA_1631(ant_diameter::T,
                            wavelength::T,
                            caz::AbstractVector{T},
                            pol::AbstractVector{T}) where T

Create ITU recommended gain profile according to ITU-R RA.1631-1 "Reference
radio astronomy antenna pattern to be used for compatibility analyses between
non-GSO systems and radio astronomy service stations based on the epfd concept".

"""
function antenna_mdl_ITU_RA_1631(gain_max::T,
    ant_diameter::T,
    wavelength::T,
    caz::AbstractVector{T},
    pol::AbstractVector{T}) where {T}
    
    # gain profile container
    gain_profile = zeros(length(pol))

    # select different parts of the gain profile
    parts = [0., 69.88/(ant_diameter / wavelength), 1., 10., 34.1, 80., 120., 180.]
    part1 = findall(i -> parts[1] < i < parts[2], pol)
    part2 = findall(i -> parts[2] <= i < parts[3], pol)
    part3 = findall(i -> parts[3] <= i < parts[4], pol)
    part4 = findall(i -> parts[4] <= i < parts[5], pol)
    part5 = findall(i -> parts[5] <= i < parts[6], pol)
    part6 = findall(i -> parts[6] <= i <= parts[7], pol)
    part7 = findall(i -> parts[7] <= i <= parts[8], pol)

    # calculate gain profile
    x1 = π * ant_diameter / (360. * wavelength) .* pol[part1]
    x2 = π * ant_diameter / (360. * wavelength) .* pol[part2]
    B = 10^3.2 * π^2 * ((π * ant_diameter / 2) / (180. * wavelength))^2
    gain_profile[1] = gain_max
    gain_profile[part1] .= gain_max .* (besselj1.(2π .* x1) ./ (π .* x1)).^2
    gain_profile[part2] .= B .* (cos.(2π .* x2 .- 3π/4 .+ .0953) ./ (π .* x2)).^2
    gain_profile[part3] .= 10. .^((29. .- 25. .* log10.(pol[part3])) ./ 10.)
    gain_profile[part4] .= 10. .^((34. .- 30. .* log10.(pol[part4])) ./ 10.)
    gain_profile[part5] .= 10. .^(-12. / 10.)
    gain_profile[part6] .= 10. .^(-7. / 10.)
    gain_profile[part7] .= 10. .^(-12. / 10.)

    # create gain dataframe
    gain_pat = DataFrame(polar=zeros(length(pol)*length(caz)), 
                         caz=zeros(length(pol)*length(caz)), 
                         gains=zeros(length(pol)*length(caz)))
    for b in eachindex(caz)
        gain_pat[((b-1)*length(pol)+1):b*length(pol), :polar] .= pol
        gain_pat[((b-1)*length(pol)+1):b*length(pol), :caz] .= caz[b]
        gain_pat[((b-1)*length(pol)+1):b*length(pol), :gains] .= gain_profile
    end
    
    return gain_pat
end



"""
    antenna_mdl_cst(gain::T,
                     caz::AbstractVector{T},
                     pol::AbstractVector{T}) where T

Create constant gain pattern for omni-directional antennas.

"""
function antenna_mdl_cst(gain::T,
    caz::AbstractVector{T},
    pol::AbstractVector{T}) where T

    gain_pat = DataFrame(caz=zeros(length(caz)*length(pol)), 
                         polar=zeros(length(caz)*length(pol)), 
                         gains=zeros(length(caz)*length(pol)))
    for b in eachindex(pol)
        gain_pat[((b-1)*length(caz)+1):b*length(caz), :caz] .= caz
        gain_pat[((b-1)*length(caz)+1):b*length(caz), :polar] .= pol[b]
        gain_pat[((b-1)*length(caz)+1):b*length(caz), :gains] .= gain
    end

    return gain_pat
end