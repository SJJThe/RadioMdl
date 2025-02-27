
"""
assumes alphas and betas are degrees. Output is converted in radians.
"""
function map_sphere(pattern::AbstractArray{T},
    alphas::AbstractArray{T},
    betas::AbstractArray{T}) where T
    
    # form 2D matrix for interpolation argument
    # add first column as last column to loop azimuth coordinates
    gain_map = hcat(reshape(pattern, length(alphas), length(betas)),
                    pattern[1:length(alphas)])
    
    # generate sampling coordinates
    a = alphas .* rad
    b = [betas; T(360)] .* rad

    return gain_map, a, b
end


"""
yields the gain pattern of an antenna, in dB, given a radiated power pattern. It
is assumed that the radiated power includes the radiation efficiency.
the angles must be in degrees

"""
function radiated_power_to_gain(rad_pow::AbstractVector{T},
    alphas::AbstractVector{T},
    betas::AbstractVector{T};
    eta_rad::Real = 1.0) where T
    
    # map the radiated power for interpolation
    rad_pow_map, a, b = map_sphere(rad_pow, alphas, betas)
    
    # integrate over the sphere
    rad_pow_avg = trapz((a, b), rad_pow_map .* sin.(a)) / (4π)

    # directivity
    dir = rad_pow ./ rad_pow_avg

    # gain
    return eta_rad .* dir
end



"""
"""
function interpolate_gain(gain::AbstractArray{T},
    alphas::AbstractArray{T},
    betas::AbstractArray{T}) where T
    
    # map the gain for interpolation
    gain_map, a, b = map_sphere(gain, alphas, betas)

    # gain function of angles in antenna coord. system
    return interpolate((a, b), gain_map, Gridded(Linear()))
end



"""
"""
function gain_to_effective_aperture(gain::Real,
    frequency::Real)
    
    wavelength = speed_c / frequency
    return gain * (wavelength^2/(4π))
end



"""
create ITU recommended gain profile
"""
function antenna_mdl_ITU(gain_max::T,
    half_beamwidth::T,
    alphas::AbstractVector{T},
    betas::AbstractVector{T}) where T

    # gain profile container
    gain_profile = zeros(length(alphas))

    # select different parts of the gain profile
    parts = [0, half_beamwidth*sqrt(17/3), 10^((49-gain_max)/25), 48, 80, 120, 180]
    part1 = findall(i -> parts[1] <= i < parts[2], alphas)
    part2 = findall(i -> parts[2] <= i < parts[3], alphas)
    part3 = findall(i -> parts[3] <= i < parts[4], alphas)
    part4 = findall(i -> parts[4] <= i < parts[5], alphas)
    part5 = findall(i -> parts[5] <= i < parts[6], alphas)
    part6 = findall(i -> parts[6] <= i <= parts[7], alphas)

    # calculate gain profile
    gain_profile[part1] .= gain_max .- 3*(alphas[part1]./half_beamwidth).^2
    gain_profile[part2] .= gain_max - 20
    gain_profile[part3] .= 29 .- 25 .*log10.(alphas[part3])
    gain_profile[part4] .= -13
    gain_profile[part5] .= -8
    gain_profile[part6] .= -13
    
    # create gain dataframe
    gain_pat = DataFrame(alphas=zeros(length(alphas)*length(betas)), 
                         betas=zeros(length(alphas)*length(betas)), 
                         gains=zeros(length(alphas)*length(betas)))
    for b in eachindex(betas)
        gain_pat[((b-1)*length(alphas)+1):b*length(alphas), :alphas] .= alphas
        gain_pat[((b-1)*length(alphas)+1):b*length(alphas), :betas] .= betas[b]
        gain_pat[((b-1)*length(alphas)+1):b*length(alphas), :gains] .= 10 .^(gain_profile./10)
    end

    return gain_pat
end

