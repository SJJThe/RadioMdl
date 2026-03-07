
#TODO: see Chapt 9 S.Paine AM package for some refraction accounting

"""
    free_space_loss(rng::T,
                    freq::T) where T

Yields the free space loss for a given range and frequency.

"""
function free_space_loss(rng::T,
    freq::T) where T

    return ( 4*pi*rng / freq_to_wave(freq) )^2
end



"""
    simple_link_budget(gain_RX::T,
                       gain_TX::T,
                       rng::T,
                       freq::T,) where T

Yields the link budget coefficient between a receiver and transmitter according
to the Friis formula:

    gain_RX * gain_TX * (wavelength / (4*pi*rng))^2

"""
function simple_link_budget(gain_RX::T,
    gain_TX::T,
    rng::T,
    freq::T) where T

    L = free_space_loss(rng, freq)

    return gain_RX * 1/L * gain_TX
end



"""
    classic_gain_link_budget(sat_coord::SphereCoord{T},
                             sat_instru::Instrument{T},
                             tel_pointing_coord::SphereCoord{T},
                             tel_instru::Instrument{T};
                             pre_load_rot_mat::Union{Matrix,Nothing} = nothing,
                             simple_approx::Bool = false,
                             beam_avoid_angle::T = 0.0,
                             turn_off::Bool = false) where T

Yields the link budget coefficient between a satellite and a telescope given
their coordinates and instruments. The link budget is computed using the Friis
formula, accounting for the gains of the satellite and telescope antennas, the
distance between them and the frequency of observation.

It is possible to give a pre-computed rotation matrix to transform from
topocentric frame to telescope's antenna frame.

If 'simple_approx' is true, the satellite coordinates are approximated in the
telescope frame by negating the azimuthal angle and keeping the polar angle
unchanged. 

#FIXME:NOT YET IMPLEMENTED
If it is 'false', the satellite coordinates are transformed from the
topocentric frame to the telescope's antenna frame by passing via an
intermediate Earth-Centered Earth-Fixed (ECEF) frame.

If 'beam_avoid_angle' is greater than 0, the function accounts for the effect of
beam avoidance by the satellite. If the satellite boresight is closer than
'beam_avoid_angle' to the telescope pointing direction, the satellite gain is
reduced by "steering" away the satellite boresight of 45 degrees. If 'turn_off'
is true, the satellite gain is set to zero instead of being steered away.

"""
function classic_gain_link_budget(sat_coord::SphereCoord{T},
    sat_instru::Instrument{T},
    tel_pointing_coord::SphereCoord{T},
    tel_antenna::Antenna{T};
    pre_load_rot_mat::Union{Matrix,Nothing} = nothing,
    simple_approx::Bool = false,
    beam_avoid_angle::T = 0.0,
    turn_off::Bool = false) where T

    # coordinate of sat in telescope frame
    sat_coord_in_tel = pass_frame_to_frame(sat_coord, tel_pointing_coord; 
                                           pre_load_rot_mat=pre_load_rot_mat)

    # telescope gain
    gain_tel = get_gain_value(tel_antenna, sat_coord_in_tel)

    # coordinate of antenna from topocentric frame to satellite frame
    if simple_approx # if sat is close to zenith
        tel_coord_in_sat = SphereCoord(-tel_pointing_coord.alpha, 
                                       tel_pointing_coord.beta, sat_coord.r)
    else
        @error "3D transform handling is not yet implemented"#tel_coord_in_sat = #FIXME: pass to ECEF as intermediate
    end

    # beam avoidance effect
    if beam_avoid_angle > 0.
        # # get boresight pointing of telescope antenna
        # tel_beam_alpha, tel_beam_beta = get_boresight_gain(tel_instru.antenna)[2:3]

        # get boresight pointing of satellite antenna
        sat_beam_alpha, sat_beam_beta = get_boresight_gain(sat_instru.antenna)[2:3]

        # transform satellite beam coords in topocentric frame
        if simple_approx
            sat_beam_coord_topo = SphereCoord(-sat_beam_alpha, sat_beam_beta, 1.)
        else
            @error "3D transform handling is not yet implemented"#sat_beam_coord_topo = #FIXME: pass to ECEF as intermediate
        end

        # co-azimuthal angles are closer than beam_avoid_angle (satellite close
        # to telescope boresight)
        if abs(sat_beam_coord_topo.alpha - tel_pointing_coord.alpha) < beam_avoid_angle
            if turn_off
                return zero(T)
            else
                tel_coord_in_sat.alpha = mod(sat_beam_alpha + 45., 360.)
            end
        # polar angles are closer than beam_avoid_angle (satellite close to
        # telescope broesight)
        elseif abs(sat_beam_coord_topo.beta - tel_pointing_coord.beta) < beam_avoid_angle
            if turn_off
                return zero(T)
            else
                tel_coord_in_sat.beta = mod(sat_beam_beta + 45., 180.)
            end
        end
    end

    # satellite gain
    gain_sat = get_gain_value(sat_instru.antenna, tel_coord_in_sat)

    #link budget
    freq_bins = freq_range(sat_instru.receiver)
    link_budget_coefs = [simple_link_budget(gain_tel, gain_sat, sat_coord.r, f) 
                         for f in freq_bins]

    return link_budget_coefs
end
