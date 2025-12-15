

"""
alphas are angles from z-axis towards x-y-plane, betas are angles from x-axis
towards y-axis, all in degree
the gain interpolator coords are in radians#FIXME: change to dgrees for
coherence
T_phy can be a scalar (constant physical temperature) or a DataFrame/YAXArray
containing time stamps and physical temperatures (for time varying T_phy)
"""
struct Antenna{T<:AbstractFloat,I<:Interpolations.GriddedInterpolation}
    ant_diameter::T # antenna diameters in m
    gain_pat::AbstractDataFrame # gain pattern
    gain_func::I # gain interpolator
    ap_eff::T # aperture efficiency
    rad_eff::T # radiation efficiency
    T_phy::Union{T,AbstractDataFrame,YAXArray{T}} # physical temperature in K
    valid_freqs::Tuple{T,T} # min and max valid frequencies for the gain model

    function Antenna(ant_diameter::T, 
        gain_pat::AbstractDataFrame,
        gain_func::I,
        ap_eff::T = .5,
        rad_eff::T = .9,
        T_phy::Union{T,AbstractDataFrame,YAXArray{T}} = 300.,
        valid_freqs::Tuple{T,T} = (zero(T),T(1e12))) where {T<:AbstractFloat,
        I<:Interpolations.GriddedInterpolation}
           
        @assert ant_diameter > 0
        for n in propertynames(gain_pat)
            @assert n in [:alphas, :betas, :gains]
        end
        @assert 0 <= rad_eff <= 1
        @assert 0 <= ap_eff <= 1
        if typeof(T_phy) <: T
            @assert T_phy > 0
        elseif typeof(T_phy) <: AbstractDataFrame
            @assert :time in propertynames(T_phy) "`:time` must be a column of T_phy"
            @assert :temp in propertynames(T_phy) "`:temp` must be a column of T_phy"
            @assert minimum(T_phy[:temp]) > 0
        elseif typeof(T_phy) <: YAXArray
            @assert hasdim(T_phy, :time) "`:time` must be a dimension of T_phy"
            @assert minimum(T_phy) > 0
        end
        @assert valid_freqs[1] < valid_freqs[2]
        
        return new{T,I}(ant_diameter, gain_pat, gain_func, rad_eff, ap_eff, T_phy, valid_freqs)
    end
end

function Antenna(ant_diameter::T, 
    gain_ant::AbstractDataFrame,
    ap_eff::T = .5,
    rad_eff::T = .9,
    T_phy::Union{T,AbstractDataFrame,YAXArray{T}} = 300.,
    valid_freqs::Tuple{T,T} = (zero(T),zero(T))) where T

    # create the gain interpolator
    alphas = subset(gain_ant, :betas => b -> b .== gain_ant[1,:betas];
                    view=true)[!,:alphas]
    betas = subset(gain_ant, :alphas => a -> a .== gain_ant[1,:alphas];
                   view=true)[!,:betas]
    gain_func = interpolate_gain(gain_ant[!,:gains], alphas, betas)

    return Antenna(ant_diameter, gain_ant, gain_func, ap_eff, rad_eff, T_phy, valid_freqs)
end

function Antenna(file_pattern_path::String,
    ant_diameter::T,
    ap_eff::T = .5,
    rad_eff::T = .9,
    T_phy::Union{T,AbstractDataFrame,YAXArray{T}} = 300.,
    valid_freqs::Tuple{T,T} = (zero(T),zero(T));
    power_tag::Symbol = :gains,
    declination_tag::Symbol = :alphas,
    azimuth_tag::Symbol = :betas) where T
    
    # load the antenna power pattern
    @assert occursin(".cut", file_pattern_path) "the power pattern file must be a .cut \
    file"
    gain_ant = power_pattern_from_cut_file(file_pattern_path)

    # rename angles columns
    rename!(gain_ant, declination_tag => :alphas)
    rename!(gain_ant, azimuth_tag => :betas)
    
    # convert into gain (do not use unique here for if different alphas for two betas)
    alphas = subset(gain_ant, :betas => b -> b .== gain_ant[1,:betas]; 
                    view=true)[!,:alphas]
    betas = subset(gain_ant, :alphas => a -> a .== gain_ant[1,:alphas]; 
                   view=true)[!,:betas]
    gain_ant[:,power_tag] = radiated_power_to_gain(gain_ant[!,power_tag], alphas, betas;
                                                   eta_rad=rad_eff)
    rename!(gain_ant, power_tag => :gains)

    return Antenna(ant_diameter, gain_ant, ap_eff, rad_eff, T_phy, valid_freqs)
end

get_antenna_diameter(a::Antenna) = a.ant_diameter
get_gain_pattern(a::Antenna) = a.gain_pat
get_gain_value(a::Antenna, alpha::Real, beta::Real) = a.gain_func(alpha, beta)
function get_directivity_value(a::Antenna, alpha::Real, beta::Real)
    return get_gain_value(a, alpha, beta) / get_rad_eff(a)
end
function get_def_angles(a::Antenna)
    return unique(get_gain_pattern(a)[!,:alphas]), unique(get_gain_pattern(a)[!,:betas])
end
get_boresight_gain(a::Antenna) = maximum(get_gain_pattern(a)[:,:gains])
function get_boresight_point(a::Antenna) 
    gain = get_gain_pattern(a)
    i = findmax(gain[:,:gains])[2]
    return gain[i,:alphas], gain[i,:betas]
end
function get_slice_gain(a::Antenna,
    beta::Real)
    gain_pat = get_gain_pattern(a)
    g_pos = subset(gain_pat, :betas => b -> b .== beta, view=true)
    g_neg = subset(gain_pat, :betas => b -> b .== beta + 180, view=true)
    alphas = [-reverse(g_neg[1:end-1,:alphas]); g_pos[2:end,:alphas]]
    gains = [reverse(g_neg[1:end-1,:gains]); g_pos[2:end,:gains]]
    return alphas, gains
end
get_ap_eff(a::Antenna) = a.ap_eff
get_rad_eff(a::Antenna) = a.rad_eff
get_T_phy(a::Antenna) = a.T_phy
get_valid_freqs(a::Antenna) = a.valid_freqs
get_hpbws(a::Antenna, wavelength::Real) = estim_hpbws(a.ant_diameter, wavelength)
function get_geometric_effective_aperture(a::Antenna)
    return get_geometric_effective_aperture(a.ap_eff, a.ant_diameter)
end
#FIXME: do the real convolution
get_antenna_temperature(a::Antenna, T_b) = get_boresight_gain(a) .* T_b ./ (4*pi)
get_antenna_radiation_loss(a::Antenna) = (1 - get_rad_eff(a)) .* get_T_phy(a)



"""
"""
struct Receiver{T<:AbstractFloat}
    freq_res::T # frequency resolution
    cent_freq::T # center frequency
    bw::T # bandwidth
    gain_amps::T # gain of amplifiers
    # receiver temperature (scalar, freq vector or freq-time matrix)
    # can includes quantum noise, LNA, etc.
    T_rx::Union{T,AbstractDataFrame,YAXArray{T}}
    freq_resp::Union{AbstractArray{T},AbstractDataFrame,YAXArray{T}} # frequency response

    function Receiver(freq_res::T,
        cent_freq::T,
        bw::T,
        gain_amps::T, 
        T_rx::Union{T,AbstractDataFrame,YAXArray{T}},
        freq_resp::Union{AbstractArray{T},AbstractDataFrame,YAXArray{T}}) where {T<:AbstractFloat}

        @assert freq_res > 0
        @assert cent_freq > 0
        @assert bw > 0
        @assert gain_amps > 0
        if typeof(T_rx) <: T
            @assert T_rx > 0
        elseif typeof(T_rx) <: AbstractDataFrame
            @assert freq_resp isa AbstractDataFrame "must work in DataFrame"
            @assert :temp in propertynames(T_rx) "`:temp` must be a column of T_rx"
            @assert :time in propertynames(T_rx) || :freq in propertynames(T_rx) "`:time` or `:freq` must be a column of T_rx"
            @assert minimum(T_rx[:temp]) > 0
            @assert :freq_resp in propertynames(freq_resp) "`:freq_resp` must be a column of freq_resp"
            @assert :freq in propertynames(freq_resp) "`:freq` must be a column of freq_resp"
            @assert (sum(freq_resp[:freq_resp].^2) - 1) <= 1e-10 "The receiver frequency response must be normalized"
            if :freq in propertynames(T_rx)
                @assert freq_resp[:freq] == T_rx[:freq] "The frequencies in freq_resp and T_rx must match"
            end
        elseif typeof(T_rx) <: YAXArray
            @assert freq_resp isa YAXArray "must work in YAXArray"
            @assert hasdim(T_rx, :time) || hasdim(T_rx, :freq) "`:time` or `:freq` must be a dimension of T_rx"
            @assert minimum(T_rx) > 0
            @assert hasdim(freq_resp, :freq) "`:freq` must be a dimension of freq_resp"
            @assert (sum(freq_resp[:freq_resp].^2) - 1) <= 1e-10 "The receiver frequency response must be normalized"
            if hasdim(T_rx, :freq)
                @assert freq_resp[:freq] == T_rx[:freq] "The frequencies in freq_resp and T_rx must match"
            end
        end

        return new{T}(freq_res, cent_freq, bw, gain_amps, T_rx, freq_resp)
    end
end

function Receiver(freq_res::T,
    cent_freq::T,
    bw::T,
    gain_amps::T, 
    T_rx::Union{T,AbstractDataFrame,YAXArray{T}}) where {T<:AbstractFloat}

    # flat frequency response
    nb_freq_chan = div(bw, freq_res)
    freq_resp = T.(ones(nb_freq_chan) ./ nb_freq_chan)

    if size(T_rx,2) > 1
        @assert size(T_rx,2) == nb_freq_chan "The second dimension of T_rx must match \
        the number of frequency channels defined by bw and freq_res"
    end

    return Receiver(freq_res, cent_freq, bw, gain_amps, T_rx, freq_resp)
end

get_freq_res(r::Receiver) = r.freq_res
get_center_freq(r::Receiver) = r.cent_freq
get_bw(r::Receiver) = r.bw
get_gain_amps(r::Receiver) = r.gain_amps
get_T_rx(r::Receiver) = r.T_rx
get_freq_resp(r::Receiver) = r.freq_resp
get_nb_freq_chan(r::Receiver) = div(r.bw, r.freq_res)
function get_center_freq_chans(r::Receiver)
    freq_chan = get_nb_freq_chan(r)
    bw_RX = get_bandwidth(r)
    rng_freq = range(-bw_RX/2 + delta_freq/2, bw_RX/2 - delta_freq/2, length=freq_chan)
    return get_center_freq(i) .+ rng_freq
end



"""
suppose the frame of antenna is oriented North-West-Up. Assumes the output of
signal_func is expressed in Kelvin (temperature).
"""
struct Instrument{T<:AbstractFloat}
    antenna::Antenna # antenna
    receiver::Receiver{T} # receiver of precision T
    coords::Dict{Symbol, T} # coordinates

    function Instrument(antenna::Antenna,
        receiver::Receiver{T},
        coords::Dict{Symbol, T} = Dict(:lat => 0., 
                                       :lon => 0.)) where {T<:AbstractFloat}
        
        ant_fmin, ant_fmax = get_valid_freqs(antenna)
        cent_freq = get_center_freq(receiver)
        bw = get_bw(receiver)
        @assert (ant_fmin <= cent_freq - bw/2) && (cent_freq + bw/2 <= ant_fmax)

        if antenna.T_phy isa AbstractDataFrame && receiver.T_rx isa AbstractDataFrame
            if :time in propertynames(receiver.T_rx)
                @assert receiver.T_rx[:time] == antenna.T_phy[:time] "the time stamps of T_rx and T_phy must match"
            end
        elseif antenna.T_phy isa YAXArray && receiver.T_rx isa YAXArray
            if hasdim(receiver.T_rx, :time)
                @assert receiver.T_rx[:time] == antenna.T_phy[:time] "the time stamps of T_rx and T_phy must match"
            end
        end

        @assert :lat in keys(coords) "`:lat` must be a key of coords"
        @assert :lon in keys(coords) "`:lon` must be a key of coords"

        return new{T}(antenna, receiver, coords)
    end
end

get_antenna(i::Instrument) = i.antenna
get_receiver(i::Instrument) = i.receiver
get_coords(i::Instrument) = i.coords

"""
"""
function instrument_psd_stat(i::Instrument{T},
    T_b::Union{T,AbstractDataFrame,YAXArray{T}},# sky brightness temperature
    integration_samp::Real = 1) where T

    # get receiver and antenna parameters
    rec = i.receiver
    ant = i.antenna

    # instrument gain coefficient
    gain = rec.gain_amps * rec.freq_res * k_boltz * impedance

    # antenna temperature
    T_a = get_antenna_temperature(ant, T_b)
    
    # instrument noise temperature
    T_n = rec.T_rx .+ get_antenna_radiation_loss(ant)

    # calculate power spectral density
    return intrument_psd_stat(gain, T_a, T_n, integration_samp)
end



"""
traj is a DataFrame where each row is indexed by a datatime and has an elevation
and azimuth angle(s). It is possible to have vectors of els and azs for a same
time (e.g. for sky mapping). All measurements are assumed to be given in SI
units (in degrees for angles, meters for distances).
"""
struct Trajectory
    traj::AbstractDataFrame # azimuth, elevation and distance info for each sampled time
    
    function Trajectory(traj::AbstractDataFrame)
        
        for n in propertynames(traj)
            @assert n in [:times, :azimuths, :elevations, :distances]
        end
        @assert length(unique(traj.times)) == length(traj.times)
        @assert typeof(traj.times[1]) == DateTime
        @assert minimum(length.(traj.azimuths) .== length(traj.azimuths[1]))
        @assert minimum(length.(traj.elevations) .== length(traj.elevations[1]))
        @assert minimum(length.(traj.distances) .== length(traj.distances[1]))
        @assert length(traj.azimuths[1]) == length(traj.elevations[1])
        @assert length(traj.azimuths[1]) == length(traj.distances[1])

        sort!(traj, :times)

        return new(traj)
    end
end
#TODO: add function to load trajectory from DataFrame with different tags
function Trajectory(file_path::String;
    time_tag::Symbol = :times,
    elevation_tag::Symbol = :altitudes,
    azimuth_tag::Symbol = :azimuths,
    distance_tag::Symbol = :distances)

    # load the trajectory as a DataFrame
    #FIXME: read also csv files without dates for real antenna positions
    if occursin(".arrow", file_path)
        traj = DataFrame(Arrow.Table(file_path))
    elseif occursin(".csv", file_path)
        traj = DataFrame(CSV.File(file_path))
    else
        error("the trajectory points are not in Arrow or CSV format")
    end
    
    # rename columns
    rename!(traj, time_tag => :times)
    rename!(traj, azimuth_tag => :azimuths)
    rename!(traj, elevation_tag => :elevations)
    rename!(traj, distance_tag => :distances)
    
    # convert time stamps into DateTime type
    @. traj[!,:times] = Dates.DateTime(traj[!,:times])

    return Trajectory(traj[!, [:times, :azimuths, :elevations, :distances]])
end

get_traj(t::Trajectory) = t.traj
function get_traj(t::Trajectory, t0::DateTime, t1::DateTime;
    skipmissing::Bool = true,
    view::Bool = true)
    
    return subset(get_traj(t), :times => t -> t0 .<= t .<= t1; 
                  skipmissing=skipmissing, view=view)
end
get_time_bounds(t::Trajectory) = (t.traj[1,:times], t.traj[end,:times])
get_time_stamps(t::Trajectory) = t.traj[:,:times]
get_azimuths(t::Trajectory) = t.traj[:,:azimuths]
get_elevations(t::Trajectory) = t.traj[:,:elevations]
get_distances(t::Trajectory) = t.traj[:,:distances]



"""
pts can contain different positions for a same time, e.g. for a sky map.
"""
struct Observation{T<:AbstractFloat}
    pts::Trajectory # trajectory of the observation
    inst::Instrument{T} # instrument used for observation
    result::AbstractArray{T} # store the results of the modeling of the observation#FIXME:DataFrame?
    #FIXME: merge with pts dataframe to get power AND position for maps...
    function Observation(pts::Trajectory,
        inst::Instrument{T},
        result::AbstractArray{T}) where T
        
        @assert length(get_time_stamps(pts)) == size(result, 1)
        @assert length(get_azimuths(pts)[1]) == size(result, 2)
        @assert get_nb_freq_chan(inst) == size(result, 3)

        return new{T}(pts, inst, result)
    end
end

function Observation(start_date::DateTime,
    stop_date::DateTime,
    trajectory::Trajectory,
    instrument::Instrument{T};
    filt_funcs::NTuple{N,Pair} = ()) where {N,T}
    
    # filter date and other from trajectory
    traj = get_traj(trajectory, start_date, stop_date; view=false)
    for filt in filt_funcs
        traj = subset(traj, filt; skipmissing=true, view=true)
    end
    isempty(traj) && error("No pointing positions found for the given time window and \
                            custom filters.")
    sort!(traj, :times)

    pts = Trajectory(traj)

    # create result storage
    len_time = length(get_time_stamps(pts))
    len_pos = length(get_azimuths(pts)[1])
    len_freq = get_nb_freq_chan(instrument)
    result = fill!(zeros(T, len_time, len_pos, len_freq), NaN)

    return Observation(pts, instrument, result)
end

get_traj(obs::Observation) = get_traj(obs.pts)
get_time_bounds(obs::Observation) = get_time_bounds(obs.pts)
get_time_stamps(obs::Observation; kwds...) = get_time_stamps(obs.pts; kwds...)
get_azimuths(obs::Observation) = get_azimuths(obs.pts)
get_elevations(obs::Observation) = get_elevations(obs.pts)
get_distances(obs::Observation) = get_distances(obs.pts)
get_instrument(obs::Observation) = obs.inst
get_result(obs::Observation) = obs.result

function estim_temp(flux::Real,
    obs::Observation)

    instru = get_instrument(obs)
    frequency = get_center_freq(instru)
    ant = get_antenna(instru)
    max_gain = get_boresight_gain(ant)
    A_eff_max = gain_to_effective_aperture(max_gain, freq_to_wave(frequency))

    return estim_temp(flux, A_eff_max)
end




"""
assumes the positions of sats are time-synced with the time samples of
observation.
Suppose the frame of satellite antenna is oriented North-East-Nadir. The antenna
pointing can be any direction from Nadir, defined in the map gain.
"""
struct Constellation{T<:AbstractFloat}
    sats::AbstractDataFrame
    tmt::Instrument{T}
    lnk_bdgt_mdl::Function

    function Constellation(sats::AbstractDataFrame,
        tmt::Instrument{T},
        lnk_bdgt_mdl::Function) where T

        # check lnk_bdgt_mdl signature is correct
        @assert hasmethod(lnk_bdgt_mdl, (T, T, Instrument{T}, T, T, T, Instrument{T}, T))

        return new{T}(sats, tmt, lnk_bdgt_mdl)
    end
end

function Constellation(sats::AbstractDataFrame,
    observation::Observation,
    sat_tmt::Instrument{T},
    lnk_bdgt_mdl::Function = sat_link_budget;
    filt_funcs::NTuple{N,Pair} = ()) where {N,T}

    # observation window
    start_date , stop_date = get_time_bounds(observation)
    
    # apply the custom filters
    sats = subset(sats, :times => t -> start_date .<= t .<= stop_date; 
                  skipmissing=true, view=true)
    for filt in filt_funcs
        sats = subset(sats, filt; skipmissing=true, view=true)
    end

    # check = minimum(get_time_stamps(observation) .== unique(sats[!,:times]))
    #=@assert check=# @warn "Observation time stamps and Constellation time stamps needs\
                              to be aligned."

    return Constellation(sats, sat_tmt, lnk_bdgt_mdl)
end

function Constellation(file_path::String, 
    observation::Observation,
    sat_tmt::Instrument{T},
    lnk_bdgt_mdl::Function = sat_link_budget;
    name_tag::Symbol = :sat,
    time_tag::Symbol = :time_stamps,
    elevation_tag::Symbol = :altitudes,
    azimuth_tag::Symbol = :azimuths,
    distance_tag::Symbol = :distances,
    filt_funcs::NTuple{N,Pair} = ()) where {N,T}

    # load the trajectory as a DataFrame
    sats = DataFrame(Arrow.columntable(Arrow.Table(file_path)))

    # rename columns
    rename!(sats, time_tag => :times)
    rename!(sats, name_tag => :sat)
    rename!(sats, azimuth_tag => :azimuths)
    rename!(sats, elevation_tag => :elevations)
    rename!(sats, distance_tag => :distances)
    
    @. sats[!,:times] = Dates.DateTime(sats[!,:times])

    sort!(sats, :times)

    return Constellation(sats, observation, sat_tmt, lnk_bdgt_mdl; 
                         filt_funcs=filt_funcs)
end

get_antenna(c::Constellation) = get_antenna(c.tmt)
get_transmitter(c::Constellation) = c.tmt
get_sats_name(c::Constellation) = unique(c.sats[!,:sat])
get_lnk_bdgt_mdl(c::Constellation) = c.lnk_bdgt_mdl

function get_sat_traj(c::Constellation,
    s::String)

    return subset(c.sats, :sat => n -> n .== s; view=true)
end

function get_sats_names_at_time(c::Constellation,
    t::DateTime)
    
    return subset(c.sats, :times => ts -> ts .== t; view=true)[!,:sat]
end
