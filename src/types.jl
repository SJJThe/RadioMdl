
"""
    Antenna(ant_diameter::T,
            gain_pat::SphereMap{T},
            ap_eff::T,
            rad_eff::T,
            T_phy::Union{T,AbstractVector{T},DimArray{T}},
            valid_freqs::Tuple{<:Real,<:Real}) where T<:AbstractFloat
        
Yields an 'Antenna' structure that defines the antenna properties, that is its
diameter 'ant_diameter', gain pattern 'gain_pat' as a 'SphereMap' ( see
['SphereMap'](@ref) ), aperture efficiency 'ap_eff', radiation efficiency
'rad_eff', physical temperature 'T_phy' as a constant, a vector or a 'DimArray'
and valid frequency range 'valid_freqs'.

It is possible to use different get functions:
'''
# returns the gain at (alpha, beta)
get_gain_value(A::Antenna, alpha::Real, beta::Real)

# returns the directivity value at (alpha, beta)
get_directivity_value(A::Antenna, alpha::Real, beta::Real) 

# returns the angle girds where the gain pattern is sampled
get_angle_grids(A::Antenna)

# returns the boresight gain
get_boresight_gain(A::Antenna)

# return the half power beamwidths
get_hpbws(A::Antenna, wavelength::Real)

# return the geometric effective aperture
get_geometric_effective_aperture(A::Antenna)

# return the antenna radiation loss
get_antenna_radiation_loss(A::Antenna)
'''

---
    Antenna(ant_diameter::T, 
            gain_pat::AbstractDataFrame,
            ap_eff::T,
            rad_eff::T,
            T_phy::Union{T,AbstractVector{T},DimArray{T}},
            valid_freqs::Tuple{<:Real,<:Real}) where T<:AbstractFloat

Yields an 'Antenna' structure based on the dataframe 'gain_pat'. The columns of
the dataframe 'gain_pat' must contain the columns `:az`, `:polar` and `:gains`.

---
    Antenna(cut_file_path::String,
            ant_diameter::T,
            ap_eff::T,
            rad_eff::T,
            T_phy::Union{T,AbstractVector{T},DimArray{T}},
            valid_freqs::Tuple{<:Real,<:Real}) where T<:AbstractFloat

Yields an 'Antenna' structure based on the cut file located at 'cut_file_path',
using the ['power_pattern_from_cut_file'](@ref) function.

"""
struct Antenna{T}
    ant_diameter::T # antenna diameters in m
    gain_pat::SphereMap{T} # gain pattern as a spherical map
    ap_eff::T # aperture efficiency
    rad_eff::T # radiation efficiency
    # physical temperature in K (constant or time dependent)
    T_phy::Union{T,AbstractVector{T},DimArray{T}}
    valid_freqs::Tuple{<:Real,<:Real} # min and max valid frequencies for the gain model

    function Antenna(ant_diameter::T, 
        gain_pat::SphereMap{T},
        ap_eff::T,
        rad_eff::T,
        T_phy::Union{T,AbstractVector{T},DimArray{T}},
        valid_freqs::Tuple{<:Real,<:Real}) where T

        @assert ant_diameter > 0
        @assert 0 <= rad_eff <= 1
        @assert 0 <= ap_eff <= 1
        @assert all(T_phy .>= 0)
        if typeof(T_phy) <: DimArray
            @assert dims(T_phy) == [:times] "the dimension of T_phy must be :times"
        end
        @assert valid_freqs[1] < valid_freqs[2]

        return new{T}(ant_diameter, gain_pat, ap_eff, rad_eff, T_phy, valid_freqs)
    end
end

function Antenna(ant_diameter::T, 
    gain_pat::AbstractDataFrame,
    ap_eff::T,
    rad_eff::T,
    T_phy::Union{T,AbstractVector{T},DimArray{T}},
    valid_freqs::Tuple{<:Real,<:Real}) where T

    # create the sphere map 
    @assert :az in propertynames(gain_pat) "the gain pattern must have a column \
                                                `:az`"
    @assert :polar in propertynames(gain_pat) "the gain pattern must have a column \
                                               `:polar`"
    @assert :gains in propertynames(gain_pat) "the gain pattern must have a column \
                                               `:gains`"
    SM = SphereMap(gain_pat; map_col=:gains)
    
    return Antenna(ant_diameter, SM, ap_eff, rad_eff, T_phy, valid_freqs)
end

function Antenna(cut_file_path::String,
    ant_diameter::T,
    ap_eff::T,
    rad_eff::T,
    T_phy::Union{T,AbstractVector{T},DimArray{T}},
    valid_freqs::Tuple{<:Real,<:Real}) where T

    # load the antenna power pattern
    pattern = power_pattern_from_cut_file(cut_file_path)

    # convert radiated power to gain
    pattern[:,:power] = radiated_power_to_gain(pattern, rad_eff)
    rename!(pattern, :power => :gains)

    return Antenna(ant_diameter, pattern, ap_eff, rad_eff, T_phy, valid_freqs)
end

get_gain_value(A::Antenna, alpha::Real, beta::Real) = A.gain_pat(alpha, beta)

function get_directivity_value(A::Antenna,
    alpha::Real,
    beta::Real)
    
    return get_gain_value(A, alpha, beta) / A.rad_eff
end

get_angle_grids(A::Antenna) = get_angle_grids(A.gain_pat)

function get_boresight_gain(A::Antenna)

    gain = A.gain_pat
    i = findmax(gain.map)[2]

    return gain.map[i], gain.alpha_grid[i[2]], gain.beta_grid[i[1]]
end

get_hpbws(A::Antenna, wavelength::Real) = estim_hpbws(A.ant_diameter, wavelength)

function get_geometric_effective_aperture(A::Antenna)
    return get_geometric_effective_aperture(A.ap_eff, A.ant_diameter)
end

get_antenna_radiation_loss(A::Antenna) = (1 - A.rad_eff) .* A.T_phy

Base.show(io::IO, A::Antenna{T}) where {T} = begin 
    print(io, "Antenna{$T}:\n")
    print(io, "diameter: $(A.ant_diameter)\n")
    print(io, "gain pattern: $(A.gain_pat)\n")
    print(io, "aperture efficiency: $(A.ap_eff)\n")
    print(io, "radiation efficiency: $(A.rad_eff)\n")
    print(io, "physical temperature: $(A.T_phy)\n")
    print(io, "valid frequency range: $(A.valid_freqs)")
end


"""
    Receiver(freq_res::T,
             cent_freq::T,
             bw::T,
             gain_amps::T,
             T_rx::Union{T,AbstractArray{T},DimArray{T}},
             freq_resp::AbstractVector{T}) where T

Yields a 'Receiver' structure that defines a receiver composed of a frequency
resolution 'freq_res', a center frequency 'cent_freq', a bandwidth 'bw', an
amplifier gain 'gain_amps', a receiver temperature 'T_rx' and a frequency
response 'freq_resp'. 

'T_rx' can be a scalar, an array which first dimension
would represent temporal evolution and second dimension the frequnecy evolution
or a 'DimArray' of dimension ':times' and/or ':freqs'.

'freq_resp' must be normalized such that √(1/N * ∑_i∈[1,N] freq_resp[i]^2) ≈ 1,
thus normalizing the receiver frequency response energy.

It is possible to use the 'get_nb_freq_chan' function to get the number of
frequency channels. 'freq_range' can be used with a 'Receiver' as argument to
get the frequency range vector of the receiver.

---
    Receiver(freq_res::T,
             cent_freq::T,
             bw::T,
             gain_amps::T,
             T_rx::Union{T,AbstractArray{T},DimArray{T}}) where T

Yields a 'Receiver' structure with a flat frequency response.

"""
struct Receiver{T<:AbstractFloat}
    freq_res::T # frequency resolution
    cent_freq::T # center frequency
    bw::T # bandwidth
    gain_amps::T # gain of amplifiers
    # receiver temperature (scalar, freq transposed vector or time-freq matrix
    # or DimArray)
    T_rx::Union{T,AbstractArray{T},DimArray{T}}
    freq_resp::AbstractVector{T} # frequency response

    function Receiver(freq_res::T, 
        cent_freq::T,
        bw::T,
        gain_amps::T,
        T_rx::Union{T,AbstractArray{T},DimArray{T}},
        freq_resp::AbstractVector{T}) where T<:AbstractFloat
        
        @assert freq_res > 0
        @assert cent_freq > 0
        @assert bw > 0
        @assert gain_amps > 0
        @assert all(T_rx .>= 0)
        if typeof(T_rx) <: AbstractArray && size(T_rx,2) > 1
            @assert size(T_rx,2) == length(freq_resp) == 
                    div(bw, freq_res) "the second dimension of T_rx must match the \
                                       length of freq_resp and match the number of \
                                       frequency channels defined by bw and freq_res"
        end
        if typeof(T_rx) <: DimArray
            @assert dims(T_rx) in [:times, :freqs] "the dimension of T_rx must be \
                                                     :times and/or :freqs"
            if :freqs in dims(T_rx)
                @assert size(T_rx,:freqs) == length(freq_resp)
            end
        end
        @assert all(sqrt.(sum(freq_resp.^2, 
                             dims=1)./length(freq_resp)) .- 1 .<= 1e-10) "\
                the receiver frequency response must be normalized"

        return new{T}(freq_res, cent_freq, bw, gain_amps, T_rx, freq_resp)
    end
end

function Receiver(freq_res::T,
    cent_freq::T,
    bw::T,
    gain_amps::T, 
    T_rx::Union{T,AbstractArray{T},DimArray{T}}) where {T<:AbstractFloat}

    # flat frequency response
    nb_freq_chan = div(bw, freq_res)
    freq_resp = T.(ones(nb_freq_chan))

    return Receiver(freq_res, cent_freq, bw, gain_amps, T_rx, freq_resp)
end

Base.show(io::IO, R::Receiver) = begin
    print(io, "Receiver:\n")
    print(io, "frequency resolution: $(R.freq_res)\n")
    print(io, "center frequency: $(R.cent_freq)\n")
    print(io, "bandwidth: $(R.bw)\n")
    print(io, "gain amplifiers: $(R.gain_amps)")
end

get_nb_freq_chan(R::Receiver) = Int(div(R.bw, R.freq_res))

freq_range(R::Receiver) = freq_range(R.freq_res, R.cent_freq, get_nb_freq_chan(R), R.bw)



"""
    Instrument(antenna::Antenna{T},
               receiver::Receiver{T},
               coords::Dict{Symbol,Union{T,AbstractVector{T}}} = 
                       Dict(:lat=>0.,:lon=>0.,:alt=>0.)) where T<:AbstractFloat

Yields an 'Instrument' structure that defines an instrument composed of an
antenna and a receiver. The position of the instrument is defined by the
'coords' composed longitude, latitude and altitude coordinates. Note that
instrument's coordinates can evolve with time.

'get_psd_gain_coeff' can be used with an 'Instrument' as argument to get the
instrument gain coefficient (amplifier gain, frequency resolution, Boltzman
constant, impedance and frequency respopnse).

"""
struct Instrument{T<:AbstractFloat,U<:Union{T,AbstractVector{T}}}
    antenna::Antenna{T} # antenna
    #FIXME: Union Transmitter with new struct here?
    receiver::Receiver{T} # receiver of precision T
    coords::Dict{Symbol,U} # coordinates

    function Instrument(antenna::Antenna{T},
        receiver::Receiver{T},
        coords::Dict{Symbol,U} = 
                Dict(:lat => 0.,:lon => 0.,:alt=>0.)) where {T<:AbstractFloat,U<:Union{T,AbstractVector{T}}}

        ant_fmin, ant_fmax = antenna.valid_freqs
        cent_freq = receiver.cent_freq
        bw = receiver.bw
        @assert (ant_fmin <= cent_freq - bw/2) && 
                (cent_freq + bw/2 <= ant_fmax) "the receiver does not cover the antenna \
                                                valid frequency range"
        if size(receiver.T_rx,1) > 1 && size(antenna.T_phy,1) > 1
            if typeof(receiver.T_rx) <: DimArray
                @assert size(receiver.T_rx,:times) == size(antenna.T_phy,1) "the \
                        first dimension (time) of T_rx and T_phy must match"
            else
                @assert size(receiver.T_rx,1) == size(antenna.T_phy,1) "\
                        the first dimension (time) of T_rx and T_phy must match"
            end
        end
        if length(coords[:lat]) > 1 && (size(receiver.T_rx,1) > 1 || 
                                        size(antenna.T_phy,1) > 1)
            max_T_length = max(size(receiver.T_rx,1), size(antenna.T_phy,1))
            @assert length(coords[:lat]) == max_T_length "the length of coords must \
                                                          match the first dimension \
                                                          (time) of T_rx and T_phy, if \
                                                          they are time dependent"
        end
        @assert :lat in keys(coords) "`:lat` must be a key of coords"
        @assert :lon in keys(coords) "`:lon` must be a key of coords"
        @assert :alt in keys(coords) "`:alt` must be a key of coords"

        return new{T,U}(antenna, receiver, coords)
    end
end

Base.show(io::IO, I::Instrument{T}) where {T} = begin
    print(io, "Instrument{$T}:\n")
    print(io, "antenna: $(I.antenna)\n")
    print(io, "receiver: $(I.receiver)\n")
    if length(I.coords[:lat]) > 1
        print(io, "mobile instrument")
    else
        print(io, "instrument at $(I.coords[:lat]), $(I.coords[:lon]), \
                   $(I.coords[:alt])")
    end
end

function get_psd_gain_coeff(I::Instrument)
    return (I.receiver.gain_amps * I.receiver.freq_res * k_boltz * 
            impedance) ./ I.receiver.freq_resp
end



"""
    Trajectory(traj::AbstractDataFrame;
               time_tag::Symbol = :times,
               azimuth_tag::Symbol = :azimuths,
               elevation_tag::Symbol = :altitudes)

Yields a 'Trajectory' structure. The columns of the DataFrame 'traj'  contain the columns
`:times`, `:azimuths` and `:elevations`. It is possible to have vectors of els
and azs for a same time (e.g. for sky mapping). Angles are in degrees.

'get_time_bounds' can be used to get the first and last times of the trajectory.

From a 'Trajectory' structure 'T', t is possible to get the trajectory over two
dates 't0' and 't1' using 'get_traj_between_dates':
'''
get_traj(T, t0, t1; skipmissing=true, view=true)
'''
that yeilds a subset of the full trajectory.

---
    Trajectory(file_path::String;
               kwds...)

Yields a 'Trajectory' structure from the file located at 'file_path'. Only '.arrow'
and '.csv' files are supported for now.

"""
struct Trajectory
    traj::AbstractDataFrame
    
    function Trajectory(traj::AbstractDataFrame;
        time_tag::Symbol = :times,
        azimuth_tag::Symbol = :azimuths,
        elevation_tag::Symbol = :elevations)
        
        # rename columns
        col_names = propertynames(traj)
        for tag in [time_tag, azimuth_tag, elevation_tag]
            @assert tag in col_names "`traj` must contain the column defined by `$tag`"
        end
        select!(traj, time_tag => :times, azimuth_tag => :azimuths, 
                elevation_tag => :elevations)

        @assert length(unique(traj.times)) == length(traj.times) "the trajectory \
                contains repeated times"
        @assert all(length.(traj.azimuths) .== length(traj.azimuths[1])) "\
                the trajectory does not contain the same number of azimuths for \
                different times"
        @assert all(length.(traj.elevations) .== length(traj.elevations[1])) "\
                the trajectory does not contain the same number of elevations for \
                different times"
                
        # convert time stamps into DateTime type
        if typeof(traj.times[1]) != DateTime
            @. traj[!,:times] = Dates.DateTime(traj[!,:times])
        end

        sort!(traj, :times)

        return new(traj)
    end
end

function Trajectory(file_path::String;
    kwds...)

    # load the trajectory as a DataFrame
    #FIXME: read also csv files without dates for real antenna positions
    if occursin(".arrow", file_path)
        traj = DataFrame(Arrow.Table(file_path))
    elseif occursin(".csv", file_path)
        traj = DataFrame(CSV.File(file_path))
    else
        error("the trajectory points are not in Arrow or CSV format")
    end

    return Trajectory(traj; kwds...)
end

get_time_bounds(T::Trajectory) = (T.traj[1,:times], T.traj[end,:times])

function get_traj(T::Trajectory, 
    t0::DateTime, 
    t1::DateTime;
    skipmissing::Bool = true,
    view::Bool = true)
    
    return subset(T.traj, :times => T -> t0 .<= T .<= t1; skipmissing=skipmissing,
                  view=view)
end



"""
Generic type for background models.
"""
abstract type AbstractBackground end



"""
    BackgroundModel(bkg_map::SphereMap)

Yields a 'BackgroundModel' structure.

When applying a 'SphereMap', 'M', to a 'BackgroundModel', 'B'
'''
B(M)
'''
the background model is multiplied by the 'map' of the 'M' for each
elements of its 'alpha_grid' and 'beta_grid'. The result is then an array
sampled over the same grids as 'M'.

#TODO: finish doc

"""
struct BackgroundModel{T} <: AbstractBackground
    bkg_map::SphereMap{T} # background model
end

(B::BackgroundModel)(alpha::T, beta::T) where T = B.bkg_map(alpha, beta)

function (B::BackgroundModel{T})(ant::Antenna) where T

    samp_alpha_grid = ant.gain_pat.alpha_grid
    samp_beta_grid = ant.gain_pat.beta_grid
    pattern_map = ant.gain_pat.map
    result_map = similar(pattern_map)
    fill!(result_map, zero(T))
    @inbounds for b_id in eachindex(samp_beta_grid)
        @simd for a_id in eachindex(samp_alpha_grid)
            a = samp_alpha_grid[a_id]
            b = samp_beta_grid[b_id]
            result_map[a_id, b_id] = pattern_map[a_id,b_id] * B(a, b)
        end
    end

    return result_map
end

function (B1::BackgroundModel)(B2::BackgroundModel)
    
    res = min.(get_angle_resolution(B1.bkg_map), get_angle_resolution(B2.bkg_map))
    grids1 = get_angle_grids(B1.bkg_map)
    grids2 = get_angle_grids(B2.bkg_map)
    alpha_min = min(minimum(grids1[1]), minimum(grids2[1]))
    alpha_max = max(maximum(grids1[1]), maximum(grids2[1]))
    beta_min = min(minimum(grids1[2]), minimum(grids2[2]))
    beta_max = max(maximum(grids1[2]), maximum(grids2[2]))
    T = promote_type(eltype(B1.bkg_map.map), eltype(B2.bkg_map.map))

    samp_alpha_grid = collect(alpha_min:res[1]:alpha_max)
    samp_beta_grid = collect(beta_min:res[2]:beta_max)
    result_map = zeros(T, length(samp_alpha_grid), length(samp_beta_grid))
    @inbounds for b_id in eachindex(samp_beta_grid)
        @simd for a_id in eachindex(samp_alpha_grid)
            a = samp_alpha_grid[a_id]
            b = samp_beta_grid[b_id]
            result_map[a_id, b_id] = B1(a, b) * B2(a, b)
        end
    end

    return BackgroundModel(SphereMap(samp_alpha_grid, samp_beta_grid, result_map))
end



"""
    ground_model(alpha_grid::AbstractVector{T},
                 beta_grid::AbstractVector{T},
                 T_ground::T;
                 mask::AbstractMatrix{Bool} = ones(Bool, length(alpha_grid), 
                                                   length(beta_grid))) where T

Yields a 'BackgroundModel' structure representing a ground model with
constant temperature 'T_ground'. The 'mask' keyword argument can be used to
mask out more than the horizon (e.g. for mountains). The mask must be a
Boolean matrix with size equal to (length(beta_grid), length(alpha_grid)).

"""
function ground_model(alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T},
    T_ground::T;
    mask::AbstractMatrix{Bool} = ones(Bool, length(alpha_grid), 
                                      length(beta_grid))) where T


    @assert size(mask) == (length(beta_grid), length(alpha_grid))

    bellow_horizon = beta_grid .> (π/2)
    mask .= mask .& .!bellow_horizon
    
    return BackgroundModel(SphereMap(alpha_grid, beta_grid, mask .* T_ground))
end



"""
    atmosphere_model(alpha_grid::AbstractVector{T},
                     beta_grid::AbstractVector{T},
                     T_bkg::T,
                     T_eff::T,
                     zenith_opacity::T;) where T

Yields a 'BackgroundModel' structure representing an atmosphere model with
constant temperature 'T_eff' and opacity 'zenith_opacity'. 'T_bkg' is the
background temperature behind the atmosphere (e.g. CMB).

"""
function atmosphere_model(alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T},
    T_bkg::T,
    T_eff::T,
    zenith_opacity::T;) where T
    
    atm_els = atmos_opacity_impact(T_bkg .- T_eff, zenith_opacity, beta_grid) .+ T_eff
    atm_map = reduce(hcat, [atm_els for _ in 1:length(alpha_grid)])

    return BackgroundModel(SphereMap(alpha_grid, beta_grid, atm_map))
end



"""
    Observation(antenna_traj::Trajectory,
        instrument::Instrument{T},
        result::Union{AbstractArray{T},DimArray{T}}) where T

Yields an 'Observation' structure that stores the trajectory of the antenna
during the observation, the instrument used and the observation results.

If 'result' is an 'AbstractArray' it must have its first dimension equal to the
number of time stamps in 'antenna_traj'. The second dimension must be equal to
the number of frequency channels in 'instrument.receiver'. The third dimension
must be equal to the number of azimuths per time stamp in 'antenna_traj'. This
means it is possible to create a 'Trajectory' with multiple azimuths per time stamp
(e.g. for sky mapping) and store the observation results for all azimuths.

If 'result' is a 'DimArray' it must have its dimensions named :times, :freqs and
:azimuths. The size of the :times dimension must be equal to the number of time
stamps in 'antenna_traj'. The size of the :freqs dimension must be equal to the number of
frequency channels in 'instrument.receiver'. The size of the :azimuths dimension must be
equal to the number of azimuths per time stamp in 'antenna_traj'.

---
    Observation(start_date::DateTime,
        stop_date::DateTime,
        trajectory::Trajectory,
        instrument::Instrument{T};
        filt_funcs::NTuple{N,Pair} = ()) where {N,T}

Yields an 'Observation' structure that stores the trajectory of the antenna
during the observation, the instrument used and the observation results. The
trajectory is filtered to only keep the points between 'start_date' and
'stop_date'. Additional filters can be applied using the 'filt_funcs' keyword
argument. Each element of 'filt_funcs' must be a Pair where the key is the
column name to filter and the value is a function that returns a boolean array
to filter the column.

"""#FIXME: make result elements be vectors if vectors in azimuth trajectory 
# => get rid of azimuths dim + avoid storing a lot of NaNs in DimArray for
# (time,azimuth) combinations not existing in trajectory
struct Observation{T<:AbstractFloat}
    antenna_traj::Trajectory # antenna trajectory during observation
    instrument::Instrument{T} # instrument used for observation
    result::Union{AbstractMatrix,DimArray} # store the results of observation

    function Observation(antenna_traj::Trajectory,
        instrument::Instrument{T},
        result::Union{AbstractMatrix,DimArray}) where T
        
        if typeof(result) <: AbstractMatrix{T}
            @assert size(result, 1) == length(antenna_traj.traj.times) "\
                    the first dimension of result must match the number of time \
                    stamps of the antenna trajectory"
            @assert size(result, 2) == get_nb_freq_chan(instrument.receiver) "\
                    the second dimension of result must match the number of frequency \
                    channels of the instrument"
            @assert size(result, 3) == length(antenna_traj.traj.azimuths[1]) "\
                    the third dimension of result must match the number of azimuths \
                    per time stamp of the antenna trajectory"
            if size(instrument.antenna.T_phy, 1) > 1
                @assert size(instrument.antenna.T_phy, 1) == size(result, 1) "\
                        the first dimension of result must match the number of time \
                        stamps of the instrument physical temperature 'T_phy'"
            end
            if size(instrument.receiver.T_rx, 1) > 1
                @assert size(instrument.receiver.T_rx, 1) == size(result, 1) "\
                        the first dimension of result must match the number of time \
                        stamps of the instrument receiver temperature 'T_rx'"
            end
        elseif typeof(result) <: DimArray{T}
            @assert name.(dims(result)) == [:times, :freqs, :azimuths] "\
                    the dimensions of result must be named :times, :freqs and :azimuths"
            @assert size(result, :times) == length(antenna_traj.traj.times) "\
                    the :times dimension of result must match the number of time \
                    stamps of the antenna trajectory"
            @assert size(result, :freqs) == get_nb_freq_chan(instrument.receiver) "\
                    the :freqs dimension of result must match the number of frequency \
                    channels of the instrument"
            @assert size(result, :azimuths) == length(antenna_traj.traj.azimuths[1]) "\
                    the :azimuths dimension of result must match the number of azimuths \
                    per time stamp of the antenna trajectory"
            if size(instrument.antenna.T_phy, 1) > 1
                @assert size(instrument.antenna.T_phy, :times) == size(result, :times) "\
                        the :times dimension of result must match the number of time \
                        stamps of the instrument physical temperature 'T_phy'"
            end
            if size(instrument.receiver.T_rx, 1) > 1
                @assert size(instrument.receiver.T_rx, :times) == size(result, :times) "\
                        the :times dimension of result must match the number of time \
                        stamps of the instrument receiver temperature 'T_rx'"
            end
        end

        return new{T}(antenna_traj, instrument, result)
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
    time_stamps = get_time_stamps(pts)
    freq_bins = freq_range(instrument.receiver)
    azimuths = get_azimuths(pts)
    result_array = fill!(zeros(T, length(time_stamps), length(freq_bins), 
                         length(azimuths)), NaN)
    result = DimArray(result_array, (Dim{:times}(time_stamps), Dim{:freqs}(freq_bins), 
                                     Dim{:azimuths}(azimuths)))

    return Observation(pts, instrument, result)
end

function get_antenna_temperature(A::Antenna{T}, 
    T_b::SphereMap{T}) where T

    # FIXME: do the real convolution (project background in antenna frame)
    # return get_boresight_gain(A)[1] .* T_b ./ (4*pi)#FIXME: WRONG!!!!!!!!!!

    
end



"""
    instrument_psd_stat(i::Instrument{T},
                        T_b::Union{T,<:AbstractArray{T},DimArray{T}},
                        integration_samp::Real = 1) where T

Yields the power spectral density and its variance for the given instrument and
sky brightness temperature 'T_b'.

"""
function instrument_psd_stat(I::Instrument{T},
    T_b::Union{T,<:AbstractArray{T},DimArray{T}},
    integration_samp::Real = 1) where T

    # get receiver and antenna parameters
    rec = I.receiver
    ant = I.antenna

    # instrument gain coefficient
    gain = get_psd_gain_coeff(I)

    # antenna temperature
    T_a = get_antenna_temperature(ant, T_b)
    
    # instrument noise temperature
    T_n = rec.T_rx .+ get_antenna_radiation_loss(ant)

    # calculate power spectral density
    return intrument_psd_stat(gain, T_a, T_n, integration_samp)
end



"""
    Satellite(name::String,
              trajectory::Trajectory,
              instrument::Instrument{T}) where T

Yields a 'Satellite' structure that defines a satellite with its name, trajectory
and instrument. Assumes the frame of satellite antenna is oriented
North-East-Nadir. The antenna pointing can be any direction from Nadir, defined
in the 'gain_pat'. 

"""
struct Satellite{T<:AbstractFloat}
    name::String # name of satellite
    trajectory::Trajectory # trajectory of satellite
    instrument::Instrument{T} # instrument of satellite
end



"""

Forms a list of Satellite structs from a DataFrame containing satellite
information.
"""
function form_satellites_list(sats_info::AbstractDataFrame) #FIXME: fill up

end



"""
"""
struct Constellation{T<:AbstractFloat}
    constellation_name::String # name of constellation
    sats::AbstractVector{Satellite{T}} # satellites in constellation
    lnk_bdgt_mdl::Function # link budget model function

    function Constellation(constellation_name::String,
        sats::AbstractVector{Satellite{T}},
        lnk_bdgt_mdl::Function) where T

        # check lnk_bdgt_mdl signature is correct
        @assert hasmethod(lnk_bdgt_mdl, (T, T, Instrument{T}, T, T, T, Instrument{T}, T))

        return new{T}(constellation_name, sats, lnk_bdgt_mdl)
    end
end

function Constellation(start_date::DateTime, #FIXME: fill up
    stop_date::DateTime,
    constellation_name::String,
    sats::AbstractDataFrame,
    sat_inst::Instrument{T},
    lnk_bdgt_mdl::Function = sat_link_budget;
    filt_funcs::NTuple{N,Pair} = ()) where {N,T}

    # sats = Satellite{T}[]
    # for sat_info in sat_infos
    #     # create trajectory
    #     traj = Trajectory(sat_info.traj_file_path; traj_kwds...)

    #     # create instrument
    #     inst = Instrument{T}(sat_info.antenna, sat_info.receiver, sat_info.coords)

    #     # create satellite
    #     push!(sats, Satellite{T}(sat_info.name, traj, inst))
    # end

    return Constellation(constellation_name, sats, lnk_bdgt_mdl)
end

function Constellation(file_path::String,  #FIXME: fill up
    observation::Observation,
    sat_tmt::Instrument{T},
    lnk_bdgt_mdl::Function = sat_link_budget;
    name_tag::Symbol = :sat,
    time_tag::Symbol = :time_stamps,
    elevation_tag::Symbol = :altitudes,
    azimuth_tag::Symbol = :azimuths,
    distance_tag::Symbol = :distances,
    filt_funcs::NTuple{N,Pair} = ()) where {N,T}

    # # load the trajectory as a DataFrame
    # sats = DataFrame(Arrow.columntable(Arrow.Table(file_path)))

    # # rename columns
    # rename!(sats, time_tag => :times)
    # rename!(sats, name_tag => :sat)
    # rename!(sats, azimuth_tag => :azimuths)
    # rename!(sats, elevation_tag => :elevations)
    # rename!(sats, distance_tag => :distances)
    
    # @. sats[!,:times] = Dates.DateTime(sats[!,:times])

    # sort!(sats, :times)

    # return Constellation(sats, observation, sat_tmt, lnk_bdgt_mdl; 
    #                      filt_funcs=filt_funcs)
end

function get_sat_traj(C::Constellation, #FIXME: fill up
    s::String)

    # return subset(C.sats, :sat => n -> n .== s; view=true)
end

function get_sats_names_at_time(C::Constellation, #FIXME: fill up
    t::DateTime)
    
    # return subset(C.sats, :times => ts -> ts .== t; view=true)[!,:sat]
end

