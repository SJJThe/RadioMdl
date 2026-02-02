
"""
    power_pattern_from_cut_file(file_path::String;
                                free_sp_imp::Real = 377,
                                verb::Bool = false)
                                
Yields the radiated power pattern, in W, of an antenna, times the radiation
efficiency already included in the `.cut` file containing co- and
cross-polarization E-field. 

Headers in the file are below a line starting with `Field`. It is composed of
the starting value of the declination angle, the step and number of samples of
the declination angle and the value of the azimuthal angle.

"""
function power_pattern_from_cut_file(file_path::String;
    free_sp_imp::Real = 377,
    verb::Bool = false)
    
    @assert occursin(".cut", file_path) "the power pattern file must be a .cut file"

    # parse file
    Es = readdlm(file_path)
    pattern = DataFrame(polar=Float64[], az=Float64[], power=Float64[])
    k = 1
    dec_step = 0.
    while k <= size(Es,1)
        header_line = k + findfirst(x -> x == "Field", Es[k:end,:])[1]
        header = Es[header_line,:]
        verb && println(header)
        dec_start = header[1]
        dec_step = header[2]
        nb_dec = header[3]
        for t in 1:nb_dec
            polar = dec_start + (t-1)*dec_step
            θ = header[4]
            # power pattern, given in dBW, is the sum of the magnitude (squared
            # modulus) of the co- and cross-polarization complex electric field,
            # devided by twice the free-space impedance
            @inbounds u = sum(Es[header_line+t,1:4].^2)/(2*free_sp_imp)
            push!(pattern, [polar, θ, u])
        end
        k = header_line+nb_dec+1
    end
    decimal_places = max(0, -floor(Int, log10(abs(dec_step - round(dec_step)))))
    pattern[!,:polar] .= round.(pattern[!,:polar]; digits=decimal_places)
    
    # !!!!!!!!!!!!!!! THIS IS ONLY THE CASE WITH DANIEL'S FORMAT !!!!!!!!!!!!!!!
    
    @warn "This function assumes TICRA generated files"

    # check that polar ∈ [-180,180] and az ∈ [0, 180[
    subset!(pattern, :polar => p -> -180. .<= p .<= 180., 
            :az => a -> 0. .<= a .< 180.)

    # at this point, when the telescope is pointed at the horizon, az = 0 gives
    # an horizontal slice, with polar > 0 oriented towards co-azimuth angles..

    # change interval so that az ∈ [0,360[ and polar ∈ [0,180]
    pattern[pattern.polar .<= 0.,:az] .+= 180.
    pattern[pattern.polar .< 0.,:polar] .*= -1.
    append!(pattern, [(;polar = zero(eltype(pattern.polar)), az = i, 
                       power = pattern[pattern.polar .== 0.,:power][1])
                      for i in pattern[pattern.polar .== maximum(pattern.polar) .&& 
                                       pattern.az .< 180.,:az]])
    
    # move the origin of az so that, when telescope points at the horizon, the
    # first slice (for the new az = 0) is vertical with polar > 0 oriented towards
    # the ground.
    pattern[:,:az] = mod.(pattern[!,:az] .- 90., 360.)
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    sort!(pattern, [:az, :polar])
    
    return pattern
end



""" Celestrak active satellites url """
const celestrak_url = "https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=csv"



"""
"""
function fetch_satellites_info(;
    csv_path::String = celestrak_url,
    name_filters::Union{Nothing,String,AbstractVector{<:String}} = nothing,
    avoid_names::Union{Nothing,String,AbstractVector{<:String}} = nothing,
    verb::Bool = false,
    save::Bool = false)

    if "https" in csv_path
        csv_sats_info = String(HTTP.get(csv_path).body)
        sats_catalog = CSV.read(IOBuffer(csv_sats_info), DataFrame)
        # save the fetched csv to reduce access to website
        if save
            CSV.write("sats_catalog_$(replace(string(now()), ":" => "_")).csv", 
                      sats_catalog)
        end
    else
        sats_catalog = CSV.read(csv_path, DataFrame)
    end

    if !isnothing(name_filters)
        for filt in name_filters
            filter!(row -> occursin(filt, row.OBJECT_NAME), sats_catalog)
        end
    end
    if !isnothing(avoid_names)
        for filt in avoid_names
            filter!(row -> !(occursin(filt, row.OBJECT_NAME)), sats_catalog)
        end
    end

    if verb
        println("Found $(size(sats_catalog,1)) satellites with the given filters")
    end

    return sats_catalog
end

