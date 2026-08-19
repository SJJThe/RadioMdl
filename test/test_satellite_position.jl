
using .Threads
using DataFrames
using Dates
using PyPlot
const plt = PyPlot
using Statistics

using Revise
using RadioMdl

coords = Dict(:lat => 42.6129479883915, :lon => -71.49379366344017, 
              :alt => 86.7689687917009)

dateformat = "yyyy-mm-dd\\THH:MM:SS.sss"
start_obs = DateTime("2026-08-09T10:20:13.000", dateformat)
stop_obs = DateTime("2026-08-09T10:22:13.000", dateformat)

name_filters = ["STARLINK"]
avoid_names = ["[DTC]"]
sats_catalog = fetch_satellites_info(; name_filters=name_filters, 
                                     avoid_names=avoid_names, save=false, verb=true)

min_elevation_filter = 5.

sats_pos_julia = compute_sats_traj(sats_catalog, start_obs, stop_obs, coords, Second(1);
                                   save=false, el_min=min_elevation_filter)

sats_pos_python = compute_sats_traj_py(sats_catalog, start_obs, stop_obs, coords, 
                                       Second(1); save=false,
                                       el_min=min_elevation_filter)

sats_julia = unique(sats_pos_julia[!,:sat])
sats_python = unique(sats_pos_python[!,:sat])

common_sats = [s for s in sats_python if s in sats_julia]

mean_diffs = zeros(length(common_sats))
std_diffs = zeros(length(common_sats))
@threads for i in eachindex(common_sats)
    sat = common_sats[i]
    sat_pos_jl = subset(sats_pos_julia, :sat => s -> s .== sat; view=true)
    sat_pos_py = subset(sats_pos_python, :sat => s -> s .== sat; view=true)

    ang_sep = Float64[]
    for (ind_t_jl, t_jl) in enumerate(sat_pos_jl[!,:times])
        t_diff, ind_t_py = findmin(abs.(t_jl .- sat_pos_py[!,:times]))
        if t_diff .<= Millisecond(110)
            sat_vec_jl = [cosd(sat_pos_jl[ind_t_jl,:elevations]) * 
                          cosd(sat_pos_jl[ind_t_jl,:azimuths]),
                          cosd(sat_pos_jl[ind_t_jl,:elevations]) * 
                          sind(sat_pos_jl[ind_t_jl,:azimuths]),
                          sind(sat_pos_jl[ind_t_jl,:elevations])]
            sat_vec_py = [cosd(sat_pos_py[ind_t_py,:elevations]) * 
                          cosd(sat_pos_py[ind_t_py,:azimuths]),
                          cosd(sat_pos_py[ind_t_py,:elevations]) * 
                          sind(sat_pos_py[ind_t_py,:azimuths]),
                          sind(sat_pos_py[ind_t_py,:elevations])]

            push!(ang_sep, acosd(sum(sat_vec_jl .* sat_vec_py)))
        end
    end

    mean_diffs[i] = mean(ang_sep)
    std_diffs[i] = std(ang_sep)
end

fig, axs = plt.subplots()
axs.errorbar(1:length(common_sats), mean_diffs, yerr=std_diffs, fmt="o")
axs.set_xlabel("Common satellites")
axs.set_ylabel("Angular Separation [deg]")
axs.grid()
fig.tight_layout()

