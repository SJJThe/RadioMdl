module RadioMdl

export adc_noise_temperature,
       Antenna,
       antenna_mdl_cst,
       antenna_mdl_ITU,
       Constellation,
       estim_casA_flux,
       estim_hpbws,
       estim_max_effective_aperture,
       estim_temp,
       estim_virgoA_flux,
       flux_to_temperature,
       free_space_loss,
       freq_to_wave,
       friis_noise_temp,
       gain_to_effective_aperture,
       get_antenna_radiation_loss,
       get_antenna_temperature,
       get_boresight_gain,
       get_center_freq_chans,
       get_def_angles,
       get_directivity_value,
       get_gain_value,
       get_geometric_effective_aperture,
       get_result,
       get_sats_name,
       get_sats_names_at_time,
       get_sat_traj,
       get_slice_gain,
       get_time_bounds,
       get_time_stamps,
       get_traj,
       Instrument,
       instrument_psd_stat,
       k_boltz,
       map_sphere,
       model_observed_temp!,
       Observation,
       power_to_temperature,
       Receiver,
       sat_link_budget,
       simple_link_budget,
       speed_c,
       Trajectory,
       temperature_to_flux,
       temperature_to_power,
       wave_to_freq

using Arrow
using DataFrames
using Dates
using DelimitedFiles
using DimensionalData
using Interpolations
using Trapz
using YAXArrays
using .Threads

""" Boltzman's constant in J/K"""
const k_boltz = 1.380649e-23

""" degree to radian conversion factor """
const rad = π/180

""" speed of light in m/s """
const speed_c = 299792458

""" impedance match in Ohm """
const impedance = 50

include("io.jl")
include("astro_mdl.jl")
include("antenna_pattern.jl")
include("coord_frames.jl")
include("types.jl")
include("sat_mdl.jl")
include("obs_mdl.jl")

end # module RadioMdl
