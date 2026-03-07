module RadioMdl

export AbstractBkgModel,
       Antenna,
       adc_noise_temperature,
       add_coords,
       antenna_mdl_cst,
       antenna_mdl_ITU_RA_1631,
       antenna_mdl_ITU_SA_509_3,
       atmos_opacity_impact,
       atmosphere_model,
       Constellation,
       classic_gain_link_budget,
       estim_casA_flux,
       estim_hpbws,
       estim_virgoA_flux,
       flux_to_temperature,
       form_satellites_list,
       freq_range,
       freq_to_wave,
       friis_noise_temp,
       gain_to_effective_aperture,
       galactic_model,
       get_angle_grids,
       get_antenna_radiation_loss,
       get_antenna_temperature,
       get_beam_solid_angle,
       get_boresight_gain,
       get_coords,
       get_directivity_value,
       get_effective_aperture,
       get_gain_value,
       get_geometric_effective_aperture,
       get_nb_freq_chan,
       get_psd_gain_coeff,
       get_sat,
       get_sat_EIRP_density,
       get_sat_traj,
       get_sats_name,
       get_sats_names_at_time,
       get_time_bounds,
       ground_model,
       Instrument,
       instrument_psd_stat,
       k_boltz,
       MovingExtendSrcTemp,
       model_observ_psd!,
       Observation,
       offset_angle_trajectory!,
       PointLikeSrcFlux,
       power_to_temperature,
       Receiver,
       radiated_power_to_gain!,
       read_VGOS_antenna_traj,
       Satellite,
       SkyMdl,
       SphereCoord,
       SphereMap,
       TiFreqArray,
       Trajectory,
       temperature_to_flux,
       temperature_to_power,
       wave_to_freq

# export adc_noise_temperature,
    #    Antenna,
    #    antenna_mdl_cst,
    #    antenna_mdl_ITU,
    #    Constellation,
    #    estim_casA_flux,
    #    estim_hpbws,
    #    estim_temp,
    #    estim_virgoA_flux,
    #    fetch_satellites_info,
    #    flux_to_temperature,
    #    free_space_loss,
    #    freq_range,
    #    freq_to_wave,
    #    friis_noise_temp,
    #    gain_to_effective_aperture,
    #    get_antenna_radiation_loss,
    #    get_antenna_temperature,
    #    get_boresight_gain,
    #    get_center_freq_chans,
    #    get_def_angles,
    #    get_directivity_value,
    #    get_gain_value,
    #    get_geometric_effective_aperture,
    #    get_result,
    #    get_sats_name,
    #    get_sats_names_at_time,
    #    get_sat_traj,
    #    get_slice_gain,
    #    get_time_bounds,
    #    get_time_stamps,
    #    get_traj,
    #    Instrument,
    #    instrument_psd_stat,
    #    k_boltz,
    #    map_sphere,
    #    model_observed_temp!,
    #    Observation,
    #    power_to_temperature,
    #    Receiver,
    #    sat_link_budget,
    #    simple_link_budget,
    #    speed_c,
    #    Trajectory,
    #    temperature_to_flux,
    #    temperature_to_power,
    #    wave_to_freq

using Arrow
using CSV
using DataFrames
using Dates
using DelimitedFiles
using DimensionalData
using HTTP
using Interpolations
using InterpolationKernels
using LinearInterpolators
using SpecialFunctions
using Trapz
using .Threads

using LinearAlgebra: dot, BLAS
BLAS.set_num_threads(1)

""" Boltzman's constant in J/K"""
const k_boltz = 1.380649e-23

""" speed of light in m/s """
const speed_c = 299792458

""" impedance match in Ohm """
const impedance = 50


include("io.jl")
include("coord_frames.jl")
using .CoordFrames
include("antenna_pattern.jl")
include("rx_mdl.jl")
include("astro_mdl.jl")
include("types.jl")
include("sat_mdl.jl")
include("obs_mdl.jl")

end # module RadioMdl
