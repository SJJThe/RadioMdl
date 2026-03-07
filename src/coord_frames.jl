
"""
    module CoordFrames

Defines the generic coordinate frame '(X,Y,Z)' used in this package as:

         Z
         |    p (α,β)
         | β /'
         |--/ '
         | /  '
         |/___'______Y
        / `   '
       /---`  '
      /  α  ` '
     /       `'
    X

with 'p' the point of spherical coordinates (α,β) in degrees, or alternatively
of cartesian coordinates (x,y,z).

"""
module CoordFrames

using DataFrames
using Interpolations

export add_coords,
       cart_to_sphe_coord,
       get_angle_grids,
       get_angle_resolution,
       interpolate_sphere_map,
       map_sphere_coords,
       pass_frame_to_frame,
       rot_mat,
       SphereCoord,
       SphereMap,
       spher_to_cart_coord



"""
    SphereCoord(alpha::T,
                beta::T, 
                r::T = one(T)) where T<:AbstractFloat

Defines the spherical coordinates '(α,β,r)' of a point in the '(X,Y,Z)' frame,
with 'α' ∈ [0, 360], 'β' ∈ [0, 180].

See ['CoordFrames'](@ref)

"""
struct SphereCoord{T<:AbstractFloat}
    alpha::T
    beta::T
    r::T

    function SphereCoord(alpha::T,
        beta::T,
        r::T = one(T)) where T<:AbstractFloat
        
        new{T}(mod(alpha, T(360)), mod(beta, T(180)), r)
    end
end
#TODO: add kwd for if elevation and azimuth are given


"""
    add_coords(a::SphereCoord{T},
               b::SphereCoord{T};
               ranges::Bool = false) where T

Returns the sum of two spherical coordinates '(α,β,r)' in the '(X,Y,Z)' frame.

"""
function add_coords(a::SphereCoord{T},
    b::SphereCoord{T};
    ranges::Bool = false) where T

    return SphereCoord(a.alpha + b.alpha, a.beta + b.beta, ranges ? a.r + b.r : one(T))
end



"""
    SphereMap(alpha_grid::AbstractVector{T},
              beta_grid::AbstractVector{T}, 
              spheremap::AbstractMatrix{T}, 
              interp_map::I) where {T,I<:Interpolations.GriddedInterpolation}

Defines a spherical mapping in a '(X,Y,Z)' coordinate frame, where 'alpha_grid'
and 'beta_grid' are sampled grids respectively defining the angles between 'X'
and 'Y' and 'Z' and 'X', in degrees ( see ['CoordFrames'](@ref) ). The
'spheremap' must be a matrix of size (length(alpha_grid), length(beta_grid)).
The 'interp_map' is an interpolator created from the 'spheremap'( see
['Interpolations.jl'](@ref) ).

The resulting 'SphereMap', 'SM', structure can be used to find the interpolated value
of the spheremap at any point on the sphere:
```
SM(180, 90) # returns the interpolated value of the spheremap at α = 180 deg, β = 90 deg

# or
c = SphereCoord(180, 90)
SM(c)
```

It is possible to add a 'Real' value, 'a' to the 'SphereMap', 'SM':
'''
SM + a
'''

It is possible to add to 'SphereMaps', 'S1' and 'S2':
'''
S1 + S2
'''
which will look for the common sampling grid between the two maps and yield a
new 'SphereMap' with the sum of the two maps interpolated on this common grid
(see ['get_common_grid'](@ref)). 

It is possible to apply the 'SphereMap', 'S1', to another
'SphereMap', 'S2'
'''
S1(S2)
'''
which will yield a new 'SphereMap'. The two maps are interpolated to new
'alpha_grid' and 'beta_grid' based on the smallest resolution in the 2
dimensions and their min and max values (see ['get_common_grid'](@ref)).

It is possible to use the 'get_angle_grids' method to retrieve the sampling
grids of the 'SphereMap'. 
'get_angle_resolution' yields the angular resolution of the 'SphereMap' in both
directions.

---
    SphereMap(alpha_grid::AbstractVector{T},
              beta_grid::AbstractVector{T},
              spheremap::AbstractMatrix{T}) where T<:AbstractFloat

Yields a SphereMap struct with an interpolator created from the given 'spheremap',
using the ['interpolate_sphere_map'](@ref) function.

---
    SphereMap(spheremap::AbstractDataFrame;
              kwds...)

Yields a 'SphereMap' struct created from a DataFrame 'spheremap' containing the
columns of angles and map values. The columns can be specified using the
keywords 'alpha_col', 'beta_col' and 'map_col' (default values are
':az', ':polar' and ':power' respectively)( see ['map_sphere_coords'](@ref] ).
The angles must be in degrees.

"""
struct SphereMap{T,V<:AbstractVector{T},M<:AbstractMatrix{T},
                 I<:Interpolations.GriddedInterpolation}
    alpha_grid::V
    beta_grid::V
    spheremap::M
    interp_map::I

    function SphereMap(alpha_grid::V,
        beta_grid::V,
        spheremap::M,
        interp_map::I) where {T,V<:AbstractVector{T},M<:AbstractMatrix{T},
                              I<:Interpolations.GriddedInterpolation}
        
        @assert all(0 .<= alpha_grid .<= 360.)
        @assert all(0 .<= beta_grid .<= 180.)
        @assert length(alpha_grid) > 1
        @assert length(beta_grid) > 1
        @assert length(alpha_grid) == size(spheremap, 1) == size(interp_map, 1) - 1
        @assert length(beta_grid) == size(spheremap, 2) == size(interp_map, 2)
        @assert eltype(interp_map) == eltype(spheremap)

        new{T,V,M,I}(alpha_grid, beta_grid, spheremap, interp_map)
    end
end

function SphereMap(alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T},
    spheremap::AbstractMatrix{T}) where T

    # create the interpolator
    interp_map = interpolate_sphere_map(spheremap, alpha_grid, beta_grid)

    return SphereMap(alpha_grid, beta_grid, spheremap, interp_map)
end

function SphereMap(spheremap::AbstractDataFrame;
    kwds...) 

    map_values, alpha_grid, beta_grid = map_sphere_coords(spheremap; kwds...)

    return SphereMap(alpha_grid, beta_grid, map_values)
end

get_angle_grids(SM::SphereMap) = (SM.alpha_grid, SM.beta_grid)

get_angle_resolution(SM::SphereMap) = (SM.alpha_grid[2] - SM.alpha_grid[1],
                                       SM.beta_grid[2] - SM.beta_grid[1])



"""
    get_common_grid(S1::SphereMap, S2::SphereMap)

Yields the common sampling grid between the two SphereMaps 'S1' and 'S2'. The
common grid is defined as the grid with the smallest resolution in the 2
dimensions and with min and max values corresponding to the min and max values
between the two maps.

"""
function get_common_grid(S1::SphereMap, S2::SphereMap)

    # define smallest resolution in the 2 dimensions
    res = min.(get_angle_resolution(S1), get_angle_resolution(S2))

    # define min and max values of joint SphereMap alphas and betas
    grids1 = get_angle_grids(S1)
    grids2 = get_angle_grids(S2)
    alpha_min = min(minimum(grids1[1]), minimum(grids2[1]))
    alpha_max = max(maximum(grids1[1]), maximum(grids2[1]))
    beta_min = min(minimum(grids1[2]), minimum(grids2[2]))
    beta_max = max(maximum(grids1[2]), maximum(grids2[2]))

    # create new sampling grid 
    samp_alpha_grid = collect(alpha_min:res[1]:alpha_max)
    samp_beta_grid = collect(beta_min:res[2]:beta_max)

    return samp_alpha_grid, samp_beta_grid
end

Base.:+(SM::SphereMap{T}, a::Real) where {T} = begin    
    SphereMap(SM.alpha_grid, SM.beta_grid, SM.spheremap .+ T(a))
end

Base.:+(S1::SphereMap, S2::SphereMap) = begin

    # get common grid
    samp_alpha_grid, samp_beta_grid = get_common_grid(S1, S2)

    # define type of resulting map and initialize
    T = promote_type(eltype(S1.spheremap), eltype(S2.spheremap))
    result_map = zeros(T, length(samp_alpha_grid), length(samp_beta_grid))
    @inbounds for b_id in eachindex(samp_beta_grid)
        @simd for a_id in eachindex(samp_alpha_grid)
            a = samp_alpha_grid[a_id]
            b = samp_beta_grid[b_id]

            # interpolate over both BackgroundModels
            result_map[a_id, b_id] = S1(a, b) + S2(a, b)
        end
    end

    return SphereMap(samp_alpha_grid, samp_beta_grid, result_map)
end

Base.show(io::IO, SM::SphereMap{T}) where {T} = begin 
    print(io, "SphereMap{$T}:\n")
    print(io, "sample grid size: $(size(SM.spheremap))\n")
    print(io, "angles resolutions (degrees): $(get_angle_resolution(SM))")
end

(SM::SphereMap)(alpha::Real, beta::Real) = SM.interp_map(alpha, beta)

(SM::SphereMap)(C::SphereCoord) = SM.interp_map(C.alpha, C.beta)

function (S1::SphereMap)(S2::SphereMap)
    
    # get common grid
    samp_alpha_grid, samp_beta_grid = get_common_grid(S1, S2)
    
    # define type of resulting map and initialize
    T = promote_type(eltype(S1.spheremap), eltype(S2.spheremap))
    result_map = zeros(T, length(samp_alpha_grid), length(samp_beta_grid))
    @inbounds for b_id in eachindex(samp_beta_grid)
        @simd for a_id in eachindex(samp_alpha_grid)
            a = samp_alpha_grid[a_id]
            b = samp_beta_grid[b_id]

            # interpolate over both BackgroundModels
            result_map[a_id, b_id] = S1(a, b) * S2(a, b)
        end
    end

    return SphereMap(samp_alpha_grid, samp_beta_grid, result_map)
end



"""
    pass_frame_to_frame(SM::SphereMap{T},
                        Gamma::Real,
                        Psi::Real,
                        new_alpha_grid::AbstractVector{T} = SM.alpha_grid, 
                        new_beta_grid::AbstractVector{T} = SM.beta_grid;
                        kwds...) where T

Yields a new map where the 'SM.spheremap' has been converted to a new frame
defined by the rotation angles 'Gamma' and 'Psi' (in degrees) in the original
frame.

It is possible to give new angle grids to resample the resulting map.

'new_alpha_grid' and 'new_beta_grid' can be 'AbstractVector's or
'AbstractRange's.

---
    pass_frame_to_frame(SM::SphereMap{T},
                        coordsphere::SphereCoord{T},
                        args...;
                        kwds...) where T

Uses 'pass_frame_to_frame' on a 'SphereMap' given a 'SphereCoord' to define the
orientation of the new frame.

"""
function pass_frame_to_frame(SM::SphereMap{T},
    Gamma::Real,#FIXME:change name varaibles
    Psi::Real,#FIXME:change name varaibles
    new_alpha_grid::AbstractVector{T} = SM.alpha_grid, 
    new_beta_grid::AbstractVector{T} = SM.beta_grid;
    kwds...) where T
    
    itp = SM.interp_map
    new_map = Matrix{T}(undef, length(new_alpha_grid), length(new_beta_grid))
    @inbounds for b in eachindex(new_beta_grid)#TODO: maybe manually vectorized
        @simd for a in eachindex(new_alpha_grid)
            # For each point in NEW frame, find where it came from in OLD frame
            alpha_orig, beta_orig = pass_frame_to_frame(new_alpha_grid[a], 
                                                        new_beta_grid[b], Gamma, Psi;
                                                        inverse = true,
                                                        kwds...)

            # sample original pattern
            new_map[a,b] = itp(alpha_orig, beta_orig)
        end
    end

    return new_map
end

function  pass_frame_to_frame(SM::SphereMap{T},
    Gamma::Real,
    Psi::Real,
    new_alpha_grid::AbstractRange{T}, 
    new_beta_grid::AbstractRange{T};
    kwds...) where T
    
    return pass_frame_to_frame(SM, Gamma, Psi, collect(new_alpha_grid),
                               collect(new_beta_grid); kwds...)
end

function pass_frame_to_frame(SM::SphereMap{T},
    spherecoord::SphereCoord{T},
    args...;
    kwds...) where T

    return pass_frame_to_frame(SM, spherecoord.alpha, spherecoord.beta, args...; kwds...)
end



"""
    map_sphere_coords(spheremap::AbstractDataFrame;
                      alpha_col::Symbol = :az,
                      beta_col::Symbol = :polar,
                      map_col::Symbol = :gains)

Yields the 2D-matrix of 'spheremap', formed from the column 'map_col', with rows and
columns corresponding to the angles in 'alpha_col' and 'beta_col', respectively. 

"""
function map_sphere_coords(spheremap::AbstractDataFrame;
    alpha_col::Symbol = :az,
    beta_col::Symbol = :polar,
    map_col::Symbol = :gains)
    
    alpha_grid = unique(spheremap[!, alpha_col])
    beta_grid = unique(spheremap[!, beta_col])

    # form 2D matrix for interpolation argument
    # add first column as last column to loop horizontal coordinates
    map_vec = spheremap[!, map_col]
    map_values = zeros(eltype(map_vec), length(alpha_grid), length(beta_grid))
    for (i, a) in enumerate(alpha_grid)
        beta_col_a = subset(spheremap, alpha_col => az -> az .== a; view=true)[!,map_col]
        map_values[i,:] .= beta_col_a
    end

    return map_values, alpha_grid, beta_grid
end



"""
    interpolate_sphere_map(spheremap::AbstractMatrix{T},
                           alphas::AbstractVector{T},
                           betas::AbstractVector{T}) where T

Yields the interpolation of the 'spheremap' matrix given the angles 'alphas' and
'betas' in degrees.

"""
function interpolate_sphere_map(spheremap::AbstractMatrix{T},
    alphas::AbstractVector{T},
    betas::AbstractVector{T}) where T
    
    # close the sphere in alpha direction
    if !(T(360.) in alphas)
        alphas = vcat(alphas, T(360.))
        spheremap = vcat(spheremap, spheremap[1,:]')
    end

    return interpolate((alphas, betas), spheremap, (Gridded(Linear()),Gridded(Linear())))
end


"""
    spher_to_cart_coord(alpha::Real,
                        beta::Real, 
                        r::Real = 1.0)

Converts spherical coordinates (in degrees) into cartesian coordinates.

---
    spher_to_cart_coord(spherecoord::SphereCoord)

Uses 'spher_to_cart_coord' on a 'SphereCoord'.

"""
function spher_to_cart_coord(alpha::Real, 
    beta::Real, 
    r::Real = 1.0)

    return [r*sind(beta)*cosd(alpha), r*sind(beta)*sind(alpha), r*cosd(beta)]
end

spher_to_cart_coord(spherecoord::SphereCoord) = spher_to_cart_coord(spherecoord.alpha, 
                                                                    spherecoord.beta, 
                                                                    spherecoord.r)



"""
    cart_to_sphe_coord(x::Real, 
                       y::Real,
                       z::Real)

Converts cartesian coordinates into spherical coordinates (in degrees).

"""
function cart_to_sphe_coord(x::Real, 
    y::Real,
    z::Real)

    r = hypot(x, y, z)
    alpha = mod(atand(y, x), 360.)
    beta  = acosd(z / r)

    return (alpha, beta, r)
end



"""
    rot_mat(alpha::Real,
            beta::Real)

yields the 3D rotation matrix given two spherical coordinates in degrees.

---
    rot_mat(spherecoord::SphereCoord)

Uses 'rot_mat' on a 'SphereCoord'.

"""
function rot_mat(alpha::Real,
    beta::Real)

    # rotation of alpha around Z
    R_z_alpha = [ cosd(alpha) -sind(alpha) 0
                  sind(alpha) cosd(alpha)  0
                        0          0     1 ]
    # rotation of beta around Y
    R_y_beta = [ cosd(beta)  0  sind(beta)
                      0     1       0      
                -sind(beta)  0  cosd(beta) ]
                      
    return R_z_alpha * R_y_beta
end

rot_mat(spherecoord::SphereCoord) = rot_mat(spherecoord.alpha, spherecoord.beta)



"""
    pass_frame_to_frame(alpha::Real,
                        beta::Real,
                        Gamma::Real,
                        Psi::Real;
                        pre_load_rot_mat::Union{Matrix,Nothing} = nothing,
                        inverse::Bool = false)

Given the spherical coordinates (α,β) in degrees of a point in the ('X','Y','Z')
frame, yields the spherical coordinates (ϕ,θ) in degrees of the same point in
the rotated frame ('p','g','b') defined by ('Gamma','Psi'), also in degrees.

        Z                                   b     p (ϕ,θ)
        |    b                              |    x
        | Ψ x     p (α,β)                   | θ /'
        |--/'    x                          |--/ '
        | / '                               | /  '
        |/__'_______Y         =>            |/___'______g
       / `  '                              / `   '
      /---` '                             /---`  '
     /  Γ  `'                            /  ϕ  ` '
    X                                   /       `'
                                       p

The rotation matrix of passage to antenna frame can be pre-computed and given as
an argument 'pre_load_rot_mat' to avoid redundant computations if many points
are transformed with the same rotation angles.

If 'inverse' is true, the rotation matrix is transposed to go from the new
frame to the old frame instead.

---
    pass_frame_to_frame(spherecoord::SphereCoord,
                        new_spherecoord::SphereCoord;
                        kwds...)

Uses 'pass_frame_to_frame' on 'SphereCoord' structures.

"""
function pass_frame_to_frame(alpha::Real,#TODO: vectorize
    beta::Real,
    Gamma::Real,
    Psi::Real;
    pre_load_rot_mat::Union{Matrix,Nothing} = nothing,
    inverse::Bool = false)

    # rotation matrix for new frame
    R = isnothing(pre_load_rot_mat) ? rot_mat(Gamma, Psi) : pre_load_rot_mat
    
    # cartesian coord of point in (X,Y,Z)
    p_XYZ = spher_to_cart_coord(alpha, beta)

    # cartesian coord of point in (p,g,b)
    obj_pgb = (inverse ? R : R') * p_XYZ

    # spherical coord of point in (p,g,b)
    obj_pgb_sph = cart_to_sphe_coord(obj_pgb[1], obj_pgb[2], obj_pgb[3])

    return obj_pgb_sph[1], obj_pgb_sph[2]
end

function pass_frame_to_frame(spherecoord::SphereCoord,
    rot_spherecoord::SphereCoord;
    kwds...)

    new_spherecoord = pass_frame_to_frame(spherecoord.alpha, spherecoord.beta, 
                                          rot_spherecoord.alpha, rot_spherecoord.beta;
                                          kwds...)

    return SphereCoord(new_spherecoord[1], new_spherecoord[2], spherecoord.r)
end

end # module

