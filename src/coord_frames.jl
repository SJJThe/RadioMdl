
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

export cart_to_sphe_coord,
       get_angle_grids,
       get_angle_resolution,
       interpolate_sphere_map,
       map_sphere_coords,
       pass_frame_to_frame,
       rot_mat,
       SphereMap,
       spher_to_cart_coord



# """
# """
# struct SphereCoords{T,V<:AbstractVector{T}}
#     alpha::V
#     beta::V
#     r::V

#     function SphereCoord(alpha::V,
#         beta::V,
#         r::V) where {T,V<:AbstractVector{T}}

#         @assert length(alpha) == length(beta) == length(r)

#         new(alpha, beta, r)
#     end
# end



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
```

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
        @assert length(alpha_grid) == size(spheremap, 1) == size(interp_map, 1) - 1
        @assert length(beta_grid) == size(spheremap, 2) == size(interp_map, 2)
        @assert eltype(interp_map) == eltype(spheremap)

        new{T,V,M,I}(alpha_grid, beta_grid, spheremap, interp_map)
    end
end

function SphereMap(spheremap::AbstractMatrix{T},
    alpha_grid::AbstractVector{T},
    beta_grid::AbstractVector{T}) where T

    # create the interpolator
    interp_map = interpolate_sphere_map(spheremap, alpha_grid, beta_grid)

    return SphereMap(alpha_grid, beta_grid, spheremap, interp_map)
end

function SphereMap(spheremap::AbstractDataFrame;
    kwds...) 

    map_values, alpha_grid, beta_grid = map_sphere_coords(spheremap; kwds...)

    return SphereMap(map_values, alpha_grid, beta_grid)
end

(SM::SphereMap)(alpha::Real, beta::Real) = SM.interp_map(alpha, beta)

get_angle_grids(SM::SphereMap) = (SM.alpha_grid, SM.beta_grid)

get_angle_resolution(SM::SphereMap) = (SM.alpha_grid[2] - SM.alpha_grid[1],
                                       SM.beta_grid[2] - SM.beta_grid[1])

Base.show(io::IO, SM::SphereMap{T}) where {T} = begin 
    print(io, "SphereMap{$T}:\n")
    print(io, "sample grid size: $(size(SM.spheremap))\n")
    print(io, "angles resolutions (degrees): $(get_angle_resolution(SM))")
end


"""
    pass_frame_to_frame(SM::SphereMap{T},
                        Gamma::Real,
                        Psi::Real,
                        new_alpha_grid::AbstractVector{T} = SM.alpha_grid, 
                        new_beta_grid::AbstractVector{T} = SM.beta_grid) where T

Yields a new 'SphereMap' structure where the 'spheremap' has been converted to
a new frame defined by the rotation angles 'Gamma' and 'Psi' (in degrees) in the
original frame.

It is possible to give new angle grids to resample the resulting map.

'new_alpha_grid' and 'new_beta_grid' can be 'AbstractVector's or
'AbstractRange's.

"""
function pass_frame_to_frame(SM::SphereMap{T},
    Gamma::Real,
    Psi::Real,
    new_alpha_grid::AbstractVector{T} = SM.alpha_grid, 
    new_beta_grid::AbstractVector{T} = SM.beta_grid) where T
    
    itp = SM.interp_map
    new_map = Matrix{T}(undef, length(new_alpha_grid), length(new_beta_grid))
    @inbounds for b in eachindex(new_beta_grid)
        @simd for a in eachindex(new_alpha_grid)
            # For each point in NEW frame, find where it came from in OLD frame
            alpha_orig, beta_orig = pass_frame_to_frame(new_alpha_grid[a], 
                                                        new_beta_grid[b], -Gamma, -Psi)

            # sample original pattern
            new_map[a,b] = itp(alpha_orig, beta_orig)
        end
    end

    return SphereMap(new_map, new_alpha_grid, new_beta_grid)
end

function  pass_frame_to_frame(SM::SphereMap{T},
    Gamma::Real,
    Psi::Real,
    new_alpha_grid::AbstractRange{T}, 
    new_beta_grid::AbstractRange{T}) where T
    
    return pass_frame_to_frame(SM, Gamma, Psi, collect(new_alpha_grid), 
                               collect(new_beta_grid))
end



"""
    map_sphere_coords(spheremap::AbstractDataFrame;
                      alpha_col::Symbol = :az,
                      beta_col::Symbol = :polar,
                      map_col::Symbol = :power)

Yields the 2D-matrix of 'spheremap', formed from the column 'map_col', with rows and
columns corresponding to the angles in 'alpha_col' and 'beta_col', respectively. 

"""
function map_sphere_coords(spheremap::AbstractDataFrame;
    alpha_col::Symbol = :az,
    beta_col::Symbol = :polar,
    map_col::Symbol = :power)
    
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
    
    # close the sphere in beta direction
    alphas_loop = vcat(alphas, T(360.))
    spheremap_loop = vcat(spheremap, spheremap[1,:]')

    return interpolate((alphas_loop, betas), spheremap_loop, Gridded(Linear()))
end


"""
    spher_to_cart_coord(alpha::Real,
                        beta::Real, 
                        r::Real = 1.0)

Converts spherical coordinates (in radians) into cartesian coordinates.

"""
function spher_to_cart_coord(alpha::Real, 
    beta::Real, 
    r::Real = 1.0)

    return [r*sin(beta)*cos(alpha), r*sin(beta)*sin(alpha), r*cos(beta)]
end



"""
    cart_to_sphe_coord(x::Real, 
                       y::Real,
                       z::Real)

Converts cartesian coordinates into spherical coordinates (in radians).

"""
function cart_to_sphe_coord(x::Real, 
    y::Real,
    z::Real)

    r = hypot(x, y, z)
    alpha = mod(atan(y, x), 2π) # ∈ (0, 2π]
    beta  = acos(z / r) # ∈ [0, π]

    return (alpha, beta, r)
end



"""
    rot_mat(alpha::Real,
            beta::Real)

yields the 3D rotation matrix given two spherical coordinates in radians.

"""
function rot_mat(alpha::Real,
    beta::Real)

    # rotation of alpha around Z
    R_z_alpha = [ cos(alpha) -sin(alpha) 0
                  sin(alpha) cos(alpha)  0
                        0          0     1 ]
    # rotation of beta around Y
    R_y_beta = [ cos(beta)  0  sin(beta)
                      0     1       0      
                -sin(beta)  0  cos(beta) ]
                      
    return R_z_alpha * R_y_beta
end



"""
    pass_frame_to_frame(alpha::Real,
                        beta::Real,
                        Gamma::Real,
                        Psi::Real)

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

"""
function pass_frame_to_frame(alpha::Real,
    beta::Real,
    Gamma::Real,
    Psi::Real)

    #convert to radians
    alpha = deg2rad(alpha)
    beta = deg2rad(beta)
    Gamma = deg2rad(Gamma)
    Psi = deg2rad(Psi)

    # rotation matrix for new frame
    R = rot_mat(Gamma, Psi)
    
    # cartesian coord of point in (X,Y,Z)
    p_XYZ = spher_to_cart_coord(alpha, beta)

    # cartesian coord of point in (p,g,b)
    obj_pgb = R' * p_XYZ

    # spherical coord of point in (p,g,b)
    obj_pgb_sph = cart_to_sphe_coord(obj_pgb[1], obj_pgb[2], obj_pgb[3])

    return rad2deg(obj_pgb_sph[1]), rad2deg(obj_pgb_sph[2])
end

end # module

