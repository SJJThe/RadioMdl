
using LinearAlgebra
using PyPlot
const plt = PyPlot
using Test

using Revise
using RadioMdl
using RadioMdl.CoordFrames


## Basic Tests 
"""
Test suite for spherical coordinate transformations
"""

println("="^60)
println("Testing Spherical Coordinate Transformations")
println("="^60)


# Test 1: Spherical to Cartesian - Cardinal points
println("\n[Test 1] Spherical to Cartesian - Cardinal Points")
println("-"^60)

# Point on +Z axis (α=0 °, β=0 °)
p1 = spher_to_cart_coord(0, 0, 1)
println("(α=0 °, β=0 °, r=1) -> ", p1)
@test p1 ≈ [0, 0, 1] atol=1e-10
println("✓ Should be [0, 0, 1] - Point on +Z axis")

# Point on +X axis (α=0 °, β=90 °)
p2 = spher_to_cart_coord(0, π/2, 1)
println("(α=0 °, β=90 °, r=1) -> ", p2)
@test p2 ≈ [1, 0, 0] atol=1e-10
println("✓ Should be [1, 0, 0] - Point on +X axis")

# Point on +Y axis (α=90 °, β=90 °)
p3 = spher_to_cart_coord(π/2, π/2, 1)
println("(α=90 °, β=90 °, r=1) -> ", p3)
@test p3 ≈ [0, 1, 0] atol=1e-10
println("✓ Should be [0, 1, 0] - Point on +Y axis")

# Point on -Z axis (α=0 °, β=180 °)
p4 = spher_to_cart_coord(0, π, 1)
println("(α=0 °, β=180 °, r=1) -> ", p4)
@test p4 ≈ [0, 0, -1] atol=1e-10
println("✓ Should be [0, 0, -1] - Point on -Z axis")


# Test 2: Cartesian to Spherical - Round trip
println("\n[Test 2] Cartesian to Spherical - Round Trip")
println("-"^60)

test_points = [
    (0, 0, 1),      # +Z
    (1, 0, 0),      # +X
    (0, 1, 0),      # +Y
    (0, 0, -1),     # -Z
    (1, 1, 0),      # XY diagonal
    (1, 1, 1),      # Space diagonal
]

for (x, y, z) in test_points
    println("\nOriginal: ($x, $y, $z)")
    
    # Convert to spherical
    α, β, r = cart_to_sphe_coord(x, y, z)
    println("  Spherical: α=$(rad2deg(α)) °, β=$(rad2deg(β)) °, r=$r")
    
    # Convert back to cartesian
    x2, y2, z2 = spher_to_cart_coord(α, β, r)
    println("  Back to Cartesian: ($x2, $y2, $z2)")
    
    @test [x2, y2, z2] ≈ [x, y, z] atol=1e-10
    println("  ✓ Round trip successful")
end


# Test 3: Rotation Matrix
println("\n[Test 3] Rotation Matrix Properties")
println("-"^60)

# Identity rotation (α=0, β=0)
R_id = rot_mat(0, 0)
println("Rotation matrix for (α=0 °, β=0 °):")
println(R_id)
@test R_id ≈ I(3) atol=1e-10
println("✓ Should be identity matrix")

# Rotation of 90 ° around Z
R_z90 = rot_mat(π/2, 0)
println("\nRotation matrix for (α=90 °, β=0 °):")
println(R_z90)
test_point = [1, 0, 0]
rotated = R_z90 * test_point
println("Point [1,0,0] rotated -> ", rotated)
@test rotated ≈ [0, 1, 0] atol=1e-10
println("✓ [1,0,0] should become [0,1,0]")


# Test 4: Frame Transformation
println("\n[Test 4] Frame Transformation - pass_frame_to_frame")
println("-"^60)

# Test 4a: No rotation (Γ=0 °, Ψ=0 °)
println("\nNo rotation (Γ=0 °, Ψ=0 °):")
α_in, β_in = 45.0, 60.0
α_out, β_out = pass_frame_to_frame(α_in, β_in, 0.0, 0.0)
println("Input: (α=$α_in °, β=$β_in °)")
println("Output: (α=$α_out °, β=$β_out °)")
@test α_out ≈ α_in atol=1e-8
@test β_out ≈ β_in atol=1e-8
println("✓ Should remain unchanged")

# Test 4b: Point on +Z axis, rotate frame 90 ° around Z
println("\nPoint on +Z (α=0 °, β=0 °), rotate frame Γ=90 ° around Z:")
α_out, β_out = pass_frame_to_frame(0.0, 0.0, 90.0, 0.0)
println("Output: (α=$α_out °, β=$β_out °)")
@test β_out ≈ 0.0 atol=1e-8
println("✓ Point on Z-axis should remain on Z-axis (β should stay 0 °)")

# Test 4c: Point on +X axis, rotate frame 90 ° around Z
println("\nPoint on +X (α=0 °, β=90 °), rotate frame Γ=90 ° around Z:")
α_out, β_out = pass_frame_to_frame(0.0, 90.0, 90.0, 0.0)
println("Output: (α=$α_out °, β=$β_out °)")
@test β_out ≈ 90.0 atol=1e-8
@test abs(α_out - (-90.0)) < 1e-8 || abs(α_out - 270.0) < 1e-8
println("✓ Point should rotate 90 ° in azimuth")

# Test 4d: Point on +Z axis, tilt frame 90 ° 
println("\nPoint on +Z (α=0 °, β=0 °), tilt frame Ψ=90 °:")
α_out, β_out = pass_frame_to_frame(0.0, 0.0, 0.0, 90.0)
println("Output: (α=$α_out °, β=$β_out °)")
@test β_out ≈ 90.0 atol=1e-8
println("✓ +Z axis should appear at β=90 ° in tilted frame")

println("\n" * "="^60)
println("All tests completed!")
println("="^60)



## Pattern Transformation Tests
println("="^60)
println("Testing Pattern Transformation with pass_frame_to_frame")
println("="^60)

# Helper function to create test patterns
function create_test_pattern(alphas, betas, pattern_type::Symbol)
    n_alpha = length(alphas)
    n_beta = length(betas)
    pattern = Matrix{Float64}(undef, n_alpha, n_beta)
    
    if pattern_type == :constant
        # Constant pattern
        pattern .= 1.0
        
    elseif pattern_type == :azimuthal
        # Pattern that varies only with alpha (azimuth)
        for i in 1:n_beta, j in 1:n_alpha
            pattern[j, i] = cosd(alphas[j])
        end
        
    elseif pattern_type == :elevation
        # Pattern that varies only with beta (elevation)
        for i in 1:n_beta, j in 1:n_alpha
            pattern[j, i] = sind(betas[i])
        end
        
    elseif pattern_type == :dipole
        # Simple dipole pattern: sin²(beta)
        for i in 1:n_beta, j in 1:n_alpha
            pattern[j, i] = sind(betas[i])^2
        end
        
    elseif pattern_type == :directional
        # Directional beam pointing at alpha=0°, beta=90°
        for i in 1:n_beta, j in 1:n_alpha
            # Gaussian-like beam
            delta_alpha = min(abs(alphas[j] - 0), abs(alphas[j] - 360))
            delta_beta = abs(betas[i] - 90)
            pattern[j, i] = exp(-(delta_alpha^2 + delta_beta^2) / 500)
        end
    end
    
    return SphereMap(pattern, alphas, betas)
end
#=
alpha_grid = collect(0:1.:360)
beta_grid = collect(0:1.:180)
pattern_test = create_test_pattern(alpha_grid, beta_grid, :directional)
map_val = reduce(hcat, [[pattern_test(alpha_grid[a],beta_grid[b]) 
                         for b in eachindex(beta_grid)] 
                        for a in eachindex(alpha_grid)])
fig, axs = plt.subplots(subplot_kw=Dict("projection"=>"polar"))
contourf(deg2rad.(alpha_grid), beta_grid, map_val; cmap="viridis")
=#

# Test 1: Constant Pattern (should remain constant after any rotation)
println("\n[Test 1] Constant Pattern - Invariance")
println("-"^60)

alphas = collect(0.0:10.0:350.0)
betas = collect(0.0:10.0:180.0)
pattern_const = create_test_pattern(alphas, betas, :constant)

println("Original pattern: all values = 1.0")
println("Applying rotation: Γ=45 °, Ψ=30 °")

pattern_transformed = pass_frame_to_frame(pattern_const, 45.0, 30.0)

println("Min value after transform: $(minimum(pattern_transformed))")
println("Max value after transform: $(maximum(pattern_transformed))")
@test all(isapprox.(pattern_transformed, 1.0, atol=1e-6))
println("✓ Constant pattern remains constant after rotation")


# Test 2: Dipole Pattern - Symmetry
println("\n[Test 2] Dipole Pattern - Azimuthal Symmetry")
println("-"^60)

pattern_dipole = create_test_pattern(alphas, betas, :dipole)

println("Rotating around Z-axis (Γ=90 °, Ψ=0 °)")
println("Dipole pattern should remain unchanged (azimuthally symmetric)")

pattern_rotZ = pass_frame_to_frame(pattern_dipole, 90.0, 0.0)

pat_diff = maximum(abs.(pattern_rotZ - pattern_dipole.spheremap))
println("Maximum difference: $pat_diff")
@test pat_diff < 0.1  # Allow some interpolation error
println("✓ Azimuthally symmetric pattern unchanged by Z-rotation")


# Test 3: Identity Transformation
println("\n[Test 3] Identity Transformation (Γ=0 °, Ψ=0 °)")
println("-"^60)

pattern_dir = create_test_pattern(alphas, betas, :directional)

println("Applying no rotation to directional pattern")
pattern_identity = pass_frame_to_frame(pattern_dir, 0.0, 0.0)

pat_diff = maximum(abs.(pattern_identity - pattern_dir.spheremap))
println("Maximum difference: $pat_diff")
@test pat_diff < 0.01
println("✓ Identity transformation preserves pattern")


# Test 4: Directional Beam Rotation
println("\n[Test 4] Directional Beam - Beam Peak Movement")
println("-"^60)

pattern_beam = create_test_pattern(alphas, betas, :directional)

# Find original peak location
idx = argmax(pattern_beam.spheremap)
i_max, j_max = Tuple(idx)
println("Original peak at: α=$(alphas[i_max]) °, β=$(betas[j_max]) °")
println("Original peak value: $(pattern_beam.spheremap[idx])")

# Rotate 90° around Z-axis
println("\nRotating 90° around Z-axis (Γ=90 °, Ψ=0 °)")
pattern_rot90 = pass_frame_to_frame(pattern_beam, 90.0, 0.0)

idx_new = argmax(pattern_rot90)
i_new, j_new = Tuple(idx_new)
println("New peak at approximately: α=$(alphas[i_new]) °, β=$(betas[j_new]) °")
println("New peak value: $(pattern_rot90[idx_new])")

# Peak should move in azimuth by ~90°
expected_alpha = mod(0.0 - 90.0, 360.0)
@test abs(alphas[i_new] - expected_alpha) < 20.0  # Within 20° due to discretization
println("✓ Beam peak moved approximately 90 ° in azimuth")


# Test 5: Energy Conservation
println("\n[Test 5] Energy Conservation")
println("-"^60)

pattern_test = create_test_pattern(alphas, betas, :directional)

println("Computing total energy (sum of pattern²)")
energy_original = sum(pattern_test.spheremap.^2)
println("Original energy: $energy_original")

pattern_rot = pass_frame_to_frame(pattern_test, 37.0, 25.0)

energy_rotated = sum(pattern_rot.^2)
println("Rotated energy: $energy_rotated")

relative_diff = abs(energy_rotated - energy_original) / energy_original
println("Relative difference: $(relative_diff * 100)%")
@test relative_diff < 0.05  # Within 5%
println("✓ Total energy approximately conserved")


# Test 6: Double Rotation (should approximate back to original)
println("\n[Test 6] Rotation and Inverse Rotation")
println("-"^60)

pattern_test = create_test_pattern(alphas, betas, :directional)
Gamma, Psi = 45.0, 30.0

println("Applying rotation: Γ=$Gamma °, Ψ=$Psi °")
pattern_rot = pass_frame_to_frame(pattern_test, Gamma, Psi)

println("Applying inverse rotation: Γ=$(-Gamma) °, Ψ=$(-Psi) °")
pattern_back = pass_frame_to_frame(pattern_rot, -Gamma, -Psi)

pat_diff = sqrt(sum((pattern_back - pattern_test.spheremap).^2) / 
           length(pattern_test.spheremap))
println("RMS difference from original: $pat_diff")
@test pat_diff < 0.15  # Allow for interpolation errors
println("✓ Rotation + inverse approximately recovers original")


# Test 7: Different Grid Sizes
println("\n[Test 7] Resample to Different Grid Size")
println("-"^60)

alphas_new = collect(0.0:15.0:345.0)
betas_new = collect(0.0:15.0:180.0)

pattern_orig = create_test_pattern(alphas, betas, :elevation)

println("Original grid: $(length(alphas))×$(length(betas))")
println("New grid: $(length(alphas_new))×$(length(betas_new))")

pattern_resampled = pass_frame_to_frame(pattern_orig, 60.0, 45.0, alphas_new, betas_new)

println("Output size: $(size(pattern_resampled))")
@test size(pattern_resampled) == (length(alphas_new), length(betas_new))
@test all(isfinite.(pattern_resampled))
println("✓ Resampled pattern has correct dimensions and finite values")


# Test 8: Elevation Pattern Rotation
println("\n[Test 8] Elevation Pattern - Tilt Effect")
println("-"^60)

pattern_elev = create_test_pattern(alphas, betas, :elevation)

println("Pattern = sin(β), varies only with β (elevation)")
idx_90 = findfirst(x -> x == 90.0, betas)
println("Original value at β=0°: $(pattern_elev.spheremap[1, 1]) \
         (should be ≈0, sin(0°)=0)")
println("Original value at β=90°: $(pattern_elev.spheremap[idx_90, 1]) \
         (should be ≈1, sin(90°)=1)")

# Tilt by 90° - what was at β=90° should now appear at β=0°
println("\nTilting frame by 90° (Γ=0°, Ψ=90°)")
println("β=0° in new frame should show what was at β=90° in old frame")
pattern_tilted = pass_frame_to_frame(pattern_elev, 0.0, 90.0)

println("New value at β=0°: $(pattern_tilted[1, 1])")
println("Expected: ≈1 (from sin(90°)=1)")

# Check that β=0° in new frame has high value (from old β=90°)
@test pattern_tilted[1, 1] > 0.8  # Should be close to sin(90°) = 1
println("✓ Tilted pattern shows expected redistribution")

println("\n" * "="^60)
println("All Pattern Transformation Tests Completed!")
println("="^60)


#=
## Realistic pattern tests
# Load antenna pattern
file_pattern_path = "/home/samthe/Documents/Data/Radio/Westford_Antenna/\
                     Gain_pattern_Ku_band/single_cut_res.cut"

pattern = RadioMdl.power_pattern_from_cut_file(file_pattern_path)

fig, axs = plt.subplots(subplot_kw=Dict("projection"=>"polar"))
contourf(deg2rad.(unique(pattern.az)), unique(pattern.polar),
         10. .*log10.(reshape(pattern.power, (length(unique(pattern.polar)),
                                 length(unique(pattern.az))))); cmap="viridis")
colorbar(label="Power (dBW)")

# Create SphereMap
smap = SphereMap(pattern; alpha_col=:az, beta_col=:polar, map_col=:power)

alpha_grid = 0:1:360
beta_grid = 0:1:180
map_val = reduce(hcat, [[smap(alpha_grid[a],beta_grid[b]) 
                         for b in eachindex(beta_grid)] 
                        for a in eachindex(alpha_grid)])
fig, axs = plt.subplots(subplot_kw=Dict("projection"=>"polar"))
contourf(deg2rad.(alpha_grid), beta_grid, 10. .*log10.(map_val); cmap="viridis")
colorbar(label="Power (dBW)")

# Pass pattern to other frame
alpha_grid_new = 0:1.:360
beta_grid_new = 0:1.:180
other_smap = pass_frame_to_frame(smap, 135., 45., alpha_grid_new, beta_grid_new)

fig, axs = plt.subplots(subplot_kw=Dict("projection"=>"polar"))
contourf(deg2rad.(alpha_grid_new), beta_grid_new, 10. .*log10.(other_smap');
         cmap="viridis")
colorbar(label="Power (dBW)")
=#

