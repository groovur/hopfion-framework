import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def compute_gaussian_density(point, flux_lines, sigma=0.3):
    """
    Treat each flux line segment as a Gaussian source.
    Total density = sum of contributions from all sources.
    """
    rho = 0.0
    for line in flux_lines:
        # Each point on the line contributes
        for segment_point in line:#[::10]:  # Sample every 10th point
            r = np.linalg.norm(point - segment_point)
            rho += np.exp(-r**2 / (2*sigma**2))

    return rho

phi = (1 + np.sqrt(5)) / 2
phi6 = phi**6
beta = 1.0  # Coupling strength (not β·ρ)
num_points = 300
sigma_gauss = 0.3  # Gaussian width for density

lines = []
for i in range(3):
    t = np.linspace(-2, 2, num_points)
    # Base direction
    phi_base = i * 2 * np.pi / 3
    x = 0.8 * t * np.cos(phi_base)
    y = 0.8 * t * np.sin(phi_base)
    z = t
    # Add twist near center (t near 0)
    twist = 0.5 * np.exp(-t**2 / 0.5) * np.sin(6*t + i*2*np.pi/3)
    x += twist * np.sin(phi_base + np.pi/2)
    y += twist * np.cos(phi_base + np.pi/2)
    pair1 = np.column_stack((x, y, z))
    pair2 = np.column_stack((-x, -y, -z))
    lines.append(pair1)
    lines.append(pair2)

# Flatten for S_eff calculation
all_points = np.vstack(lines)

# Position-dependent density
r = np.sqrt(np.sum(all_points**2, axis=1) + 1e-10)

# Option 1: Density falls off from center (like mass distribution)
#rho = np.exp(-r**2 / 4.0)  # Gaussian density profile

# Option 2: Constant background (current approach)
rho = np.ones_like(r) * 0.1

# Option 3: Density from knot itself (distance to nearest flux line)
# (More complex - need to compute distance to all lines)

# Option 4: Gaussian Density
# Compute position-dependent density (Gaussian sum from flux lines)
#rho = np.array([compute_gaussian_density(p, lines, sigma=sigma_gauss) for p in all_points])
#rho = rho / rho.max() if rho.max() > 0 else rho  # normalize (avoid div0)

# Calculate angle θ (polar angle from z-axis)
cos_theta = all_points[:, 2] / r
theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))

# Suppression with position-dependent density
S_eff = (1 / phi6) * (np.sin(theta)**4) / (1 + beta * rho)

# Normalize for color mapping
S_norm = (S_eff - S_eff.min()) / (S_eff.max() - S_eff.min() + 1e-10)

# Plot
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

colors = plt.cm.coolwarm(S_norm)

start_idx = 0
for line in lines:
    end_idx = start_idx + num_points
    line_points = all_points[start_idx:end_idx]
    line_colors = colors[start_idx:end_idx]
    for j in range(num_points-1):
        ax.plot(line_points[j:j+2, 0], line_points[j:j+2, 1], line_points[j:j+2, 2],
                color=line_colors[j], linewidth=3, alpha=0.9)
    start_idx = end_idx

# Labels
pair_labels = ['Color 1', 'Color 1 (anti)', 'Color 2', 'Color 2 (anti)', 
               'Color 3', 'Color 3 (anti)']
for i, line in enumerate(lines):
    end_point = line[-1] * 1.1
    ax.text(end_point[0], end_point[1], end_point[2], pair_labels[i],
            color='black', fontsize=11, fontweight='bold')

sm = plt.cm.ScalarMappable(cmap='coolwarm', 
                           norm=plt.Normalize(S_eff.min(), S_eff.max()))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.6, aspect=10)
cbar.set_label('Normalized Suppression S_eff(θ, ρ)', fontsize=12)

ax.set_xlabel('X', fontsize=12)
ax.set_ylabel('Y', fontsize=12)
ax.set_zlabel('Z', fontsize=12)
ax.set_title('3-Flux Knot: S_eff(θ,ρ) = sin⁴θ / [φ⁶(1+βρ)]\n' + 
             f'β = {beta:.1f}, φ⁶ ≈ {phi6:.3f}', fontsize=14)
ax.view_init(elev=25, azim=60)

plt.tight_layout()
plt.show()
