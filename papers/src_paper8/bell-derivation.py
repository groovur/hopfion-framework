import numpy as np

def compute_w_plus(theta, f):
    """Probability weight for +1 outcome: [cos(θ/2)^4]^f"""
    return (np.cos(theta / 2)**4)**f

def compute_w_minus(theta, f):
    """Probability weight for -1 outcome: [sin(θ/2)^4]^f"""
    return (np.sin(theta / 2)**4)**f

def random_unit_vectors_mid_plane_bias(num_samples, bias_sigma=np.pi/6):
    """
    Generate unit vectors on the sphere with bias toward the xy-plane (equator, θ=π/2).
    bias_sigma: Gaussian width around θ=π/2 (smaller = stronger bias).
    """
    # Sample theta with Gaussian bias around π/2
    theta_mean = np.pi / 2
    theta = np.random.normal(loc=theta_mean, scale=bias_sigma, size=num_samples)

    # Wrap theta to [0, π]
    theta = np.mod(theta, np.pi)

    # Sample phi uniformly
    phi = 2 * np.pi * np.random.rand(num_samples)

    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.stack((x, y, z), axis=1)

def random_unit_vectors_no_plane_bias(num_samples):
    """Full isotropic sphere sampling"""
    cos_theta = 2 * np.random.rand(num_samples) - 1  # uniform in [-1,1]
    sin_theta = np.sqrt(1 - cos_theta**2)
    phi = 2 * np.pi * np.random.rand(num_samples)
    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta
    return np.stack((x, y, z), axis=1)

def compute_correlation_with_dynamic_rescale(alpha_deg, beta_deg, f_base=2.0, beta=5.0, num_samples=100000, dim=2, bias_sigma=np.pi/6):
    alpha = np.deg2rad(alpha_deg)
    beta_angle = np.deg2rad(beta_deg)

    a_hat = np.array([np.cos(alpha), np.sin(alpha), 0.0])
    b_hat = np.array([np.cos(beta_angle), np.sin(beta_angle), 0.0])

    # Generate n_hat depending on dim
    if dim == 2:
        # Plane-constrained: theta = π/2, phi uniform
        phi = 2 * np.pi * np.random.rand(num_samples)
        n_hat = np.stack((np.cos(phi), np.sin(phi), np.zeros(num_samples)), axis=1)
    elif dim == 3:
        # Full isotropic 3D
        n_hat = random_unit_vectors_no_plane_bias(num_samples)
    elif dim == 'biased':
        # 3D with mid-plane bias
        n_hat = random_unit_vectors_mid_plane_bias(num_samples, bias_sigma=bias_sigma)
    else:
        raise ValueError("dim must be 2, 3, or 'biased'")

    cos_theta_A = np.dot(n_hat, a_hat)
    cos_theta_B = np.dot(n_hat, b_hat)
    theta_A = np.arccos(np.clip(cos_theta_A, -1.0, 1.0))
    theta_B = np.arccos(np.clip(cos_theta_B, -1.0, 1.0))

    # Dynamic softening: higher f when theta small (aligned), but denom softens
    # Clip f to prevent extreme values
    f_A = np.clip(f_base / (1 + beta * np.cos(theta_A)), 0.1, 10.0)
    f_B = np.clip(f_base / (1 + beta * np.cos(theta_B)), 0.1, 10.0)

    w_A_plus  = compute_w_plus(theta_A, f_A)
    w_A_minus = compute_w_minus(theta_A, f_A)
    norm_A = w_A_plus + w_A_minus + 1e-12  # epsilon for zero norm
    p_A_plus = w_A_plus / norm_A

    w_B_plus  = compute_w_plus(theta_B, f_B)
    w_B_minus = compute_w_minus(theta_B, f_B)
    norm_B = w_B_plus + w_B_minus + 1e-12
    p_B_plus = w_B_plus / norm_B

    p_same = p_A_plus * (1 - p_B_plus) + (1 - p_A_plus) * p_B_plus
    p_diff = p_A_plus * p_B_plus + (1 - p_A_plus) * (1 - p_B_plus)

    E_raw = np.mean(p_same - p_diff)

    # Dynamic rescaling based on actual geometry
    sin4_theta = np.sin(theta_A)**4  # or average over A/B
    avg_sin4 = np.mean(sin4_theta)
    rescale_factor = 1 / avg_sin4 if avg_sin4 > 0.01 else 2.0
    
    E_rescaled = E_raw * rescale_factor
    return E_raw, E_rescaled


def compute_chsh(f_base=2.0, beta=5.0, num_samples=100000, dim=2, bias_sigma=np.pi/6):
    """
    Standard CHSH angles (degrees): a=0, a'=45, b=22.5, b'=67.5
    CHSH = |E(a,b) + E(a,b')| + |E(a',b) - E(a',b')|

    dim: 2 = plane-constrained
         3 = full isotropic 3D
         'biased' = 3D with mid-plane bias (use bias_sigma)
    """
    angles = [0, 22.5, 45, 67.5]

    E_ab_raw, E_ab_res = compute_correlation_with_dynamic_rescale(angles[0], angles[1], f_base=f_base, beta=beta, num_samples=num_samples, dim=dim, bias_sigma=bias_sigma)
    E_abp_raw, E_abp_res = compute_correlation_with_dynamic_rescale(angles[0], angles[3], f_base=f_base, beta=beta, num_samples=num_samples, dim=dim, bias_sigma=bias_sigma)
    E_apb_raw, E_apb_res = compute_correlation_with_dynamic_rescale(angles[2], angles[1], f_base=f_base, beta=beta, num_samples=num_samples, dim=dim, bias_sigma=bias_sigma)
    E_apbp_raw, E_apbp_res = compute_correlation_with_dynamic_rescale(angles[2], angles[3], f_base=f_base, beta=beta, num_samples=num_samples, dim=dim, bias_sigma=bias_sigma)

    chsh_rescaled = abs(E_ab_res + E_abp_res) + abs(E_apb_res - E_apbp_res)
    chsh_raw = abs(E_ab_raw + E_abp_raw) + abs(E_apb_raw - E_apbp_raw)

    return chsh_raw, chsh_rescaled


# Example usage - compare 2D, 3D random, 3D biased
print("CHSH comparisons (N=100000, f_base=2.0):")
for beta in [0.1, 0.3, 0.5, 0.7, 1.0]:
    print(f"\nbeta = {beta:.1f}")
    
    # 2D plane-constrained
    chsh_raw_2d, chsh_res_2d = compute_chsh(f_base=2.0, beta=beta, dim=2)
    print(f"  2D plane: Raw = {chsh_raw_2d:.3f}, Rescaled = {chsh_res_2d:.3f}")
    
    # 3D full isotropic
    chsh_raw_3d, chsh_res_3d = compute_chsh(f_base=2.0, beta=beta, dim=3)
    print(f"  3D random: Raw = {chsh_raw_3d:.3f}, Rescaled = {chsh_res_3d:.3f}")
    
    # 3D mid-plane biased (sigma = π/6 ~30) 
    chsh_raw_bias, chsh_res_bias = compute_chsh(f_base=2.0, beta=beta, dim='biased', bias_sigma=np.pi/6)
    print(f"  3D biased (σ=π/6): Raw = {chsh_raw_bias:.3f}, Rescaled = {chsh_res_bias:.3f}")

    # 3D mid-plane biased (sigma = π/12 ~15) 
    chsh_raw_bias, chsh_res_bias = compute_chsh(f_base=2.0, beta=beta, dim='biased', bias_sigma=np.pi/12)
    print(f"  3D biased (σ=π/12): Raw = {chsh_raw_bias:.3f}, Rescaled = {chsh_res_bias:.3f}")
    
    # 3D mid-plane biased (sigma = π/18 ~10) 
    chsh_raw_bias, chsh_res_bias = compute_chsh(f_base=2.0, beta=beta, dim='biased', bias_sigma=np.pi/18)
    print(f"  3D biased (σ=π/18): Raw = {chsh_raw_bias:.3f}, Rescaled = {chsh_res_bias:.3f}")

