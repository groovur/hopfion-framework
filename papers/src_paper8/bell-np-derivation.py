import numpy as np

def random_unit_vectors(num_samples):
    """Full isotropic 3D sampling of unit vectors n_hat"""
    cos_theta = 2 * np.random.rand(num_samples) - 1  # uniform in cosθ [-1,1]
    sin_theta = np.sqrt(1 - cos_theta**2)
    phi = 2 * np.pi * np.random.rand(num_samples)
    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta
    return np.stack((x, y, z), axis=1)

def chiral_weights_nonpert(theta, lambda_, beta_base=5.0, gamma=1.0):
    """
    Non-perturbative weights with angle-dependent density feedback.
    
    Physical picture: Density feedback is stronger in high-suppression 
    regions (perpendicular, θ≈90°) where softening is most needed. 
    Aligned directions (θ≈0) already have low suppression (sin⁴θ≈0), 
    so require less density feedback.
    
    beta_local = beta_base / (1 + γ·cos(θ))
      → Larger at θ=90° (perpendicular needs more softening)
      → Smaller at θ=0° (aligned already soft)
    
    This represents the non-perturbative (unobserved) state where 
    density feedback from the knot itself (not measurement) provides 
    baseline softening.
    """
    cos_theta = np.cos(theta)
    beta_local = beta_base / (1 + gamma * cos_theta)
    
    # Suppression for ±1 outcomes
    s_plus = (np.sin(theta)**4) / (1 + beta_local)
    s_minus = (np.cos(theta)**4) / (1 + beta_local)
    
    # Boltzmann weights
    w_plus  = np.exp(-lambda_ * s_plus)
    w_minus = np.exp(-lambda_ * s_minus)
    
    # Normalize
    norm = w_plus + w_minus + 1e-12
    p_plus = w_plus / norm
    return p_plus

def compute_correlation_nonpert(alpha_deg, beta_deg, lambda_=5.0, beta_base=5.0, gamma=1.0, num_samples=100000):
    """
    Non-perturbative correlation with denominator-based local density feedback softening.
    Full 3D isotropic sampling.
    """
    alpha = np.deg2rad(alpha_deg)
    beta_angle = np.deg2rad(beta_deg)

    a_hat = np.array([np.cos(alpha), np.sin(alpha), 0.0])
    b_hat = np.array([np.cos(beta_angle), np.sin(beta_angle), 0.0])

    n_hat = random_unit_vectors(num_samples)

    cos_theta_A = np.dot(n_hat, a_hat)
    cos_theta_B = np.dot(n_hat, b_hat)
    theta_A = np.arccos(np.clip(cos_theta_A, -1.0, 1.0))
    theta_B = np.arccos(np.clip(cos_theta_B, -1.0, 1.0))

    # Local softening probabilities
    p_A_plus = chiral_weights_nonpert(theta_A, lambda_, beta_base, gamma)
    p_B_plus = chiral_weights_nonpert(theta_B, lambda_, beta_base, gamma)

    # Singlet-like anti-correlation expectation
    E = np.mean(p_A_plus * (1 - p_B_plus) + (1 - p_A_plus) * p_B_plus - p_A_plus * p_B_plus - (1 - p_A_plus) * (1 - p_B_plus))
    return E

def compute_chsh_nonpert(lambda_=5.0, beta_base=5.0, gamma=1.0, num_samples=100000):
    angles = [0, 22.5, 45, 67.5]

    # Sample n_hat once for rescale computation
    n_hat_sample = random_unit_vectors(num_samples)
    theta_sample = np.arccos(np.abs(n_hat_sample[:, 2]))  # z-component
    avg_sin4 = np.mean(np.sin(theta_sample)**4)
    rescale_factor = 1 / avg_sin4

    # Compute raw E for each pair
    E_ab   = compute_correlation_nonpert(angles[0], angles[1], lambda_, beta_base, gamma, num_samples)
    E_abp  = compute_correlation_nonpert(angles[0], angles[3], lambda_, beta_base, gamma, num_samples)
    E_apb  = compute_correlation_nonpert(angles[2], angles[1], lambda_, beta_base, gamma, num_samples)
    E_apbp = compute_correlation_nonpert(angles[2], angles[3], lambda_, beta_base, gamma, num_samples)

    chsh_raw = abs(E_ab + E_abp) + abs(E_apb - E_apbp)
    chsh_rescaled = chsh_raw * rescale_factor

    return chsh_raw, chsh_rescaled, rescale_factor

# Example usage & visualization
if __name__ == "__main__":
    print("Non-perturbative 3D CHSH (denominator feedback, N=100000)")
    print("beta_base | Raw CHSH | Rescaled CHSH | Dynamic Rescale Factor")
    beta_values = [0.1, 0.3, 0.5, 0.7, 1.0, 2.0]
    for beta in beta_values:
        chsh_raw, chsh_res, rescale = compute_chsh_nonpert(lambda_=5.0, beta_base=beta, gamma=1.0)
        print(f"{beta:8.2f} | {chsh_raw:.4f}   | {chsh_res:.4f}       | {rescale:.4f}")
