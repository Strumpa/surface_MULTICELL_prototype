## Translation of 2D PIJ problem from A. Hébert's book MATLAB to Python

import numpy as np
from SYB2D import sybt2d, tij2d, sybrhl, indpos

def get_radii(region_volumes):
    """
    Calculate the radii of concentric circles that correspond to the given region volumes.
    Parameters:
        region_volumes (list): List of volumes for each region.
    Returns:
        list: List of radii for each region.
    """
    radii = []
    cumulative_volume = 0.0
    for vol in region_volumes:
        cumulative_volume += vol
        radius = np.sqrt(cumulative_volume / np.pi)
        if vol == region_volumes[-1]:
            return radii
        else:
            radii.append(radius)
    return radii

def al1eig(A, eps):
    """
    % find the fundamental eigenvalue and corresponding eigenvector of
    % equation (a-eval)*evect=0 using the power method.
    % function [iter,evect,eval] = al1eig(A,eps)
    % (c) 2020 Alain Hebert, Polytechnique Montreal
    """
    maxiter = 1000
    n = A.shape[0]
    evect = np.ones(n)
    eval = 0.0
    for iter in range(1, maxiter + 1):
        evect_new = np.dot(A, evect)
        eval_new = np.linalg.norm(evect_new)
        evect_new = evect_new / eval_new
        if np.linalg.norm(evect_new - evect) < eps and abs(eval_new - eval) < eps:
            return iter, evect_new, eval_new
        evect = evect_new
        eval = eval_new

    return iter, evect, eval


if __name__ == "__main__":

    reference_keff_from_book = 1.171670
    ## Problem 3.13 from A. Hebert's book p. 194
    a = np.sqrt(4.9)  # pincell side
    albedo = 1.0 # isotropic boundary condition
    nangle = 14 # number of angles
    ngauss = 4 # number of gauss quadrature points
    Vol = [0.4, 0.7, 0.4, 1.3, 2.1] # 2D volumes : cm^2
    sig_tot = [0.2, 0.0, 0.5, 0.0, 0.3] # total macroscopic cross sections : cm^-1
    sig_scattering = [0.05, 0.0, 0.05, 0.0, 0.05] # scattering macroscopic cross sections : cm^-1
    nu_sig_f = [1.4, 0.0, 0.0, 0.0, 0.0] # neutron production cross section = nu (avg number of neutrons per fission) times Sigma_f, fission macroscopic xs.

    S0 = np.diag(sig_scattering)
    Qfiss = np.diag(nu_sig_f)

    rad = get_radii(Vol)  # compute radii from volumes
    tracks = sybt2d(a, rad, nangle, ngauss)
    nsurf = int(tracks[0])
    nvol = int(tracks[1])
    surfaces = tracks[5:5+nsurf]
    volumes = tracks[5+nsurf:5+nsurf+nvol]
    print("2D CARCEL SYB Problem 3.13 from A. Hebert's book p.194")
    print("Number of Surfaces: ", nsurf)
    print("Number of Volumes: ", nvol)
    print("Surfaces: ", surfaces)
    print("Volumes: ", volumes)

    print("Computing Tij matrix...")
    Tij = tij2d(tracks, sig_tot)
    print("Tij computation complete!")
    print(f"Tij shape: {Tij.shape}")
    print(f"Tij sum BEFORE normalization: {np.sum(Tij)}")
    
    print("\nApplying Stamm'ler normalization...")
    Tij_normalized = sybrhl(tracks, sig_tot, Tij)
    print("Normalization complete!")
    print(f"Tij sum AFTER normalization: {np.sum(Tij_normalized)}")
    print(f"Max difference between original and normalized: {np.max(np.abs(Tij - Tij_normalized))}")
    
    # Optional: Extract probability sub-matrices (similar to PIJ_2D_CARCEL.m)
    print("\n--- Extracting Probability Sub-Matrices ---")
    
    # Reconstruct full symmetric matrix from packed format
    T_matrix = np.zeros((nsurf + nvol, nsurf + nvol))
    for i in range(1, nsurf + nvol + 1):
        for j in range(i, nsurf + nvol + 1):
            idx = indpos(i, j) - 1
            T_matrix[i-1, j-1] = Tij_normalized[idx]
            T_matrix[j-1, i-1] = Tij_normalized[idx]  # Symmetry
    
    # Verify symmetry
    if not np.allclose(T_matrix, T_matrix.T):
        print("WARNING: T_matrix is not symmetric!")
    else:
        print("T_matrix symmetry verified ✓")
    
    print(f"\nFull T_matrix (normalized):")
    print(T_matrix)
    print(f"T_matrix diagonal: {np.diag(T_matrix)}")
    
    # Extract sub-blocks
    t_SS = T_matrix[0:nsurf, 0:nsurf]  # Surface to surface
    t_Sv = T_matrix[0:nsurf, nsurf:nsurf+nvol]  # Surface to volume
    t_vS = T_matrix[nsurf:nsurf+nvol, 0:nsurf]  # Volume to surface  
    t_ij = T_matrix[nsurf:nsurf+nvol, nsurf:nsurf+nvol]  # Volume to volume
    
    print(f"t_SS shape: {t_SS.shape} - Surface-to-surface probabilities")
    print(f"t_Sv shape: {t_Sv.shape} - Surface-to-volume probabilities")
    print(f"t_vS shape: {t_vS.shape} - Volume-to-surface probabilities")
    print(f"t_ij shape: {t_ij.shape} - Volume-to-volume collision probabilities")
    
    # Check conservation for one region (row sum should be ≤ 1)
    print(f"\nConservation check for region 1:")
    region_idx = 0  # First region
    row_sum = np.sum(t_vS[region_idx, :]) + np.sum(t_ij[region_idx, :] * sig_tot)
    print(f"  Sum of probabilities from region 1: {row_sum:.6f} (should be ≤ 1.0)")

    # Recover probabilities : surface to surface 
    P_SS = np.zeros((nsurf, nsurf))
    for alpha in range(nsurf):
        for beta in range(nsurf):
            P_SS[alpha, beta] = t_SS[alpha, beta] * 4 / surfaces[alpha]
    print(f"\nRecovered surface-to-surface probabilities P_SS:")
    print(P_SS)
    print(f"P_SS sum check (row sums should be ≤ 1.0):")
    for alpha in range(nsurf):
        row_sum = np.sum(P_SS[alpha, :])
        print(f"  Row {alpha+1} sum: {row_sum:.6f}")
    # Recover probabilities : volume i to volume j
    p_ij = np.zeros((nvol, nvol))
    for i in range(nvol):
        for j in range(nvol):
            p_ij[i, j] = t_ij[i, j] / volumes[i]
    print(f"\nRecovered volume-to-volume probabilities p_ij:")
    print(p_ij)
    print(f"p_ij sum check (row sums should be ≤ 1.0):")
    for i in range(nvol):
        row_sum = np.sum(p_ij[i, :])
        print(f"  Row {i+1} sum: {row_sum:.6f}")

    # Recover probabilities : volume to surface
    P_vS = np.zeros((nvol, nsurf)) 
    for i in range(nvol):
        for alpha in range(nsurf):
            P_vS[i, alpha] = t_vS[i, alpha] / volumes[i]
    print(f"\nRecovered volume-to-surface probabilities P_vS:")
    print(P_vS)
    print(f"P_vS sum check (row sums should be ≤ 1.0):")
    for i in range(nvol):
        row_sum = np.sum(P_vS[i, :])
        print(f"  Row {i+1} sum: {row_sum:.6f}")
    
    # Recover probabilities : surface to volume
    p_Sv = np.zeros((nsurf, nvol))
    for alpha in range(nsurf):
        for j in range(nvol):
            p_Sv[alpha, j] = t_Sv[alpha, j] * 4 / surfaces[alpha]
    print(f"\nRecovered surface-to-volume probabilities p_Sv:")
    print(p_Sv)
    print(f"p_Sv sum check (row sums should be ≤ 1.0):")
    for alpha in range(nsurf):
        row_sum = np.sum(p_Sv[alpha, :])
        print(f"  Row {alpha+1} sum: {row_sum:.6f}")

    print("\n--- End of Probability Sub-Matrix Extraction ---")
    
    print("\n--- Computing Reduced Collision Probabilities ---")
    P_SS_tilde = albedo * np.linalg.inv(np.identity(nsurf) - albedo * P_SS)
    print(f"P_SS_tilde:\n{P_SS_tilde}")
    
    Pvv_tilde = p_ij + np.matmul(np.matmul(P_vS, P_SS_tilde), p_Sv)
    print(f"\nPvv_tilde:\n{Pvv_tilde}")
    print(f"Pvv_tilde row sums: {np.sum(Pvv_tilde, axis=1)}")

    print(f"\nS0 (scattering matrix):\n{S0}")
    print(f"Pvv_tilde * S0:\n{np.matmul(Pvv_tilde, S0)}")
    
    W = np.matmul(np.linalg.inv(np.identity(nvol) - np.matmul(Pvv_tilde, S0)), Pvv_tilde)

    print("\nFinal matrix W (Neutron transport operator):")
    print(W)
    print("W computation complete!")
    print("=" * 70)
    # Compute eigenvalue using power method
    A = np.matmul(W, Qfiss)
    eps = 1e-10
    iter, evect, eval = al1eig(A, eps)
    print(f"Fundamental eigenvalue (k-effective): {eval:.6f} found in {iter} iterations")
    print(f"Corresponding eigenvector (flux distribution):")
    print(evect)
    print("\n--- End of 2D CARCEL SYB Problem ---")
    
    

