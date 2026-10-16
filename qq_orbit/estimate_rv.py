import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Lunar gravitational parameter (km^3/s^2)
MU_MOON = 4902.8001

def dynamics(t, state):
    """Two-body problem dynamics equations
    Args:
        t: Time (not used in autonomous system)
        state: State vector [rx, ry, rz, vx, vy, vz]
    Returns:
        State derivative [vx, vy, vz, ax, ay, az]
    """
    r = state[:3]  # Position vector
    v = state[3:]  # Velocity vector
    r_norm = np.linalg.norm(r)  # Magnitude of position vector
    a = -MU_MOON * r / (r_norm ** 3)  # Gravitational acceleration
    return np.concatenate([v, a])

def predict_step(state, P, dt, Q):
    """EKF prediction step
    Args:
        state: Current state vector [rx, ry, rz, vx, vy, vz]
        P: Current state covariance matrix (6x6)
        dt: Time step (s)
        Q: Process noise covariance matrix (6x6)
    Returns:
        predicted_state: Predicted state vector
        predicted_P: Predicted covariance matrix
    """
    # Predict state using numerical integration
    sol = solve_ivp(dynamics, [0, dt], state, method='RK45')
    predicted_state = sol.y[:, -1]
    
    # Calculate Jacobian of the transition function
    r = predicted_state[:3]
    r_norm = np.linalg.norm(r)
    r_norm3 = r_norm ** 3
    r_norm5 = r_norm ** 5
    
    # Initialize 6x6 Jacobian matrix F
    F = np.eye(6)
    # Position-velocity coupling terms
    F[0, 3] = F[1, 4] = F[2, 5] = dt
    
    # Acceleration-position coupling terms (linearized)
    for i in range(3):
        for j in range(3):
            if i == j:
                F[i+3, j] = -MU_MOON * (r_norm**2 - 3*r[i]**2) / r_norm5 * dt
            else:
                F[i+3, j] = -MU_MOON * (-3*r[i]*r[j]) / r_norm5 * dt
    
    # Update covariance matrix
    predicted_P = F @ P @ F.T + Q
    return predicted_state, predicted_P

def update_step(predicted_state, predicted_P, z, R):
    """EKF update step
    Args:
        predicted_state: Predicted state vector from prediction step
        predicted_P: Predicted covariance matrix from prediction step
        z: Measurement vector (position) [rx, ry, rz]
        R: Measurement noise covariance matrix (3x3)
    Returns:
        updated_state: Updated state vector after measurement incorporation
        updated_P: Updated covariance matrix
    """
    # Observation model Jacobian matrix
    H = np.zeros((3, 6))
    H[0, 0] = H[1, 1] = H[2, 2] = 1.0  # Only position is measured
    
    # Calculate Kalman gain
    y = z - predicted_state[:3]  # Measurement residual
    S = H @ predicted_P @ H.T + R  # Residual covariance
    K = predicted_P @ H.T @ np.linalg.inv(S)  # Kalman gain
    
    # Update state and covariance
    updated_state = predicted_state + K @ y
    updated_P = (np.eye(6) - K @ H) @ predicted_P
    return updated_state, updated_P

def ekf_estimate_rv(measurements, dt=60):
    """Main EKF function to estimate position and velocity from measurements
    Args:
        measurements: Position measurements array with shape (n_samples, 3)
        dt: Sampling time interval (s), default is 1 second
    Returns:
        estimated_r: Estimated position array with shape (n_samples, 3)
        estimated_v: Estimated velocity array with shape (n_samples, 3)
    """
    n_samples = len(measurements)
    
    # Initial state and covariance
    initial_state = np.concatenate([measurements[0], [0.1, 0.1, 0.1]])  # [r0, v0_initial_guess]
    P = np.diag([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])  # Initial covariance
    
    # Noise parameters
    Q = np.diag([1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4])  # Process noise covariance
    R = np.diag([0.1**2, 0.1**2, 0.1**2])  # Measurement noise covariance
    
    # Store estimation results
    estimated_states = [initial_state]
    
    # Iterate through measurements
    for i in range(1, n_samples):
        # Prediction step
        predicted_state, predicted_P = predict_step(
            estimated_states[-1], P, dt, Q
        )
        
        # Update step
        updated_state, updated_P = update_step(
            predicted_state, predicted_P, measurements[i], R
        )
        
        estimated_states.append(updated_state)
        P = updated_P
    
    # Separate position and velocity from states
    estimated_states = np.array(estimated_states)
    return estimated_states[:, :3], estimated_states[:, 3:]


def main():
    """Demonstration of the EKF implementation"""
    estimated_r, estimated_v = ekf_estimate_rv(r[11000:12000:60,:])
    
    # 3. Visualize results
    fig, axs = plt.subplots(2, 3, figsize=(18, 10))
    
    # Plot position components
    plt.show()

if __name__ == "__main__":
    main()
