import numpy as np


def calculate_fd_matrix(x, risk_free_rate, mu, volatility, dt, theta):
    j = np.arange(len(x))
    alpha = 0.5 * dt * theta * (volatility**2 * j**2 - mu * j)
    beta = 1 - dt * theta * (volatility**2 * j**2 + risk_free_rate)
    kappa = 0.5 * dt * theta * (volatility**2 * j**2 + mu * j)
    return np.array([alpha, beta, kappa])