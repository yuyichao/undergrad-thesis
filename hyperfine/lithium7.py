#!/usr/bin/env python

from pylab import *
from scipy.optimize import newton

h = 6.62606957e-34
hbar = h / 2 / pi
mu_B = 9.27400968e-24
Delta_W = 803.5040866e6
g_L = 1
g_S = 2.002319
S = 0.5
L = 0
J = 0.5
g_J = g_S * (J * (J + 1) - L * (L + 1) + S * (S + 1)) / 2 / J / (J + 1)
g_I = -0.00118221306

def freq_F(si, m_F, B):
    x = mu_B * B * (g_J - g_I) / h / Delta_W
    return (-Delta_W / 2 / 4 + mu_B * g_I * m_F * B / h +
            si * Delta_W / 2 * sqrt(1 + 2 * m_F * x / (3 / 2 + .5) + x**2))

levels = {
    (2, 2): lambda B: freq_F(1, 2, B),
    (2, 1): lambda B: freq_F(1, 1, B),
    (2, 0): lambda B: freq_F(1, 0, B),
    (2, -1): lambda B: freq_F(1, -1, B),
    (2, -2): lambda B: freq_F(1, -2, B),
    (1, 1): lambda B: freq_F(-1, 1, B),
    (1, 0): lambda B: freq_F(-1, 0, B),
    (1, -1): lambda B: freq_F(-1, -1, B)
    }

def calc_field(state_e, state_g, freq, unc=None, init=10):
    init = init / 1e4
    state_e = levels[state_e]
    state_g = levels[state_g]
    trans = lambda B: state_e(B) - state_g(B)
    B0 = newton(lambda B: trans(B) - freq, init) * 1e4
    if unc is None:
        return B0
    Delta_B = abs(newton(lambda B: trans(B) - freq - unc, init) -
                  newton(lambda B: trans(B) - freq + unc, init)) * 1e4 / 2
    return B0, Delta_B
