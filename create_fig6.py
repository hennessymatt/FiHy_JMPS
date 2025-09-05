"""
Plots the instantaneous and steady-state responses for Fig 6
"""

import ucompress as uc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colours
from fihy_params.fihy71 import FiHy71

cols = colours.TABLEAU_COLORS


plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    "font.size": 20,
    "figure.autolayout": True
    })


# applied forces
forces = np.array([-10, -5, -1])

# pre-allocate
exp_steady = np.zeros(3)
exp_instant = np.zeros(3)

dir = 'exp_data/'

# load exp data
for n in range(len(forces)):
    D = np.loadtxt(f'{dir}Fig4_v2-FiHy071-load_{-forces[n]}.0N.csv', delimiter = ',')
    eps = D[:,2]
    lam = 1 + eps

    exp_instant[n] = lam[0]
    exp_steady[n] = lam[-2]

"""
Create model
"""
pars = FiHy71(nondim=True)
mech = uc.mechanics.FibreRecruitment(distribution='quartic', homogeneous=True)
perm = uc.permeability.KozenyCarman()
os = uc.osmosis.FloryHuggins()
model = uc.base_models.Hydrogel(mech, perm, os, pars)

# set up mass fractions
all_psi = [0.8, 0.96]

# create plots
fig_0, ax_0 = plt.subplots()
fig_inf, ax_inf = plt.subplots()

model_F = np.linspace(0, 10, 100)

for psi, ls in zip(all_psi, ['-', '-.']):

    # update mass fraction and calculate porosity
    pars.update("psi_0", psi)
    pars = uc.fitting.PorosityCalculator(pars).solve(have_fibres=True, update_params=True)
    model.assign(pars)

    # calculate swelling ratio
    J_h = 1 / (1 - pars.physical["phi_0"])

    # compute Flory parameter and hydration state
    chi_calc = uc.fitting.ChiCalculator(model, pars)
    chi, beta_r, beta_z, phi = chi_calc.solve(J_0 = J_h)

    # update parameters
    pars.update("chi", chi)
    pars.update("beta_r", beta_r)
    pars.update("beta_z", beta_z)
    pars.update("phi_0", phi)
    model.assign(pars)    

    # pre-allocate
    model_initial = np.zeros(100)
    model_steady = np.zeros(100)

    # compute instant and steady responses
    for n in range(100):
        pars.update('F', -model_F[n])
        problem = uc.experiments.ForceControlled(model, pars)
        
        sol = problem.initial_response()
        model_initial[n] = sol.lam_z

        sol = problem.steady_response()
        model_steady[n] = sol.lam_z


    # plot simulation data
    label = f'$\phi_0 = {pars.physical["phi_0"]:.2f}$'
    ax_0.plot(model_F, model_initial, 'k', ls = ls, lw = 3, label = label) 
    ax_inf.plot(model_F, model_steady, 'k', ls = ls, lw = 3, label = label)


# plot experimental data
ax_0.plot(-forces, exp_instant, 'o', ms = 10, markeredgecolor = 'k', mew = 1)
ax_inf.plot(-forces, exp_steady, 's', ms = 10, markeredgecolor = 'k', mew = 1)

for ax in [ax_0, ax_inf]:
    ax.set_xlabel(r'$-F$ (N)')
    ax.set_ylabel(r'$\beta_z$ (-)')
    ax.legend()
    
plt.show()