"""
Plots the peak fluid-load fraction and the proportion of engaged
fibres as a function of the initial swelling ratio.  Creates
Fig 9
"""

import ucompress as uc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colours
from fihy_params.fihy71 import FiHy71
from scipy.integrate import quad

plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    "font.size": 24,
    "figure.autolayout": True
    })


f = lambda x, lam_m: max(60 * x ** 2 * (x - 1) * (x - lam_m) / (-3 * lam_m ** 5 + 5 * lam_m ** 4 - 5 * lam_m + 3), 0)

pars = FiHy71(nondim=False)
pars.physical['phi_0'] = 0
pars.physical['psi_0'] = 0

mech = uc.mechanics.FibreRecruitment(pars, distribution='quartic', homogeneous=True)
os = uc.osmosis.FloryHuggins()
model = uc.base_models.Hydrogel(mech, None, os, pars)

NL = 10
NJ = 30

fig_flf, ax_flf = plt.subplots(1,1, figsize=(7, 6))
fig_prop, ax_prop = plt.subplots(1,1, figsize = (7, 6))

ls = ['-', '--', '-.']

"""
Fix lam_m and vary J
"""
pars.physical['lam_m'] = 2
J = np.logspace(np.log10(1.1), np.log10(100), NJ)
E_m = [4e3, 4e4, 4e5]
NE = len(E_m)
v = np.zeros(NJ)
flf = np.zeros(NJ)
pi = np.zeros(NJ)

for m in range(NE):
    pars.physical['E_m'] = E_m[m]
    print('---------------------------------------------------------')
    for n in range(NJ):

        # dehydate 
        pars.physical['beta_r'] = 1
        pars.physical['beta_z'] = 1
        pars.physical['phi_0'] = 0
        model.assign(pars)

        # calculate chi
        calc = uc.fitting.ChiCalculator(model, pars)
        chi, beta_r, beta_z, phi_0 = calc.solve(J[n])

        # compute proportion of engaged fibre
        v[n] = quad(f, 1, beta_r, args = (pars.physical["lam_m"], ))[0]

        # update params for the hydrated state
        pars.physical['chi'] = chi
        pars.physical['beta_r'] = beta_r
        pars.physical['beta_z'] = beta_z
        pars.physical['phi_0'] = phi_0
        model.assign(pars)

        # solve for the instantaneous response
        exp = uc.experiments.ForceControlled(model, pars)
        sol = exp.initial_response()
        flf[n] = sol.fluid_load_fraction

    label = f'$E_m = {int(E_m[m]/1000):d}$ kPa'
    ax_prop.plot(J, v, label = label, lw = 3, ls = ls[m])
    ax_flf.plot(J, flf, label = label, lw = 3, ls = ls[m])



fig_label = ['(a)', '(b)']
n = 0
for ax in [ax_flf, ax_prop]:
    ax.set_xlabel(r'$J_h$ (-)')
    n += 1 



ax_flf.set_ylim(0.2, 1)

ax_flf.set_ylabel(r'$-p A / F$ (-)')
ax_prop.set_ylabel('Fraction of engaged fibres')

ax_flf.legend()

# fig_flf.savefig('sensitivity_1.pdf')
# fig_prop.savefig('sensitivity_2.pdf')

plt.show()