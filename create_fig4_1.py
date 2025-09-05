"""
Solves the hydration problem and creates Fig 4 panels (a)--(d)
"""

import ucompress as uc
import numpy as np
import matplotlib.pyplot as plt
from fihy_params import FiHy75

plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    "font.size": 20,
    "figure.autolayout": True
    })

N = 50
chi_range = np.linspace(0, 0.8, N)
J_h = np.zeros(N); alpha = np.zeros(N); beta = np.zeros(N); P = np.zeros(N); phi = np.zeros(N)

mech = {
    "No recruitment": {
        "model": uc.mechanics.FibreReinforced(homogeneous=True),
        "style": "-",
        "label": "No recruit."
    },
    "Recruitment": {
        "model": uc.mechanics.FibreRecruitment(distribution='quartic', homogeneous=True),
        "style": "--",
        "label": "Recruit."
    }
}

pars = FiHy75(nondim = True)
os = uc.osmosis.FloryHuggins()
perm = uc.permeability.Constant()

plt.figure(figsize=(16, 4))
ax_J = plt.subplot(143)
ax_a = plt.subplot(141)
ax_b = plt.subplot(142)
ax_p = plt.subplot(144)

# set the initial porosity to zero for hydration step
pars.update("phi_0", 0)

# set initial guess
alpha[-1] = 1.1
beta[-1] = 1.1

for key in mech:

    model = uc.base_models.Hydrogel(mech[key]["model"], perm, os, pars)
    
    for n, chi in enumerate(chi_range):

        pars.update("chi", chi)
        model.assign(pars)

        exp = uc.experiments.Hydration(model, pars)
        alpha[n], beta[n], phi_0 = exp.steady_response(alpha[n-1], beta[n-1])
        J_h[n] = 1 / (1 - phi_0)
        phi[n] = phi_0
        P[n] = model.osmosis.eval_osmotic_pressure(J_h[n])

    style = mech[key]["style"]
    label = mech[key]["label"]

    ax_a.plot(chi_range, alpha, style, lw = 2, label = label)
    ax_b.plot(chi_range, beta, style, lw = 2)
    ax_J.plot(chi_range, phi, style, lw = 2)
    ax_p.plot(chi_range, P * pars.scaling["stress"] / 1000, style, lw = 2)

fig_label = ['(a)', '(b)', '(c)', '(d)']
n = 0
for ax in [ax_a, ax_b, ax_J, ax_p]:
    ax.set_xlabel(r'$\chi$')
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.annotate(fig_label[n], 
                xy=(0.05,0.9), 
                xycoords='axes fraction',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0)
                )
    n += 1 


ax_a.set_ylabel(r'$\alpha_\parallel$ (-)')
ax_b.set_ylabel(r'$\alpha_\perp$ (-)')
ax_J.set_ylabel(r'$\phi_0$ (-)')
ax_p.set_ylabel(r'$\Pi$ (kPa)')

# ax_a.legend()

plt.savefig('hydration_1.pdf')
# plt.show()