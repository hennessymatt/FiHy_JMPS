"""
Solves the hydration problem and creates Fig 4 panels (e)--(h)
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

N = 30
alpha_range = np.linspace(0.2, 0.8, N)
# alpha_range = np.logspace(-2, np.log10(0.95), N)
J_h = np.zeros(N); alpha = np.zeros(N); beta = np.zeros(N); P = np.zeros(N); phi = np.zeros(N)



mech = {
    # "NH": uc.mechanics.NeoHookean()
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

pars = FiHy75(nondim=True)
os = uc.osmosis.FloryHuggins()
perm = uc.permeability.Constant()

plt.figure(figsize=(16, 4))
ax_J = plt.subplot(143)
ax_a = plt.subplot(141)
ax_b = plt.subplot(142)
ax_p = plt.subplot(144)

# fig_p, ax_p = plt.subplots()

# chi_m = 0.55
# chi_f = 1.07
# chi_f = chi_m
chi = 0.57

# Initial guess
alpha[-1] = 1.1
beta[-1] = 1.5

# Set the initial porosity to zero
pars.update("phi_0", 0)

for key in mech:

    model = uc.base_models.Hydrogel(mech[key]["model"], perm, os, pars)
    
    for n, alpha_f in enumerate(alpha_range):

        # chi = chi_m + (chi_f - chi_m) * exp(-(1-alpha_f) / 0.05)
        # print(f'alpha_f = {alpha_f:.2f}, chi = {chi:.2e}')
        pars.update("chi", chi)
        pars.update("Phi_f", alpha_f)
        model.assign(pars)

        hydration = uc.experiments.Hydration(model, pars)
        alpha[n], beta[n], phi_0 = hydration.steady_response(alpha[n-1], beta[n-1])
        J_h[n] = 1 / (1 - phi_0)
        phi[n] = phi_0
        P[n] = model.osmosis.eval_osmotic_pressure(J_h[n])
    
    style = mech[key]["style"]

    ax_a.plot(alpha_range, alpha, style, lw = 2)
    ax_b.plot(alpha_range, beta, style, lw = 2)
    ax_J.plot(alpha_range, phi, style, lw = 2)
    ax_p.plot(alpha_range, P * pars.scaling["stress"] / 1000, style, lw = 2)

fig_label = ['(e)', '(f)', '(g)', '(h)']
n = 0
for ax in [ax_a, ax_b, ax_J,ax_p]:
    ax.set_xlabel(r'$\Phi_f$')
    ax.set_xticks([0.2, 0.4, 0.6, 0.8])
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

plt.savefig('hydration_2.pdf')
# plt.show()