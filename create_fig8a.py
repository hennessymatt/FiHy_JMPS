"""
This script computes the peak (instantaneous) FLF as a function
of the fibre fraction for different maximum recruitment stretches.
This creates Fig 8a.
"""

import ucompress as uc
import matplotlib.pyplot as plt
import numpy as np
from fihy_params import FiHy71

plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    "font.size": 24,
    "figure.autolayout": True
    })

N = 51
Phi_f = np.r_[0, np.linspace(0.01, 0.8, N-1)]
J_h = np.zeros(N); alpha = np.zeros(N); beta = np.zeros(N); P = np.zeros(N); phi = np.zeros(N)
flf = np.zeros(N); P = np.zeros(N); Szz = np.zeros(N); Pi = np.zeros(N); lam_r = np.zeros(N); lam_z = np.zeros(N)

lam_m = [1, 2, 3]

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

pars = FiHy71(nondim=True)
os = uc.osmosis.FloryHuggins()
perm = uc.permeability.Constant()

model_i = uc.base_models.Hydrogel(mech["No recruitment"]["model"], perm, os, pars)
model_a = uc.base_models.Hydrogel(mech["Recruitment"]["model"], perm, os, pars)

# set the initial guesses of the hydration deformation
alpha[-1] = 1.1
beta[-1] = 1.1

fig, ax_S = plt.subplots(1,1, figsize = (7, 6))
ls = ['-', '--', '-.']

for m in range(len(lam_m)):
    print('====================================================')
    print(lam_m[m])
    pars.update("lam_m", lam_m[m])

    if lam_m[m] == 1:
        model = model_i
    else:
        model = model_a

    for n in range(N):

        # Update the fibre fraction
        pars.update("Phi_f", Phi_f[n])

        # Set the initial state to be dehydrated
        pars.update("phi_0", 0)
        pars.update("beta_r", 1)
        pars.update("beta_z", 1)

        # update the model
        model.assign(pars)

        """
        Solve the hydration problem
        """
        exp = uc.experiments.Hydration(model, pars)
        alpha[n], beta[n], phi[n] = exp.steady_response(alpha[n-1], beta[n-1])

        pars.update("phi_0", phi[n])
        pars.update("beta_r", alpha[n])
        pars.update("beta_z", beta[n])

        model.assign(pars)

        """
        Solve for the instantaneous and equilibrium response
        """
        exp = uc.experiments.ForceControlled(model, pars)
        
        if n == 0:
            sol_0 = exp.initial_response(lam_z_0=0.1)
        else:
            sol_0 = exp.initial_response()


        # re-dim
        # sol_0.redimensionalise(pars)

        lam_r[n] = 1/np.sqrt(sol_0.lam_z)
        lam_z[n] = sol_0.lam_z

        _, _, Szz[n] = model.mechanics.eval_stress(lam_r[n], lam_r[n], lam_z[n])
        Pi[n] = model.osmosis.eval_osmotic_pressure(1)

        P[n] = sol_0.p[0]
        flf[n] = sol_0.fluid_load_fraction

    label = f'$\lambda_m = {lam_m[m]:.1f}$'

    F = -pars.physical["F"]
    mean_S = np.pi * Szz / F
    mean_P = np.pi * P * lam_r**2 / F
    mean_Pi = np.pi * Pi * lam_r**2 / F

    lw = 3; ms = 12
    ax_S.plot(Phi_f[1:], flf[1:], lw = lw, label = label, ls = ls[m])
    if m == len(lam_m) - 1:
        label = f'No fibres'
        ax_S.plot(Phi_f[0], mean_P[0] - mean_Pi[0], 'd', ms = ms, markeredgecolor = 'k', label = label)


"""
Plot customisation
"""

fig_label = ['(a)']
n = 0
for ax in [ax_S]:
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.annotate(fig_label[n], 
                xy=(0.05,0.9), 
                xycoords='axes fraction',
                bbox=dict(facecolor='1', edgecolor='none', pad=3.0)
                )
    n += 1 

ax_S.set_xlabel(r'$\Phi_f$ (-)')
ax_S.set_ylabel(r'$-p A / F$ (-)')


ax_S.legend()
# fig.savefig('peak_flf_1.pdf')
plt.show()