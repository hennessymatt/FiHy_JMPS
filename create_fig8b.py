"""
This script computes the peak (instantaneous) FLF as a function
of the fibre fraction for different maximum recruitment stretches.
Creates Fig 8b.
"""

import ucompress as uc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colours
from fihy_params import FiHy71

cols = colours.TABLEAU_COLORS

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

pars = FiHy71(nondim = True)
os = uc.osmosis.FloryHuggins()
perm = uc.permeability.Constant()

# update the matrix stiffness
pars.update("E_m", 400 * 1e3)

model_i = uc.base_models.Hydrogel(mech["No recruitment"]["model"], perm, os, pars)
model_a = uc.base_models.Hydrogel(mech["Recruitment"]["model"], perm, os, pars)

# set the initial guesses of the hydration deformation
alpha[-1] = 1.1
beta[-1] = 1.1

fig, ax_S = plt.subplots(1,1, figsize = (7, 6))
ls = ['-', '--', '-.']
for m,c in zip(range(len(lam_m)), cols):
    print('====================================================')
    print(lam_m[m])
    pars.update("lam_m", lam_m[m])

    if lam_m[m] == 1:
        model = model_i
    else:
        model = model_a

    for n in range(N):

        print('updating parameters', end = '...')

        # Update the fibre fraction
        pars.update("Phi_f", Phi_f[n])

        # Set the initial state to be dehydrated
        pars.update("phi_0", 0)
        pars.update("beta_r", 1)
        pars.update("beta_z", 1)

        # update the model
        model.assign(pars)
        print('done')
        """
        Solve the hydration problem
        """
        print('solving hydration problem')
        exp = uc.experiments.Hydration(model, pars)
        alpha[n], beta[n], phi[n] = exp.steady_response(alpha[n-1], beta[n-1])
        print('done')

        print('updating parameters', end = '...')
        pars.update("phi_0", phi[n])
        pars.update("beta_r", alpha[n])
        pars.update("beta_z", beta[n])

        model.assign(pars)
        print('done')

        """
        Solve for the instantaneous and equilibrium response
        """
        print('solving initial response', end = '...')
        exp = uc.experiments.ForceControlled(model, pars)
        sol_0 = exp.initial_response()
        print('done')

        # re-dim
        # sol_0.redimensionalise(pars)

        print('post-processing soln', end = '...')
        lam_r[n] = 1/np.sqrt(sol_0.lam_z)
        lam_z[n] = sol_0.lam_z

        _, _, Szz[n] = model.mechanics.eval_stress(lam_r[n], lam_r[n], lam_z[n])
        Pi[n] = model.osmosis.eval_osmotic_pressure(1)

        P[n] = sol_0.p[0]
        print('done')


    label = f'$\lambda_m = {lam_m[m]:.1f}$'

    F = -pars.physical["F"]
    mean_S = np.pi * Szz / F
    mean_P = np.pi * P * lam_r**2 / F
    mean_Pi = np.pi * Pi * lam_r**2 / F

    lw = 3; ms = 12
    ax_S.plot(Phi_f[1:], mean_P[1:] - mean_Pi[1:], lw = lw, label = label, ls = ls[m])

    if m == len(lam_m) - 1:
        label = f'No fibres'
        ax_S.plot(Phi_f[0], mean_P[0] - mean_Pi[0], 'd', ms = ms, markeredgecolor = 'k', label = label)


"""
Plot customisation
"""

fig_label = ['(b)']
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
# fig.savefig('peak_flf_2.pdf')

plt.show()