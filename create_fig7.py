"""
Solves the time-dependent unconfined compression problem and creates
Fig 7
"""

import matplotlib.pyplot as plt
import ucompress as uc
import numpy as np
import matplotlib.colors as colours
from fihy_params import FiHy71

cols = colours.TABLEAU_COLORS

plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    "font.size": 20,
    "figure.autolayout": True
    })


"""
Set the parameters
"""
# psi = 0.8
psi = 0.96
forces = [-1, -5, -10]
ls = ['-', '--', '-.']

"""
Set up the model
"""
pars = FiHy71(nondim=True)
pars.update('phi_0', 0.97)


mech = uc.mechanics.FibreRecruitment(distribution='quartic')
perm = uc.permeability.HolmesMow()
os = uc.osmosis.FloryHuggins()

model = uc.base_models.Hydrogel(mech, perm, os, pars)

"""
Solve the hydration problem
"""
J_h = 1 / (1 - pars.physical["phi_0"])
print(f'J_h = {J_h:.2f}')

chi_calc = uc.fitting.ChiCalculator(model, pars)
chi, beta_r, beta_z, phi = chi_calc.solve(J_0 = J_h)

pars.update("chi", chi)
pars.update("beta_r", beta_r)
pars.update("beta_z", beta_z)
pars.update("phi_0", phi)

"""
Solve the time-dependent problem
"""
n = 0
for f, c in zip(forces, cols):

    pars.update('t_end', 1e5)
    pars.update('F', f)
    model.assign(pars)

    # solve the model and plot
    problem = uc.experiments.ForceControlled(model, pars)
    problem.solver_opts["monitor_convergence"] = False
    problem.solver_opts["newton_tol"] = 1e-5
    problem.solver_opts["newton_max_iterations"] = 20

    sol = problem.transient_response()
    sol.redimensionalise(pars)

    label = f'$F = {f}$ N'
    plt.semilogx(sol.t, 
                 sol.lam_z, 
                 lw = 3, 
                 color = cols[c],
                 label = label,
                 ls = ls[n]
                 )
    n += 1

    # load the exp data
    D = np.loadtxt(f'exp_data/Fig4_v2-FiHy071-load_{-f}.0N.csv', delimiter = ',')
    t = D[:,0]
    eps = D[:,2]
    lam = 1 + eps

    # interpolate the data
    t_i = np.logspace(-1, 5, 30)
    lam_i = np.interp(t_i, t, lam)

    # plot the exp data
    plt.semilogx(t_i, 
                 lam_i, 
                 'o', 
                 ms = 10,
                 color = cols[c],
                 mew = 1,
                 markeredgecolor = 'k'
    )
    


plt.xlabel(r'$t$ (s)')
plt.ylabel(r'$\beta_z$ (-)')
plt.xticks([1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5])
plt.legend()

# plt.savefig(f'transient_{int(J_h)}.pdf')
plt.show()
