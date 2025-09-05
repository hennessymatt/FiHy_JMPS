"""
Plots the stress-strain data in Fig 5
"""

import ucompress as uc 
import matplotlib.pyplot as plt
import matplotlib.colors as colours
import numpy as np
from fihy_params import FiHy25, FiHy66, FiHy71, FiHy75

cols = colours.TABLEAU_COLORS
plt.rcParams.update({
    "text.usetex": True,
    # "font.family": "sans-serif",
    # "font.sans-serif": ["Helvetica"],
    "font.size": 20,
    "figure.autolayout": True
    })


"""
Create structures for storing experimental data
"""
data_dir = 'exp_data/'

data_25 = {
    'files': {
        'stress_strain_fihy025_a3.csv': 'FiHy025 a3',
        'stress_strain_fihy025_a4.csv': 'FiHy025 a4'
    },
    'params': FiHy25()
}

data_66 = {
    'files': {
        'stress_strain_fihy066_c5.csv': 'FiHy066 c5',
        'stress_strain_fihy066_b4.csv': 'FiHy066 b4'
    },
    'params': FiHy66()
}

data_71 = {
    'files': {
        'stress_strain_fihy071_b5.csv': 'FiHy071 b5',
        'stress_strain_fihy071_c4.csv': 'FiHy071 c4'
        },
    'params': FiHy71()
}

data_75 = {
    'files': {
        'stress_strain_fihy075_a4.csv': 'FiHy075 a4',
        'stress_strain_fihy075_b5.csv': 'FiHy075 b5'
        },
    'params': FiHy75()
}


data = {
     25: data_25,
     66: data_66,
     71: data_71,
     75: data_75
}

"""
Set up the modelling parameters
"""

# set up the model
pars = FiHy71()
mech = uc.mechanics.FibreRecruitment(distribution='quartic', homogeneous=True)
os = uc.osmosis.FloryHuggins()
perm = uc.permeability.Constant()

model = uc.base_models.Hydrogel(mech, perm, os, pars)

mech_wo = uc.mechanics.FibreReinforced(homogeneous=True)
model_wo = uc.base_models.Hydrogel(mech_wo, perm, os, pars)


"""
Run the comparison
"""
max_strain = [0.5, 0.3, 0.3, 0.5]
ylim = [0.5, 0.5, 0.5, 0.5]

for d, eps_max, yl in zip(data, max_strain, ylim):
    print('==================================')
    plt.figure()

    """
    Set the parameters for the fresh material
    """
    pars = data[d]["params"]
    print(f'Phi_f = {pars.physical["Phi_f"]}')

    J_0 = 1 / (1 - pars.physical["phi_0"])

    model.assign(pars)
    model_wo.assign(pars)

    """
    Solve the hydration problem for a fixed value of J_0
    """
    for m in [model, model_wo]:
        chi_calc = uc.fitting.ChiCalculator(m, pars)
        chi, beta_r, beta_z, phi_0 = chi_calc.solve(J_0)

        pars.update("beta_r", beta_r)
        pars.update("beta_z", beta_z)
        pars.update("phi_0", phi_0)
        m.assign(pars)

    """
    Instant compression of the material
    """
    # axial strain and stretches
    eps = np.linspace(0, eps_max, 30)
    lam_z = 1 - eps

    lam_r = 1 / np.sqrt(lam_z)

    for m, ls in zip([model, model_wo], ['-', '--']):
        S_r, _, S_z = m.mechanics.eval_stress(lam_r, lam_r, lam_z)
        p = lam_r * S_r

        # Calculate the total PK1 stress
        S_z_T = S_z - lam_r**2 * p

        # Plot
        plt.plot(eps, 
                -S_z_T / 1e6, 
                ls = ls,
                lw = 4, 
                color = 'k'
        )
        
    """
    Plot the exp data
    """
    marker = ['o', 's']
    n = 0
    for f in data[d]['files']:

        D = np.loadtxt(f'{data_dir}{f}', delimiter=',', skiprows=1)
        eps_data = D[:,0]
        sigma = D[:,1]

        plt.plot(eps_data, 
                sigma, 
                marker[n], 
                ms = 10, 
                label = data[d]['files'][f], 
                markeredgecolor = 'k',
                mew = 1
                )
        n += 1

    plt.xlabel('Axial strain (mm/mm)')
    plt.ylabel('Axial stress (MPa)')
    plt.xlim([0, 0.5])
    plt.ylim([0, yl])
    plt.legend()
    plt.savefig(f'stress_strain_{d}.pdf')

plt.show()