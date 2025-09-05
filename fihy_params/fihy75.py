from ucompress.parameters.base_parameters import Parameters

class FiHy75(Parameters):
    """
    A parameter class for FiHy75
    """

    def __init__(self, nondim = False):
        
        super().__init__(nondim = nondim)
        
        # Label
        self.label = 'FiHy75'

    def set_parameters(self):

        """
        Define the dimensional parameters.  Here, we use SI units
        """

        R_g = 8.314
        V_m = 18e-6
        T = 23 + 273
        G_T = R_g * T / V_m

        self.physical = {
            "R": 5.165e-3,        # initial radius (m)
            "E_m": 4000,      # stiffness of gel matrix (Pa)
            "nu_m": 0,        # Poisson's ratio of the material
            "E_f": 4e8,         
            "lam_m": 2.9713,      # fitted to stress-strain w/ quartic dist
            "Psi_f": 1-0.75,
            "Phi_f": 0.26,
            "lam_z": 0.4,     # axial strain (-)
            "F": -10,         # applied force (N)
            "chi": 0.5657,    
            "G_T":   G_T,

            # Densities [g/ml]
            "rho_w": 1,
            "rho_m": 1.2,
            "rho_f": 1.145,

            # Permeability parameters            
            "k_0": 2e-13,     # initial hydraulic conductivity (m2 / Pa / s)
            "M": 2,
            "a": 4,

            # Initial state of the material 
            "psi_0": 0.8,   # initial mass fraction of solvent (-)  
            "phi_0": 0.8259,     # initial porosity (-)
            "beta_r": 1,
            "beta_z": 1,

            # Experiment start/end times
            "t_start": 1e-2,    # start time (s)
            "t_end": 1e3,       # final time (s)


        }

        """
        Define the computational parameters
        """
        self.computational = {
            "N": 50,
            "Nt": 200,
            "t_spacing": 'log'
        }

    

    def non_dimensionalise(self):
        """
        Carries out the non-dimensionalisation of all of the physical
        parameters as well as the final simulation time.  In its 
        current form, this method only applies to simple neo-Hookean
        materials.  This method will therefore have to be overloaded
        if using a more complex model with extra parameters.
        """

        # copy the dimensional dict into the physical dict
        self.physical = self.dimensional.copy()

        # non-dimensionalise all dimensional quantities
        self.physical["R"] /= self.scaling["space"]
        self.physical["E_m"] /= self.scaling["stress"]
        self.physical["E_f"] /= self.scaling["stress"]
        self.physical["G_T"] /= self.scaling["stress"]
        self.physical["k_0"] /= self.scaling["permeability"]
        self.physical["F"] /= self.scaling["force"]
        self.physical["t_start"] /= self.scaling["time"]
        self.physical["t_end"] /= self.scaling["time"]