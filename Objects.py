import numpy as np

class Boat:

    def __init__(self, sf=1.3, mass=4324, b=1.84, beta=13, loa=7.9,
                 rho=1000, nu=1.19e-6, f=0, eps=3, a=np.sqrt(2), r=0.5):

        self.sf = sf  # Factor de seguridad
        self.mass = mass  # Masa del bote en Kg
        self.b = b  # Manga del bote en m
        self.beta = beta  # Ángulo muerto del centro del barco en °
        self.loa = loa  # Slora
        self.l_cg = loa * 0.4  # centro de gravedad longitudinal
        self.rho = rho  # Densidad del agua en kg/m3
        self.nu = nu  # Viscosidad cinematica del agua
        self.f = f  # Distancia perpendicular desde cg hasta la linea thrust
        self.eps = eps  # Ángulo del vector de empuje medido desde la quilla (epsilon)
        self.a = a  # Factor de interferencia, para un catamaran raiz de 2
        self.r = r  # Tasa de separación

    def drag_sav(self, v, mass, eta_T=0.95, eta_D=0.668):
        from scipy.optimize import fsolve
        """
        Function for calculating drag, power and boat position according to Savistsky
        method.
        Savistsky originally calculates the effective power.

        Source of model: Daniel Savistsky, 1964,
            'Hydrodynamic design of planing hulls'.
        Source with an easier explanation: Molland,
            'Ship resistance and propulsion', page 216.
        Source for catamarans adaptation: 1979,
            'Interference effect of catamaran planing hulls'.

        Input arguments:
            mass: Boat mass (kg)
            b: Mean chine beam (m)
            beta: Deadrise angle (ì§¸)
            l_cg: Longitudinal center of gravity (m)
            rho: water density (kg/m3)
            nu: Kinematic viscosity of water (-)
            v: Boat speed (m/s)
            f: Perpendicular distance from thrust line to CG (m)
            eps: Thrust line angle with respect to keel (°)
                        if f and eps are zero, thrust line passes through cg and parallel to keel
            A: Interference factor due to demihulls [-]. 1 for a monohull
            r: Separation ratio due to demihulls [-]. 1 for a monohull
            eta_T: Transmission efficiency [0-1] (-)
            eta_D: Quasi-propulsive coefficient [0-1] (-)

        Output arguments:
            tao: Angle of trim (ì§¸)
            d: Draft of the keel at transom (m)
            T: Trust. Calculated with the effective power. (kN)
            inst_kw: Installer power (kW)
            R_t: Total resistance (kN)
            S: Wetted area (m2)
            lk: Wetted keel length (m)
            lm: Mean keel wetted length (m)
            lp: Distance from keel-transom to longitudinal center of pressure (m)

        """

        c_v = v / np.sqrt(9.81 * self.b)  # Coefficient of viscous resistance.
        c_lbeta = (self.mass + mass) * 9.81 / (0.5 * self.rho * self.b ** 2 * v ** 2)  # Lift coefficient of a
        # warped surface.
        Eq1 = lambda c_l0: c_lbeta - c_l0 + 0.0065 * self.beta * c_l0 ** 0.6  # equation for
        # calculating c_l0.
        c_l0_initial_guess = 0.5
        c_l0 = fsolve(Eq1, c_l0_initial_guess)[0]  # Lift coefficient o a flat surface

        # Iterative calculation
        tao = 0
        delta_M = -1

        while delta_M <= 0.0:

            tao += 1e-2
            Eq2 = lambda lamb: c_l0 - (tao ** 1.1 * self.r ** (3 / 2)) * ((0.012 * lamb ** 0.5) / self.a \
                                                                     + 0.0055 * lamb ** (2.5) * self.a / (
                                                                                 (c_v ** 2) * self.r))  # Equation for
            # calculating lambda.
            lamb = fsolve(Eq2, 1)[0]
            lm = lamb * self.b
            lp = lm * (0.75 - 1 / (5.21 * ((c_v / lamb) ** 2) * (self.r / self.a ** 2) + 2.39))

            re = v * lm / self.nu
            cf = 0.075 / (np.log10(re) - 2) ** 2  # coefficient of frictional resistance

            if self.a == 1 and self.r == 1:  # If the boat is a monohull, then...
                S = lm * self.b / (np.cos(np.deg2rad(self.beta)))  # wetted surface area for a monohull
            else:  # If these values are different from 1, then it is a catamaran
                S = self.b ** 2 * ((self.r / np.cos(np.deg2rad(self.beta))) + 2 * lamb * np.tan(np.deg2rad(tao)))

            Df = 0.5 * self.rho * S * v ** 2 * cf  # drag due to frictional resistance

            T = ((self.mass + mass) * 9.81 * np.sin(np.deg2rad(tao)) + Df) / 1000 * np.cos(np.deg2rad(self.eps))  # thrust [kN]
            N = (self.mass + mass) * 9.81 * np.cos(np.deg2rad(tao)) / 1000 - T * np.sin(np.deg2rad(self.eps))

            PT_L = self.l_cg * np.sin(np.deg2rad(self.eps)) + self.f  # perpendicular vertical
            # distance from point p from where moments
            # are calculated to thrust line

            delta_M = -N * lp + (self.mass + mass) * 9.81 * self.l_cg / 1000 - T * PT_L * (np.cos(np.deg2rad(self.eps)) +
                                                                                  np.sin(np.deg2rad(self.eps)) *
                                                                                  np.cos(np.deg2rad(90 - self.eps)))

        lk = (1 / 2) * ((self.b / np.pi) * (np.tan(np.deg2rad(self.beta)) / np.tan(np.deg2rad(tao))) \
                        + lm)  # Keel wetted length.
        # lk results from combining equations 3 and 5 from "hydrodynamic
        # design of planing hulls" by savistky.
        d = lk * np.sin(np.deg2rad(tao))  # Draft of the keel at transom
        # kCG = VCG*np.cos(tao) + np.sin(tao)*( l_cg-v_cg*np.tan(tao) )# [m] Vertical
        # distance from keel-transom intersection to CG. Equation found
        # trough geometric calculations.☻

        R_t = T * np.cos(np.deg2rad(tao + self.eps))  # total drag [kN]
        P_e = R_t * v  # effective power [kW]

        inst_kw = P_e * (1 / eta_D) * (1 / eta_T)  # installed power [kW]

        # print(lm, lk, lamb)

        return inst_kw
