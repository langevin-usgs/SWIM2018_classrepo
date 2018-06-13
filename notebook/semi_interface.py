import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipkinc as F
from scipy.special import ellipeinc as E
from scipy.optimize import fsolve, brentq


def fint(phi, a):
    kappasq = (2 + np.sqrt(3)) / 4  # Eq. 43
    kappa = np.sqrt(kappasq)
    # theta = arccos((sqrt(3) - 1 - phi / a) / (sqrt(3) + 1 + phi / a))  # Eq. 43
    theta = np.arccos(
        -1 + 2 * np.sqrt(3) / (np.sqrt(3) + 1 + phi / a))  # Eq. 43, reworked
    # Eq. 42
    rv = (3 ** (-0.25) - 3 ** 0.25) * F(theta, kappasq) + \
         2 * 3 ** 0.25 * E(theta, kappasq) - \
         2 * 3 ** 0.25 * np.sin(theta) * np.sqrt(1 - kappasq * np.sin(theta) ** 2) / (
         1 + np.cos(theta))
    return rv

class SemiCoastBase:
    def __init__(self, k, H, c, Q0, rhof, rhos, Ls, ztop=0, sealevel=0, label=None):
        self.k = k
        self.H = H
        self.c = c
        if self.c < 1e-12:
            self.c = 1e-12
        self.Q0 = Q0
        self.rhof = rhof
        self.rhos = rhos
        if Ls is None:
            self.Ls = np.inf
        else:
            self.Ls = Ls
        self.ztop = ztop
        self.hs = sealevel + (sealevel - ztop) * (self.rhos - self.rhof) / self.rhof
        self.label = label

    def findcase(self):
        # Variables
        self.lab = np.sqrt(self.k * self.H * self.c)
        self.nu = (self.rhos - self.rhof) / self.rhof
        self.mu = self.Q0 / (self.k * self.H) * self.lab / self.H / self.nu

        # Case I or Case II
        if self.mu < np.sqrt(2 / 3):
            self.Loverlab = (18 * self.mu) ** (1 / 3)  # Eq. 25
            if self.Loverlab * self.lab <= self.Ls:
                self.case = 1
                self.phi0 = (3 * self.mu ** 2 / 2) ** (1 / 3)
            else:
                self.case = None
        else:
            self.Loverlab = np.sqrt(6)  # Eq. 29
            self.doverlab = np.log((self.mu + np.sqrt(self.mu ** 2 + 1 / 3)) / \
                            (1 + np.sqrt(2 / 3)))  # Eq.38
            if self.Loverlab * self.lab + self.doverlab * self.lab <= self.Ls:
                self.case = 2
                self.phi0 = (1 - np.sqrt(2 / 3)) / 2 * np.exp(-self.doverlab) + \
                            (1 + np.sqrt(2 / 3)) / 2 * np.exp(self.doverlab)
            else:
                self.case = None

        if self.case is None:
            self.Lsoverlab = self.Ls / self.lab
            atr = fsolve(self.find_atr, 0.1)  # Eq. 46
            mutr = np.sqrt(2 * (1 + atr ** 3) / 3)  # Eq. 47
            if self.mu < mutr:
                self.case = 3
                amax = (3 * self.mu ** 2 / 2) ** (1 / 3)
                rv = fsolve(self.find_a_case3, amax / 2, full_output=1)
                if rv[2] == 1:
                    self.a = rv[0]
                    self.phi0 = (3 * self.mu ** 2 / 2 - self.a ** 3) ** (
                    1 / 3)  # Eq. 48
                else:
                    print('Error: a was not found for case 3, mu was:' + str(self.mu))
            else:
                self.case = 4
                rv = brentq(self.find_doverlab_case4, 0, self.Lsoverlab,
                            full_output=1)
                if rv[1].converged == 1:
                    self.doverlab = rv[0]
                    gamma0 = -np.tanh(self.doverlab) + self.mu / np.cosh(
                        self.doverlab)  # Eq.53
                    self.a = (3 * gamma0 ** 2 / 2 - 1) ** (1 / 3)  # Eq.54
                    self.phi0 = 1.0 / np.cosh(self.doverlab) - self.mu * \
                                np.sinh(-self.doverlab) / np.cosh(self.doverlab)
                else:
                    print('Error: doverlab was not found for case 4')
    
    def find_atr(self, a):
        # Eq. 46. atr is the zero of this function
        rv = np.sqrt(3 * a / 2) * (fint(1, a) - fint(0, a)) + self.Lsoverlab
        return rv

    def find_a_case3(self, a):
        phi0 = (3 * self.mu ** 2 / 2 - a ** 3) ** (1 / 3)
        rv = np.sqrt(3 * a / 2) * (fint(phi0, a) - fint(0, a)) + self.Lsoverlab
        return rv

    def find_doverlab_case4(self, doverlab):
        # Function that has to become zero for case IV
        gamma0 = -np.tanh(doverlab) + self.mu / np.cosh(doverlab)  # Eq.53
        term = 3 * gamma0 ** 2 / 2 - 1  # Eq.54 (argument only)
        # print('doverlab, term', doverlab, term)
        if term < 0:
            # print('term < 0 in find_doverlab_case4')
            rv = -1
        else:
            a = term ** (1 / 3)  # Eq.54
            rv = np.sqrt(3 * a / 2) * (
            fint(1, a) - fint(0, a)) + self.Lsoverlab - doverlab
        return rv
    
class Coast:
    def __init__(self, k, H, grad, rhof, rhos, ztop=0, label=None):
        self.k = k
        self.H = H
        self.grad = np.abs(grad)
        self.ztop = ztop
        self.label = label
        self.rhof = rhof
        self.rhos = rhos
        self.Q0 = self.k * self.H * self.grad
        self.D = np.abs(self.ztop)
        self.alpha = self.rhof / (self.rhos - self.rhof)
    def head(self, x):
        x = np.atleast_1d(x)
        phitoe = 0.5 * self.k * self.H ** 2 / self.alpha
        C = -self.k * (self.D * self.H + 0.5 * self.H ** 2) / self.alpha
        phi = -self.Q0 * x
        h = np.zeros_like(phi)
        h[phi <= phitoe] = np.sqrt(2 * phi[phi <= phitoe] / (self.k * self.alpha)) + self.D / self.alpha
        h[phi > phitoe] = (phi[phi > phitoe] - C) / (self.k * self.H)
        return h
    def interface(self, x):
        h = self.head(x)
        zeta = -self.alpha * h
        zeta[zeta < -(self.D + self.H)] = -(self.D + self.H)
        return zeta
    def toe(self):
        phitoe = 0.5 * self.k * self.H ** 2 / self.alpha
        xtoe = -phitoe / self.Q0
        return xtoe
    def plot(self, xmin=None, newfig=True):
        if newfig:
            if xmin is None:
                xmin = self.toe()
            plt.figure(figsize=(8, 4))
            plt.plot([xmin, 0], [self.ztop - self.H, self.ztop - self.H], linewidth=5, color='k')
            plt.plot([xmin, 0], [self.ztop, self.ztop], linewidth=5, color='k')
            plt.xlim(xmin, 0)
            plt.ylim(self.ztop - self.H, self.ztop)
        x = np.linspace(self.toe(), 0, 100)
        z = self.interface(x)
        plt.plot(x, z, zorder=100, label=self.label)
        return plt.gca()

class SemiCoast(SemiCoastBase):
    def __init__(self, k, H, c, grad, rhof, rhos, Ls=None, ztop=0, sealevel=0, label=None):
        SemiCoastBase.__init__(self, k, H, c, grad * k * H, rhof, rhos, Ls, ztop, sealevel, label)
        self.findcase()
        self.initialize()
    def initialize(self):
        if self.case == 1:
            self.doverlab = (1 - (3 * self.mu ** 2 / 2) ** (2 / 3)) / (
            2 * self.mu)  # Eq.28
            self.xtoe = -self.doverlab * self.lab
            self.xtip = self.Loverlab * self.lab
        elif self.case == 2:
            self.xtoe = self.doverlab * self.lab
            self.xtip = (self.doverlab + self.Loverlab) * self.lab
        elif self.case == 3:
            self.doverlab = (1 - self.phi0 ** 2) / (2 * self.mu)  # Eq.27 solved for u with phi=1 and phi0
            self.xtoe = -self.doverlab * self.lab
            self.xtip = self.Lsoverlab * self.lab
        elif self.case == 4:
            self.xtoe = self.doverlab * self.lab
            self.xtip = self.Lsoverlab * self.lab
        else:
            print('Error: case was not found')

    def toe(self):
        return self.xtoe

    def tip(self):
        return self.xtip
    
    def interface(self, N=100):
        # returns (x, z) where z is interface elevation for N points
        if self.case == 1:
            phi = np.nan * np.ones(N)
            u = np.linspace(-self.doverlab, self.Loverlab, N)
            u1 = (-self.doverlab <= u) & (u <= 0)
            phi[u1] = np.sqrt(-2 * self.mu * u[u1] + self.phi0 ** 2)  # Eq. 27
            u2 = (0 < u) & (u <= self.Loverlab)
            phi[u2] = (u[u2] - self.Loverlab) ** 2 / 6
            x = u * self.lab
        elif self.case == 2:
            u = np.linspace(0, self.Loverlab, N)
            phi = (u - self.Loverlab) ** 2 / 6
            x = (self.doverlab + u) * self.lab
        elif self.case == 3:
            phi = np.linspace(1, 0, N)
            u = np.nan * np.ones(N)
            phi1 = phi >= self.phi0
            u[phi1] = (self.phi0 ** 2 - phi[phi1] ** 2) / (2 * self.mu)
            phi2 = phi < self.phi0
            u[phi2] = np.sqrt(3 * self.a / 2) * (
            fint(phi[phi2], self.a) - fint(0, self.a)) + self.Lsoverlab
            x = u * self.lab
        elif self.case == 4:
            phi = np.linspace(1, 0, N)
            u = np.sqrt(3 * self.a / 2) * (fint(phi, self.a) - fint(0, self.a)) + (
            self.Lsoverlab - self.doverlab)
            x = (self.doverlab + u) * self.lab
        h = self.nu * self.H * phi + self.hs
        eta = (h - self.hs) / self.nu
        zeta = self.ztop - eta
        return x, zeta

    #def head(self, x):
    #    # Only implemented for Case 1
    #    if self.case == 1:
    #        u = x / self.lab
    #        if u <= -self.doverlab:
    #            phi = self.grad * (-self.doverlab - u) + 1
    #        elif (u > -self.doverlab) & (u <= 0):
    #            phi = sqrt(-2 * self.mu * u + self.phi0 ** 2)
    #        elif (u > 0) & (u <= self.Loverlab):
    #            phi = (u - self.Loverlab) ** 2 / 6
    #        else:
    #            phi = nan
    #    elif self.case == 2:
    #        pass
    #    return self.nu * self.H * phi + self.hs

    def onshorex(self, h):
        # returns onshore x location of given head
        phi = (h - self.hs) / (self.nu * self.H)
        if self.case == 1:
            if phi >= 1:
                u = (1 - phi) * (self.nu * self.H) / (
                self.grad * self.lab) - self.doverlab
            elif (phi < 1) & (phi > self.phi0):
                u = (self.phi0 ** 2 - phi ** 2) / (2 * self.mu)
            else:
                u = 0
            x = u * self.lab
        elif self.case == 2:
            if phi > self.phi0:
                u = (self.phi0 - phi) * (self.nu * self.H) / (
                    self.grad * self.lab) - self.doverlab
            else:
                u = -doverlab
            x = (self.doverlab + u) * self.lab
        elif self.case == 3:
            if phi >= 1:
                u = (1 - phi) * (self.nu * self.H) / (
                self.grad * self.lab) - self.doverlab
            elif (phi < 1) & (phi > self.phi0):
                u = (self.phi0 ** 2 - phi ** 2) / (2 * self.mu)
            else:
                u = 0
            x = u * self.lab
        elif self.case == 4:
            if phi > self.phi0:
                u = (self.phi0 - phi) * (self.nu * self.H) / (
                    self.grad * self.lab) - self.doverlab
            else:
                u = -self.doverlab
            x = (self.doverlab + u) * self.lab            
        return x
    
    def plot(self, xmin=None, xmax=None, newfig=True):
        if newfig:
            if xmin is None:
                xmin = self.toe()
            if xmax is None:
                xmax = self.tip()
            plt.figure(figsize=(8, 4))
            plt.plot([xmin, xmax], [self.ztop - self.H, self.ztop - self.H], linewidth=5, color='k')
            plt.plot([min(xmin, 0), min(xmax, 0)], [self.ztop, self.ztop], linewidth=5, color='k')
            plt.plot([max(xmin, 0), min(xmax, self.Ls)], [self.ztop, self.ztop], linewidth=5, color=[.8, .8, .8])
            plt.xlim(xmin, xmax)
            plt.ylim(self.ztop - self.H, self.ztop)
        x, z = self.interface()
        plt.plot(x, z, zorder=100, label=self.label)
        return plt.gca()
    
class SemiCoastHead(SemiCoast):
    def __init__(self, k, H, c, h, x, rhof, rhos, Ls=None, ztop=0, sealevel=0, label=None):
        assert x <=0, "Input error: x must be less than zero"
        if Ls is None:
            self.Linput = np.inf
        else:
            self.Linput = Ls
        # First find position without finite L, then use the specified L
        # This is a bit clunky, but may work ok
        SemiCoast.__init__(self, k, H, c, h / np.abs(x), rhof, rhos, np.inf, ztop, sealevel, label)
        assert h > self.hs, "Input error: inland head smaller than equivalent freshwater head at top of aquifer"
        self.givenx = x
        self.givenh = h
        # find gradient using log transform of grad to avoid negative values
        # Start is exact solution for regular confined flow dived by 10 to start out low and avoid errors
        start = (self.givenh - self.hs) / (np.abs(x) + self.lab) / 10  
        #result = brentq(self.findgrad2, 1e-12, start)
        result = fsolve(self.findgrad, np.log(start), full_output=1)
        if result[2] == 1:
            self.grad = np.exp(result[0])
            self.initialize()
        else:
            print('Error: gradient could not be found with fsolve when L=inf')
        #
        if self.tip() > self.Linput:
            self.Ls = self.Linput
            self.initialize()
            result = fsolve(self.findgrad, np.log(self.grad), full_output=1)
            if result[2] == 1:
                self.grad = np.exp(result[0])
                self.initialize()
            else:
                print('Error: gradient could not be found with fsolve when L=inf')
            
    def findgrad(self, g):
        # g is log transformed
        self.grad = np.exp(g)
        self.Q0 = self.grad * self.k * self.H
        self.findcase()
        self.initialize()
        xh = self.onshorex(self.givenh)
        return xh - self.givenx
    
class IslandInterface:
    def __init__(self, k=10, D=50, c=100, rhof=1000, rhos=1025, W=1000, N=0.001, Nstreamlines=None):
        self.k = k
        self.D = D
        self.c = np.max((c, 1e-12))  # make sure c is above zero
        self.rhof = rhof
        self.rhos = rhos
        self.W = W
        self.N = N
        self.Nstreamlines = Nstreamlines
        #
        self.alpha = self.rhof / (self.rhos - self.rhof)
        self.phitoe = 0.5 * self.k * (self.alpha + 1) / (self.alpha ** 2) * self.D ** 2
        self.C = -0.5 * self.k * (self.alpha + 1) / self.alpha * self.D ** 2
        #
        self.lab = np.sqrt(self.k * self.D * self.c)
        self.nu = (self.rhos - self.rhof) / self.rhof
        self.grad = self.N * self.W / (self.k * self.D)
        self.mu =  self.grad * self.lab / self.D / self.nu
        self.h0 = self.nu * self.D * (3 * self.mu ** 2 / 2) ** (1 / 3)
        self.Loutflow = (18 * self.mu) ** (1 / 3) * self.lab
    def plot(self, xmin=None, xmax=None, plotbase=True):
        x = np.linspace(-self.W, self.W, 100)
        phicoast = 0.5 * self.k * (self.alpha + 1) * self.h0 ** 2
        phi = -self.N / 2 * (x ** 2 - self.W ** 2) + phicoast
        h =  np.nan * np.ones(len(x))
        h[phi <= self.phitoe] = np.sqrt(2 * phi[phi <= self.phitoe] / (self.k * (self.alpha + 1)))
        h[phi > self.phitoe] = np.sqrt(2 / self.k * (phi[phi > self.phitoe] - self.C)) - self.D
        zeta = -self.D * np.ones(len(x))
        zeta[phi <= self.phitoe] = -self.alpha * h[phi <= self.phitoe]
        # Outflow zone
        x2 = np.linspace(0, self.Loutflow, 100)
        phi2 = (x2 / self.lab - self.Loutflow / self.lab) ** 2 / 6
        h2 = self.nu * self.D * phi2
        zeta2 = -h2 / self.nu
        # Plot results
        plt.figure(figsize=(10, 5))
        plt.plot(x, h, 'b')
        plt.plot(x2 + self.W, h2, 'b')
        plt.plot(-x2 - self.W, h2, 'b')
        plt.plot(x, zeta, 'r')
        plt.plot(x2 + self.W, zeta2, 'r')
        plt.plot(-x2 - self.W, zeta2, 'r')
        # Plot streamlines
        xh = np.hstack(((-x2 - self.W)[::-1], x, x2 + self.W))
        Qoutflow = -self.k * (0 - zeta2) * self.D * self.nu * (x2 / self.lab - self.Loutflow / self.lab) / (3 * self.lab)
        Qh = np.hstack((-Qoutflow[::-1], self.N * x, Qoutflow))
        htop = np.hstack((h2[::-1], h, h2))
        zetabot = np.hstack((zeta2[::-1], zeta, zeta2))
        xs = np.vstack((xh, xh))
        ys = np.vstack((htop, zetabot))
        Qs = np.vstack((Qh, np.zeros(len(xh))))
        if self.Nstreamlines is not None and self.Nstreamlines > 0:
            plt.contour(xs, ys, Qs, self.Nstreamlines)
        if xmin is None:
            xmin = -self.W - self.Loutflow
        if xmax is None:
            xmax = self.W + self.Loutflow
        plt.plot([min(xmin, -self.W - self.Loutflow), -self.W], [0, 0], linewidth=5, color=[.8, .8, .8])
        plt.plot([min(xmax, self.W + self.Loutflow), self.W], [0, 0], linewidth=5, color=[.8, .8, .8])
        if plotbase: plt.plot([xmin, xmax], [-self.D, -self.D], linewidth=5, color='k')
        plt.xlim(xmin, xmax)
        return plt.gca()


## Input
## Case 1
#sc1 = SemiCoast(k=10, H=10, c=100, grad=0.0005, rhof=1000, rhos=1025, Ls=1000)
## Case 2
#sc2 = SemiCoast(k=10, H=10, c=100, grad=0.00375, rhof=1000, rhos=1025, Ls=1000)
## Case 3
#sc3 = SemiCoast(k=10, H=10, c=100, grad=0.0005, rhof=1000, rhos=1025, Ls=80)
## Case 4
#sc4 = SemiCoast(k=10, H=10, c=100, grad=0.00375, rhof=1000, rhos=1025, Ls=150)
#
#sch = SemiCoastHead(k=10, H=10, c=100, x=-1000, h=0.25, rhof=1000, rhos=1025, Ls=inf, ztop=0, sealevel=0)

    
