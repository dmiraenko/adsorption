from numpy import sqrt
from pint import UnitRegistry
from functions import *

class IAS():

    def __init__(self, funcs):
        mod = __import__("functions")
        # for f in funcs:
        #     print(f["kw"])
        self.funcs = [getattr(mod, f["fname"])(**(f["kw"])) for f in funcs]

    def simple_IAS(self, total_pressure, y):
        z = newton(lambda z: total_pressure * sum([y[i] / self.funcs[i].pure_comp_pressure(z) for i in range(len(y))]) - 1 ,\
            sum([y[i] * self.funcs[i].red_spr_pressure(total_pressure) for i in range(len(y))]), \
            fprime = lambda z : - total_pressure * sum([y[i] / ((p0 := self.funcs[i].pure_comp_pressure(z)) * self.funcs[i].surface_conc(p0)) for i in range(len(y))])
        )

        return [(total_pressure * y[i] / self.funcs[i].pure_comp_pressure(z)).to('').magnitude for i in range(len(y))]

    def fast_IAS(self, total_pressure, y):
        ct = total_pressure * sum([self.funcs[i].henry_law_const * y[i] for i in range(len(y))])
        p0 = [0.01 * ct / f.henry_law_const for f in self.funcs]
        print(p0)

        it = 0
        eps = 1e-7
        max_iter = 1000
        while(True):
            z0 = self.funcs[0].red_spr_pressure(p0[0])
            g = [sum([y[i] / p0[i] for i in range(len(y))]) - 1 / total_pressure] + \
                [z0 - self.funcs[i].red_spr_pressure(p0[i]) for i in range(1, len(y))]

            j = [self.funcs[0].surface_conc(p0[0]) / p0[0]] + \
                [y[i] / (self.funcs[i].surface_conc(p0[i]) * p0[i]) for i in range(1, len(y))]

            fsum = 0
            for i in range(1, len(y)):
                fsum += j[i]
                g[0] -= j[i] * g[i]
            
            g[0] /= - y[0] / p0[0]**2 - j[0] * fsum
            g[1:] = [(j[0] * g[0] - g[i]) * j[i] * p0[i]**2 / y[i] for i in range(1, len(y))]

            p0 = [k if (k := p0[i] - g[i]).magnitude > 0 else eps * total_pressure.units for i in range(len(y))]
            print(p0)

            e = sqrt(sum([((g[i] / p0[i]).to("").magnitude)**2 for i in range(len(y))]))
            it += 1
            if (e < eps):
                return p0
            if (it > max_iter):
                raise TimeoutError("Number of iterations exceeded the maximum.")
            

if __name__ == "__main__":

    ureg = UnitRegistry()
    total_pressure = 700 * ureg.kPa
    y = [0.25, 0.25, 0.25, 0.25] * ureg.dimensionless
    cmus = [32.47, 82.34, 584.43, 30.15] * ureg.mmol / ureg.gram
    b = [0.255, 2.7682, 97.7962, 23.703] * (1 / ureg.MPa)
    t = [0.777, 0.323, 0.134, 0.660] * ureg.dimensionless

    dat = [{"fname" : "Toth", "kw" : {"sat_surface_conc" : cmus[i], "b" : b[i], "t" : t[i]}} for i in range(len(y))]

    ias = IAS(dat)

    print(ias.fast_IAS(total_pressure, y))



