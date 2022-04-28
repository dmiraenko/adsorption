from pint import UnitRegistry
from functions import *

if __name__ == "__main__":
    ureg = UnitRegistry()
    total_pressure = 200 * ureg.kPa
    y = [0.25, 0.25, 0.25, 0.25] * ureg.dimensionless
    cmus = [32.47, 82.34, 584.43, 30.15] * ureg.mmol / ureg.gram
    b = [0.255, 2.7682, 97.7962, 23.703] * (1 / ureg.MPa)
    t = [0.777, 0.323, 0.134, 0.660] * ureg.dimensionless

    funcs = [Toth(cmus[i], b[i], t[i]) for i in range(len(y))]

    z = newton(lambda z: total_pressure * sum([y[i] / funcs[i].pure_comp_pressure(z) for i in range(len(y))]) - 1 ,\
        sum([y[i] * funcs[i].red_spr_pressure(total_pressure) for i in range(len(y))]), \
        fprime = lambda z : - total_pressure * sum([y[i] / ((p0 := funcs[i].pure_comp_pressure(z)) * funcs[i].surface_conc(p0)) for i in range(len(y))])
    )

    print([(total_pressure * y[i] / funcs[i].pure_comp_pressure(z)).to('').magnitude for i in range(len(y))])
