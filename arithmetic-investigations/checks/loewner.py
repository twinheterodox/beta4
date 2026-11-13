#!/usr/bin/env python3
"""One fixed coarse potential screen, with outward Arb checks.

No geometry or trial-measure optimization. See the manuscript and evidence guide for the uniform
determinant estimate and for the distinction between failure of these
bounds and impossibility of the construction.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
from flint import arb, ctx

ctx.prec = 160
ROOT = Path(__file__).resolve().parents[1] / "review"
C, MASS, LENGTH, HALF_WIDTH = Q(13,100), Q(213,200), Q(113,100), Q(1,5)
FIELD = Q(99,100)
POTENTIAL_LOWER = Q(2121,1250)  # 1.6968
TRIAL_ENERGY_UPPER = Q(2416291,1000000)
NODE_ENERGY_LOWER = Q(1614687,1000000)
DECAY = Q(311,125)  # 2.488


def ball(x):
    x = Q(x)
    return arb(x.numerator) / x.denominator


def xlog(x):
    # Test removable singularities on the exact rational argument.
    return arb(0) if x == 0 else ball(x) * ball(abs(x)).log()


def trial_potential(y):
    return 1 - (xlog(y+HALF_WIDTH)-xlog(y-HALF_WIDTH))/ball(2*HALF_WIDTH)


def log_integral(x, y):
    if y == 0:
        return xlog(x)-ball(x)
    return (ball(x)*ball(x*x+y*y).log()/2-ball(x)
            +ball(y)*ball(x/y).atan())


def node_potential(y):
    return -(log_integral(Q(1),y)+log_integral(LENGTH,y))/2


def primitive(x):
    return arb(0) if x == 0 else ball(x*x)*ball(abs(x)).log()/2-ball(3*x*x/4)


def interval_energy(a, b):
    return primitive(a-b)-primitive(a)-primitive(b)


def main():
    records = []
    for j in range(200):
        left, right = Q(j,100), Q(j+1,100)
        # Both potentials decrease on the positive imaginary axis.
        lower = (ball(MASS)*trial_potential(right)-node_potential(left)
                 +ball(FIELD*left)*arb.pi())
        assert lower > ball(POTENTIAL_LOWER), (j, lower)
        records.append({'left':str(left),'right':str(right),
                        'outward_lower_expression':str(lower)})
    # For y>=2: g(y)>=FIELD*pi*y-MASS*HALF_WIDTH/y, an increasing function.
    tail = ball(2*FIELD)*arb.pi()-ball(MASS*HALF_WIDTH/2)
    assert tail > ball(POTENTIAL_LOWER)
    trial_energy = ball(Q(3,2))-ball(2*HALF_WIDTH).log()
    node_energy = (interval_energy(Q(1),Q(1))+2*interval_energy(Q(1),LENGTH)
                   +interval_energy(LENGTH,LENGTH))/4
    assert trial_energy < ball(TRIAL_ENERGY_UPPER)
    assert node_energy > ball(NODE_ENERGY_LOWER)
    decay_from_rational_bounds = (2*MASS*POTENTIAL_LOWER
                                 -MASS*MASS*TRIAL_ENERGY_UPPER
                                 +NODE_ENERGY_LOWER)
    assert decay_from_rational_bounds > DECAY
    arithmetic = 16*(1+C)*(1+C/2)
    gap = DECAY-arithmetic
    assert arithmetic == Q(24069,1250) and gap == Q(-20959,1250)
    out = {
        'scope':'One geometry, one trial measure, 200 monotonic cells and one tail.',
        'status':'INCONCLUSIVE: compatible coarse estimates do not establish a positive gap.',
        'arb_precision_bits':ctx.prec,
        'parameters':{'c':str(C),'mass':str(MASS),'node_endpoint':str(LENGTH),
                      'trial_half_width':str(HALF_WIDTH),'field_fraction':str(FIELD)},
        'potential_lower':str(POTENTIAL_LOWER),
        'tail_lower_expression':str(tail),
        'trial_energy_ball':str(trial_energy),
        'trial_energy_upper':str(TRIAL_ENERGY_UPPER),
        'node_energy_ball':str(node_energy),
        'node_energy_lower':str(NODE_ENERGY_LOWER),
        'decay_from_rational_bounds':str(decay_from_rational_bounds),
        'certified_coarse_decay_rate':str(DECAY),
        'chosen_clearing_rate':str(arithmetic),
        'screened_bound_gap':str(gap),
        'decimal_rates':{'decay':float(DECAY),'chosen_clearing':float(arithmetic),
                         'bound_gap':float(gap)},
        'not_claimed':'The minimal denominator cost or optimal actual determinant decay; no impossibility theorem.',
        'cells':records,
    }
    (ROOT/'loewner-checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cells'},indent=2))


if __name__ == '__main__':
    main()
