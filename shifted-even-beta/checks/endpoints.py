"""Exact rational decay certificates for the shifted-interval examples."""
from checks import Params, F, decay
from pathlib import Path
import json


def main():
    assert __debug__, 'Run without -O'
    cases = [
        # a means that the first index is 2a; m=2a-1.
        (2, 1, 198, 7, F(21, 2000)),
        (3, 1, 330, 7, F(103, 10000)),
        (4, 1, 462, 7, F(104, 10000)),
        (10, 2, 3080, 7, F(19, 10000)),
        (10, 3, 6950, 7, F(59, 100000)),
    ]
    result = [decay(Params(B, a, s+1, r), u) for a, B, s, r, u in cases]
    # A single rational comparison also certifies every odd m for B=1:
    # S=67m, U=U0^m e0^(4(m-1)), and U0*e0^4<1.
    e0, u = F(68, 25), F(1, 50)
    f0 = ((u+18)/u)**4*((u+7)/(u+11))**67
    U0 = e0**(4*67)*(18+u)**72*4**(4*67)*(F(7)**7/F(11)**11)**67
    assert f0 < 1 and U0*e0**4 < 1
    record = {'status': 'All displayed inequalities checked over Q; decimals illustrative.',
              'examples': result,
              'B1_all_odd_m': {'S': '67m', 'u': '1/50', 'f0_lt_1': True,
                              'U0_times_e0_power4_lt_1': True,
                              'ratio_illustrative': float(U0*e0**4)}}
    Path(__file__).with_name('endpoints.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
