"""Laurent coefficient routines from the companion paper Two shorter intervals of even Dirichlet beta values."""


N = 96
ROOTS = range(-239, 240, 2)
ORDERS = [min(a, N - a, 30) + 1 for a in range(N + 1)]


def unit_factors(a):
    x = N - 2 * a
    factors = [(x - v, 3) for v in ROOTS]
    factors.extend((2 * (b - a), -order)
                   for b, order in enumerate(ORDERS) if b != a)
    if x:
        factors.append((x, 1))
    return factors


def coefficients_by_recurrence(a, q, p):
    order = q * ORDERS[a]
    shift = q if 2 * a == N else 0
    length = order - shift
    constant = 1
    logarithmic = [0] * (length - 1)
    for value, exponent in unit_factors(a):
        exponent *= q
        constant = constant * pow(value, exponent, p) % p
        inverse = pow(value, -1, p)
        term = inverse
        for j in range(length - 1):
            logarithmic[j] = (logarithmic[j] + exponent * term) % p
            term = -term * inverse % p
    coefficients = [constant]
    for k in range(1, length):
        derivative = sum(logarithmic[j] * coefficients[k - 1 - j] for j in range(k))
        coefficients.append(derivative * pow(k, -1, p) % p)
    return list(reversed([0] * shift + coefficients))


def coefficients_by_products(a, q, p):
    from flint import fmpz_mod_poly_ctx

    ctx = fmpz_mod_poly_ctx(p)
    order = q * ORDERS[a]
    x = N - 2 * a
    numerator = ctx.one()
    for v in ROOTS:
        numerator = numerator.mul_low(ctx([x - v, 1]), order)
    numerator = numerator.pow_trunc(3, order).mul_low(ctx([x, 1]), order)
    denominator = ctx.one()
    for b, pole_order in enumerate(ORDERS):
        if b != a:
            factor = ctx([2 * (b - a), 1]).pow_trunc(pole_order, order)
            denominator = denominator.mul_low(factor, order)
    unit = numerator.mul_low(denominator.inverse_series_trunc(order), order)
    jet = unit.pow_trunc(q, order)
    return [int(jet[order - i]) for i in range(1, order + 1)]
