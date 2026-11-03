"""Exact two-sided p value of a Spearman (or Pearson) correlation test, also when it underflows double precision.

scipy's test: t = r sqrt((n - 2) / (1 - r^2)), p = P(|T| > t) with n - 2 degrees of freedom
             = I_x(df / 2, 1 / 2), x = df / (df + t^2)   (regularized incomplete beta)
For tiny p, log I_x(a, b) is computed from I_x(a, b) = x^a (1 - x)^b / (a B(a, b)) * 2F1(a + b, 1; a + 1; x)
(DLMF 8.17.8), with the hypergeometric series summed term by term in double precision.
"""
import math

from scipy import stats
from scipy.special import gammaln


def log10_p(r, n):
    df = n - 2
    t2 = r * r * df / ((1 - r) * (1 + r))
    p = 2 * stats.t.sf(math.sqrt(t2), df)
    if p > 1e-280:
        return math.log10(p)
    x, a, b = df / (df + t2), df / 2, 0.5
    s, term, k = 0.0, 1.0, 0
    while term > 1e-17 * s or k < 3:
        s += term
        term *= x * (a + b + k) / (a + 1 + k)
        k += 1
    ln_p = a * math.log(x) + b * math.log1p(-x) - math.log(a) - (gammaln(a) + gammaln(b) - gammaln(a + b)) + math.log(s)
    return ln_p / math.log(10)


def p_text(r, n, mathtext=True):
    """'p = 3.2 × 10^-123456' (matplotlib mathtext) or plain text"""
    l10 = log10_p(r, n)
    e = math.floor(l10)
    m = 10 ** (l10 - e)
    if m >= 9.95:
        m, e = 1.0, e + 1
    if e >= -3:
        return f"p = {10 ** l10:.2g}"
    if mathtext:
        return rf"p = {m:.1f} $\times$ 10$^{{{e}}}$"
    return f"p = {m:.1f} × 10^{e}"


if __name__ == "__main__":
    # check against scipy where scipy does not underflow
    for r, n in ((0.05, 20000), (0.3, 3000), (0.5, 2000), (0.2, 100000)):
        exact = stats.t.sf(math.sqrt(r * r * (n - 2) / (1 - r * r)), n - 2) * 2
        df = n - 2
        t2 = r * r * df / (1 - r * r)
        x, a, b = df / (df + t2), df / 2, 0.5
        s, term, k = 0.0, 1.0, 0
        while term > 1e-17 * s or k < 3:
            s += term
            term *= x * (a + b + k) / (a + 1 + k)
            k += 1
        series = (a * math.log(x) + b * math.log1p(-x) - math.log(a) - (gammaln(a) + gammaln(b) - gammaln(a + b))
                  + math.log(s)) / math.log(10)
        sci = math.log10(exact) if exact > 0 else float("-inf")
        print(f"r={r} n={n}: scipy log10 p {sci:.6f}  series {series:.6f}")
    for r, n in ((0.850236, 906513), (0.685004, 39524), (0.6337, 906513)):
        print(r, n, p_text(r, n, mathtext=False))
