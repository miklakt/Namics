#!/usr/bin/env python3
"""Compare Namics against analytic Flory-Huggins and Scheutjens-Fleer SCF results.

Run:  python3 tests/analytic_benchmarks.py [section ...]
Sections: slit graft adsorption depletion confinement binodal tension (default: all)

Lattice units, simple cubic (lambda0 = 4/6, lambda1 = 1/6). Single-chain cases use a dilute chain
(free with phi_b = 1e-10, or restricted with tiny n) so the mean field is negligible and the
propagator is an exact lattice random walk in the imposed walls/fields.
"""
import json, math, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin" / "namics"
WORK = Path(tempfile.mkdtemp(prefix="namics_analytic_"))
L0, L1 = 4 / 6, 1 / 6
FAILED = []

TAIL = """newton : isaac : method : pseudohessian
newton : isaac : tolerance : 1e-10
output : json : append : false
"""
SOLVENT = """mon : S : freedom : free
mol : sol : composition : (S)1
mol : sol : freedom : solvent
"""


def lattice(n, lo, up, geometry="planar"):
    return (f"lat : L : gradients : 1\nlat : L : geometry : {geometry}\nlat : L : lattice_type : simple_cubic\n"
            f"lat : L : n_layers : {n}\nlat : L : lowerbound : {lo}\nlat : L : upperbound : {up}\n")


def run(text):
    """Run one input; return (last problem dict, stdout) or (None, stdout) when no output was written."""
    path = WORK / "case.in"
    path.write_text(text)
    out = WORK / "case.output.json"
    out.unlink(missing_ok=True)
    proc = subprocess.run([str(BIN), str(path)], capture_output=True, text=True, cwd=WORK)
    if not out.exists():
        return None, proc.stdout
    return json.loads(out.read_text())["problems"][-1], proc.stdout


def check(label, value, reference, rtol):
    err = abs(value / reference - 1)
    ok = err <= rtol
    if not ok:
        FAILED.append(label)
    print(f"  [{'ok' if ok else 'FAIL'}] {label}: namics={value:.10g} reference={reference:.10g} rel.err={err:.1e} (tol {rtol:.0e})")


def dilute_chain_gn(N, n, *, chi_wall=0.0, geometry="planar", lo="surface", up="mirror", wall="lowerbound",
                    ensemble="free", extra_output=""):
    amount = "mol : pol : phibulk : 1e-10" if ensemble == "free" else "mol : pol : n : 1e-12"
    text = lattice(n, lo, up, geometry) + f"""mon : W : freedom : frozen
mon : W : frozen_range : {wall}
mon : A : freedom : free
mon : A : chi_W : {chi_wall}
mol : pol : composition : (A){N}
mol : pol : freedom : {ensemble}
{amount}
""" + SOLVENT + TAIL + "json : mol : pol : GN\n" + extra_output
    d, log = run(text)
    if d is None:
        raise RuntimeError(log[-1500:])
    return d["mol"]["pol"]


def slope_ln_gn(n_short, n_long, **kw):
    g1 = dilute_chain_gn(n_short, **kw)["GN"]
    g2 = dilute_chain_gn(n_long, **kw)["GN"]
    return (math.log(g2) - math.log(g1)) / (n_long - n_short)


# ---------------------------------------------------------------- single chain, exact lattice results

def section_slit():
    print("== ideal chain in a slit of D layers between absorbing walls (exact): Lambda = l0 + 2 l1 cos(pi/(D+1))")
    for D in (5, 10, 20):
        def gn(N):
            text = lattice(D, "surface", "surface") + f"""mon : W : freedom : frozen
mon : W : frozen_range : lowerbound
mon : V : freedom : frozen
mon : V : frozen_range : upperbound
mon : A : freedom : free
mol : pol : composition : (A){N}
mol : pol : freedom : free
mol : pol : phibulk : 1e-10
""" + SOLVENT + TAIL + "json : mol : pol : GN\n"
            return run(text)[0]["mol"]["pol"]["GN"]
        s = (math.log(gn(2000)) - math.log(gn(1000))) / 1000
        check(f"D={D} ln Lambda", s, math.log(L0 + 2 * L1 * math.cos(math.pi / (D + 1))), 1e-6)
        # full transfer-matrix GN at finite N (independent implementation)
        q = np.ones(D)
        for _ in range(299):
            p = np.pad(q, 1)
            q = L0 * q + L1 * (p[:-2] + p[2:])
        check(f"D={D} GN(N=300) vs transfer matrix", gn(300), q.sum(), 1e-6)


def section_graft():
    print("== chain grafted at an absorbing wall (exact image method on the lattice)")
    def trinomial(steps, half):
        p = np.zeros(2 * half + 1)
        p[half] = 1
        for _ in range(steps):
            p = L0 * p + L1 * (np.roll(p, 1) + np.roll(p, -1))
        return p
    for N in (50, 200, 800):
        n = N + 20
        text = lattice(n, "surface", "mirror") + f"""mon : W : freedom : frozen
mon : W : frozen_range : lowerbound
mon : X : freedom : pinned
mon : X : pinned_range : firstlayer
mon : A : freedom : free
mon : E : freedom : free
mol : pol : composition : (X)1(A){N - 2}(E)1
mol : pol : freedom : restricted
mol : pol : n : 1e-10
""" + SOLVENT + TAIL + "json : mol : pol : GN\njson : mon : E : phi\n"
        d, _ = run(text)
        z = np.arange(1, n + 1)
        p = trinomial(N - 1, n + 2)
        c = n + 2
        image = np.array([p[c + zz - 1] - p[c + zz + 1] for zz in z])   # start z=1, image source at z=-1
        phi_end = np.array(d["mon"]["E"]["phi"])
        check(f"N={N} survival GN", d["mol"]["pol"]["GN"], image.sum(), 1e-5)
        check(f"N={N} <z_end>", (phi_end * z).sum() / phi_end.sum(), (image * z).sum() / image.sum(), 1e-6)
        print(f"         (continuum Rayleigh limit sqrt(pi l1 N) = {math.sqrt(math.pi * L1 * N):.4f})")


def section_adsorption():
    print("== adsorbing wall, ground state (exact SF lattice): l1 a^2 + l0 a = l1/(w-1), Lambda = w (l0 + l1 a), w = exp(-chi l1)")
    for chi in (-1.5, -2.0, -3.0):
        w = math.exp(-chi * L1)
        a = (-L0 + math.sqrt(L0 ** 2 + 4 * L1 ** 2 / (w - 1))) / (2 * L1)
        s = slope_ln_gn(1000, 2000, n=300, chi_wall=chi, ensemble="restricted")
        check(f"chi_W={chi} ln Lambda", s, math.log(w * (L0 + L1 * a)), 1e-5)
    chi_c = math.log(1 - L1) / L1
    lo = slope_ln_gn(4000, 8000, n=300, chi_wall=chi_c + 0.05, ensemble="restricted")
    hi = slope_ln_gn(4000, 8000, n=300, chi_wall=chi_c - 0.05, ensemble="restricted")
    ok = lo < 0 < hi
    if not ok:
        FAILED.append("critical adsorption bracket")
    print(f"  [{'ok' if ok else 'FAIL'}] critical chi_W = ln(1-l1)/l1 = {chi_c:.5f}: dlnGN/dN = {lo:+.1e} at chi_c+0.05, {hi:+.1e} at chi_c-0.05")


# ---------------------------------------------------------------- single chain, asymptotic (Gaussian) results

def section_depletion():
    print("== depletion of free ideal chains at a hard wall: d(delta)/d(sqrt N) -> 2 sqrt(l1/pi)")
    delta = {}
    for N in (400, 1600, 6400):
        n = int(12 * math.sqrt(N * L1)) + 20
        delta[N] = -dilute_chain_gn(N, n, extra_output="json : mol : pol : theta_exc\n")["theta_exc"] / 1e-10
    check("slope N=1600..6400", (delta[6400] - delta[1600]) / (80 - 40), 2 * math.sqrt(L1 / math.pi), 1e-4)


def section_confinement():
    print("== confinement free energy per segment -> l1 k^2 / R_eff^2 (fit k^2 and wall offset from R = 40, 80)")
    for geometry, k2, label in (("planar", math.pi ** 2 / 4, "half-slit (mirror midplane), k=pi/2"),
                                ("cylindrical", 2.404825557695773 ** 2, "cylinder, k=j01"),
                                ("spherical", math.pi ** 2, "sphere, k=pi")):
        y = []
        for R in (40, 80):
            y.append(-slope_ln_gn(8 * R * R, 16 * R * R, n=R, geometry=geometry, lo="mirror", up="surface", wall="upperbound"))
        slope = (1 / math.sqrt(y[1]) - 1 / math.sqrt(y[0])) / 40
        check(f"{label}: k^2", 1 / (slope ** 2 * L1), k2, 1e-3)


# ---------------------------------------------------------------- Flory-Huggins bulk and interface

def fh_binodal(NA, NB, chi):
    f = lambda p: p / NA * math.log(p) + (1 - p) / NB * math.log(1 - p) + chi * p * (1 - p)
    f1 = lambda p: (math.log(p) + 1) / NA - (math.log(1 - p) + 1) / NB + chi * (1 - 2 * p)
    F = lambda x: np.array([f1(x[0]) - f1(x[1]), (f(x[0]) - x[0] * f1(x[0])) - (f(x[1]) - x[1] * f1(x[1]))])
    # seed from the lower convex hull of f (avoids the trivial root phi1 = phi2)
    p = np.unique(np.concatenate([np.logspace(-14, -0.3, 4000), 1 - np.logspace(-14, -0.3, 4000)]))
    fv = np.array([f(q) for q in p])
    hull = []
    for i in range(len(p)):
        while len(hull) >= 2 and (fv[hull[-1]] - fv[hull[-2]]) * (p[i] - p[hull[-2]]) >= (fv[i] - fv[hull[-2]]) * (p[hull[-1]] - p[hull[-2]]):
            hull.pop()
        hull.append(i)
    k = max(range(len(hull) - 1), key=lambda j: p[hull[j + 1]] - p[hull[j]])
    x = np.array([p[hull[k]], p[hull[k + 1]]])
    for _ in range(100):
        J = np.zeros((2, 2))
        h = 1e-9 * np.minimum(x, 1 - x)
        for j in range(2):
            e = np.zeros(2)
            e[j] = h[j]
            J[:, j] = (F(x + e) - F(x - e)) / (2 * h[j])
        dx = np.linalg.solve(J, -F(x))
        while np.any(x + dx <= 0) or np.any(x + dx >= 1):
            dx /= 2
        x += dx
        if np.abs(dx).max() < 1e-15:
            break
    return x


def two_phase_box(NA, NB, chi, n, theta_A):
    """A/B box between mirrors, interface seeded mid-box. Returns (phi_A, omega) or None."""
    z = np.arange(1, n + 1)
    for restricted in ("A", "B"):
        for amp, width in ((0.5, 3.0), (0.1, 3.0), (0.5, 8.0), (0.1, 8.0), (1.0, 3.0), (2.0, 3.0), (0.25, 3.0)):
            step = np.tanh((z - n / 2) / width)
            (WORK / "guess.json").write_text(json.dumps({"initial_guess": {"u": {"mon": {"A": list(amp * step), "B": list(-amp * step)}}}}))
            fa, fb = ("restricted", "solvent") if restricted == "A" else ("solvent", "restricted")
            theta = theta_A if restricted == "A" else n - theta_A
            text = lattice(n, "mirror", "mirror") + f"""mon : A : freedom : free
mon : A : chi_B : {chi}
mon : B : freedom : free
mol : A : composition : (A){NA}
mol : A : freedom : {fa}
mol : B : composition : (B){NB}
mol : B : freedom : {fb}
mol : {restricted} : theta : {theta}
sys : noname : initial_guess : file
sys : noname : guess_inputfile : guess.json
""" + TAIL + "json : mol : A : phi\njson : sys : noname : grand_potential_density\n"
            d, _ = run(text)
            if d is None:
                continue
            phi = np.array(d["mol"]["A"]["phi"])
            om = np.array(d["sys"]["noname"]["grand_potential_density"])
            # reject non-physical or non-separated solutions (Namics may write output for unconverged runs)
            if phi.max() > 1 + 1e-9 or phi.min() < -1e-9 or abs(phi[0] - phi[-1]) < 1e-3 or abs(om[0] - om[-1]) > 1e-2:
                continue
            return phi, om
    return None


def section_binodal():
    print("== FH coexistence from a phase-separated box; interface placed where both plateaus have equal pressure")
    print("   (on a lattice the interface feels a one-layer-periodic pinning force; at fixed theta it can hold a")
    print("    small pressure difference, so theta is tuned by a sub-layer shift until omega(1) = omega(n))")
    for NA, NB, chi, n in ((20, 20, 0.15, 200), (50, 1, 0.8, 200), (100, 1, 0.7, 300),
                           (10, 30, 0.25, 200), (30, 10, 0.25, 200), (5, 1, 1.5, 200)):
        b = fh_binodal(NA, NB, chi)
        th0, dphi = n * (b[0] + b[1]) / 2, b[1] - b[0]
        cache = {}
        def g(k):
            if k not in cache:
                r = two_phase_box(NA, NB, chi, n, th0 + k * dphi)
                cache[k] = (math.nan, None) if r is None else (r[1][0] - r[1][-1], r[0])
            return cache[k]
        ks = np.linspace(0, 1, 9)
        vals = [g(k)[0] for k in ks]
        best = (math.inf, None)
        for i in [i for i in range(8) if vals[i] * vals[i + 1] <= 0]:
            a, c, fa, fc = ks[i], ks[i + 1], vals[i], vals[i + 1]
            for it in range(80):
                m = (a + c) / 2 if it % 3 == 2 or fc == fa else c - fc * (c - a) / (fc - fa)
                fm, phi = g(m)
                if abs(fm) < abs(best[0]):
                    best = (fm, phi)
                if abs(fm) < 1e-11 or c - a < 1e-12:
                    break
                if fa * fm < 0:
                    c, fc = m, fm
                else:
                    a, fa = m, fm
            if abs(best[0]) < 1e-10:
                break
        if best[1] is None:
            FAILED.append(f"binodal {NA}/{NB}")
            print(f"  [FAIL] NA={NA} NB={NB} chi={chi}: no converged two-phase solution")
            continue
        lo, hi = sorted((best[1][-1], best[1][0]))
        check(f"NA={NA} NB={NB} chi={chi} phi'", lo, b[0], 1e-7)
        check(f"NA={NA} NB={NB} chi={chi} phi''", hi, b[1], 1e-7)


def section_tension():
    print("== symmetric blend interfacial tension vs Helfand-Tagami sqrt(chi/6) (limit chiN -> inf, then chi -> 0)")
    ratios = {}
    for chi in (0.1, 0.075, 0.05):
        xs, gs = [], []
        for chiN in (10, 14, 18, 22):        # chiN > ~30 loses phi_B = 1 - phi_A to cancellation in double precision
            N = int(round(chiN / chi))
            n = int(10 * math.sqrt(N / 6)) + 60
            r = two_phase_box(N, N, chi, n, n / 2)
            if r is None:
                continue
            xs.append(1 / (chi * N))
            gs.append((r[1] - r[1][0]).sum())    # gamma = sum_z [omega(z) - omega_bulk]
        if len(xs) < 3:
            FAILED.append(f"tension chi={chi}")
            print(f"  [FAIL] chi={chi}: too few converged points")
            continue
        g_inf = np.linalg.lstsq(np.vander(np.array(xs), 3), np.array(gs), rcond=None)[0][-1]
        ratios[chi] = g_inf / math.sqrt(chi / 6)
        print(f"  chi={chi}: gamma(chiN=10..22)={[round(v, 5) for v in gs]} -> gamma_inf={g_inf:.5f}; HT={math.sqrt(chi / 6):.5f}; ratio={ratios[chi]:.4f}")
    if len(ratios) >= 2:
        chis = sorted(ratios)
        fit = np.polyfit(chis, [ratios[c] for c in chis], 1)     # lattice correction is O(chi)
        check("gamma_inf / gamma_HT extrapolated to chi -> 0", fit[1], 1.0, 1e-2)


SECTIONS = {"slit": section_slit, "graft": section_graft, "adsorption": section_adsorption,
            "depletion": section_depletion, "confinement": section_confinement,
            "binodal": section_binodal, "tension": section_tension}

if __name__ == "__main__":
    if not BIN.exists():
        sys.exit(f"build first: {BIN} not found")
    for name in sys.argv[1:] or SECTIONS:
        SECTIONS[name]()
    print(f"\n{len(FAILED)} failed" + (": " + ", ".join(FAILED) if FAILED else ""))
    sys.exit(1 if FAILED else 0)
