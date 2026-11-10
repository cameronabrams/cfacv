#!/usr/bin/env python3
"""Check cfacv's online collective-variable values against the trajectory.

cfacv computes its CVs inside NAMD's calcforces callback and reports them to
a TAMDoutputFile.  This recomputes the same CVs from the DCD the same run
wrote and compares, which is the "Implementation" validation gate: the online
values must agree with an independent measurement before any result computed
with them is interpreted.

Centres and CV definitions are read from the same files cfacv itself reads,
so the comparison cannot drift from the run it is checking.  Centres may be
multi-atom groups; those are reduced to a mass-weighted centre of mass, as
cfacv does.

Usage:
    posthoc_cvcheck.py --psf alad_wb.psf --dcd cvcheck.dcd \\
                       --cvout cvcheck.cvout --label label.pdb \\
                       --cv cv_check.inp --dcdfreq 10

Requires MDAnalysis.  Exits non-zero if any CV disagrees beyond --tol.
"""
import argparse, sys
import numpy as np

def read_centers(label_pdb):
    """Group atoms by the beta column, as cfacv's read_centersPDB does:
    beta is the 1-based group index, 0.00 means 'not in any centre'."""
    groups, serial = {}, 0
    for line in open(label_pdb):
        if not line.startswith("ATOM"):
            continue
        serial += 1                      # cfacv counts ATOM records in order
        beta = line[60:66].strip()
        try:
            b = float(beta)
        except ValueError:
            continue
        if b == 0.0:
            continue
        groups.setdefault(int(b) - 1, []).append(serial - 1)   # 0-based index
    return [groups[k] for k in sorted(groups)]

def read_cvs(cv_inp):
    cvs = []
    for line in open(cv_inp):
        p = line.split()
        if len(p) < 2 or p[0].startswith("#"):
            continue
        cvs.append((p[0].upper(), tuple(int(x) for x in p[1:])))
    return cvs

def bond(a, b):
    return np.linalg.norm(a - b)

def angle(a, b, c):
    u, v = a - b, c - b
    cosine = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return np.arccos(np.clip(cosine, -1.0, 1.0))

def dihed(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))

def rgyr(pts, masses):
    m = masses[:, None]
    com = (pts * m).sum(0) / masses.sum()
    return np.sqrt((masses * ((pts - com) ** 2).sum(1)).sum() / masses.sum())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psf", required=True)
    ap.add_argument("--dcd", required=True)
    ap.add_argument("--cvout", required=True)
    ap.add_argument("--label", default="label.pdb")
    ap.add_argument("--cv", default="cv.inp")
    ap.add_argument("--dcdfreq", type=int, required=True)
    ap.add_argument("--tol", type=float, default=1e-3,
                    help="max allowed |online - post hoc| (default 1e-3)")
    a = ap.parse_args()

    try:
        import MDAnalysis as mda
    except ImportError:
        sys.exit("MDAnalysis is required: pip install MDAnalysis")

    centers = read_centers(a.label)
    cvs = read_cvs(a.cv)
    print("centres: %d   CVs: %d" % (len(centers), len(cvs)))
    for i, g in enumerate(centers):
        print("  centre %d: %d atom(s)" % (i, len(g)))

    online = []
    for line in open(a.cvout):
        p = line.split()
        if p[1:2] == ["Th"]:
            online.append([int(p[2])] + [float(x) for x in p[3:]])
    if not online:
        sys.exit("no 'Th' records in %s -- set TAMDoutputlevel to include bit 2" % a.cvout)
    online = np.array(online)
    steps, onv = online[:, 0].astype(int), online[:, 1:]
    if onv.shape[1] != len(cvs):
        sys.exit("cvout has %d columns but %s defines %d CVs"
                 % (onv.shape[1], a.cv, len(cvs)))

    u = mda.Universe(a.psf, a.dcd)
    masses = u.atoms.masses

    post, matched = [], []
    for k, ts in enumerate(u.trajectory):
        step = (k + 1) * a.dcdfreq
        if step not in steps:
            continue
        pos = ts.positions
        # reduce each centre to its mass-weighted centre of mass
        c = np.array([(pos[g] * masses[g][:, None]).sum(0) / masses[g].sum()
                      for g in centers])
        row = []
        for typ, idx in cvs:
            pts = [c[i] for i in idx]
            if typ == "BOND":    row.append(bond(*pts))
            elif typ == "ANGLE": row.append(angle(*pts))
            elif typ == "DIHED": row.append(dihed(*pts))
            elif typ == "RGYR":  row.append(rgyr(pos[centers[idx[0]]], masses[centers[idx[0]]]))
            elif typ.startswith("CARTESIAN_"):
                row.append(pts[0]["XYZ".index(typ[-1])])
            else:
                sys.exit("unhandled CV type %s" % typ)
        post.append(row); matched.append(step)

    if not post:
        sys.exit("no DCD frame matched a CV report step; check --dcdfreq")
    post, matched = np.array(post), np.array(matched)
    onv_m = onv[np.isin(steps, matched)]
    print("compared %d frames\n" % len(matched))

    hdr = "  %-5s %-20s %-5s %13s %15s %12s %11s"
    print(hdr % ("cv", "definition", "unit", "online mean", "post-hoc mean",
                 "max |diff|", "RMS diff"))
    bad = []
    for j, (typ, idx) in enumerate(cvs):
        d = onv_m[:, j] - post[:, j]
        if typ == "DIHED":                      # compare on the circle
            d = (d + np.pi) % (2 * np.pi) - np.pi
        mx, rms = np.abs(d).max(), np.sqrt((d ** 2).mean())
        unit = "rad" if typ in ("DIHED", "ANGLE") else "A"
        if mx >= a.tol:
            bad.append((j, typ, mx))
        print(hdr % ("cv%d" % j, "%s%s" % (typ, idx), unit,
                     "%.5f" % onv_m[:, j].mean(), "%.5f" % post[:, j].mean(),
                     "%.2e" % mx, "%.2e" % rms) + ("   <-- MISMATCH" if mx >= a.tol else ""))
    print()
    if bad:
        print("DISAGREEMENT beyond tol=%g:" % a.tol)
        for j, typ, mx in bad:
            print("  cv%d (%s): max |diff| %.3e" % (j, typ, mx))
        return 1
    print("ALL CV TYPES AGREE (tol=%g)" % a.tol)
    return 0

if __name__ == "__main__":
    sys.exit(main())
