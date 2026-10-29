"""Which heads fail the energy cross-check in cf_par, and why?"""
import os, sys, numpy as np, h5py
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_minimal_wrap_core as core
p = sys.argv[1]
src = open(core.__file__).read()
# re-run the core's reconstruction but capture the arrays instead of asserting
src = src.replace("        assert ok.all(), (", "        globals()['_DIAG'] = dict(ok=ok, mine=mine, ref=ref, rho=rho, z=z, st=hs[fi,mi,hi], hst=hst[fi,mi,hi], loaded=loaded, fi=fi)\n        assert True or ok.all(), (")
ns = {"__name__": "diag"}
exec(compile(src, core.__file__, "exec"), ns)
ns["build_run"](p, 1)
D = ns["_DIAG"]
bad = ~D["ok"]
print("failing heads: %d / %d" % (bad.sum(), bad.size))
err = np.abs(D["mine"] - D["ref"])
print("%5s %6s %10s %10s %12s %14s %12s" % ("frame","state","rho nm","z nm","ref energy","abs err","err/ref"))
for k in np.argsort(-err)[:12]:
    print("%5d %6d %10.4f %10.5f %12.4e %14.3e %12.2e" % (D["fi"][k], D["st"][k], 1000*D["rho"][k],
          1000*D["z"][k], D["ref"][k], err[k], err[k]/max(1e-30, abs(D["ref"][k]))))
r = D["ref"]
print("\nenergy scale max |ref| = %.4f ; abs err of failing heads / that scale: max %.2e"
      % (np.abs(r).max(), err[bad].max()/np.abs(r).max()))
print("failing heads by state:", {int(s): int(((D['st']==s)&bad).sum()) for s in (1,2)})
print("failing heads: loaded extension range %.4f..%.4f nm" % (1000*D["loaded"][bad].min(), 1000*D["loaded"][bad].max()))
