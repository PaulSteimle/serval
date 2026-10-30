# ATM weighting — progress

Branch: `ATM_weighting` (branched from `CARM_NIR_tel_calib`, incl. its uncommitted telluric-correction work).

Goal: replace the static telluric line mask by weights derived from the per-epoch fitted
telluric model, so the RV χ² uses what the atmosphere actually lets through.

## Status

- [x] Analysis of current telluric flow (this file)
- [ ] Model-error term in template, coadd, RV loop
- [ ] Per-epoch `T_min` cut replacing static mask
- [ ] CLI options `-atmerr`, `-atmtmin`
- [ ] Calibrate η from residuals vs. optical depth
- [ ] A/B test (e.g. Barnard's Star, CARM_NIR / CARM_VIS / HPF) vs. 0.25-limit mask

## Current flow (as found)

Telluric model is applied in three places in `src/serval.py`:

| Stage | Line (approx.) | What happens |
|---|---|---|
| Template (`spt`) | ~1446 | `fit_atm_par` on `atm_cal_order`, `f /= yatm`, `e /= yatm` |
| Coadd | ~1633 | same per spectrum and order (`calc_atm(..., order=o)`) |
| RV loop | ~2183 | same, `atm_par` stored on `spoklist[n]` |

In all three: `atm._calc_atm` sets `T < 0.02` to NaN → those pixels get `flag.atm`.

### Weighting already present (implicitly)
`polyreg` weights with `1/e²`. Since `e` is divided by `T`, the weight per pixel is
`(T / e_phot)²` → telluric pixels are already downweighted by `T²`, but only for photon noise.

### Hard cut by static mask
- `atmmask_file` resolves to `inst.atmspec_mask` when `atmspec` is set
  (`telluric_mask_CARM_NIR_0.25_limit.dat`, `telluric_mask_hpf_0.25_limit.dat`;
  CARM_VIS has no `atmspec_mask` → `telluric_mask_carm_short.dat`).
- RV loop: `b2[tellmask(w2)>0.01] |= flag.atm` → excluded from `pind`, i.e. removed from χ².
- Coadd: such pixels go into `tellind` → fixed weight `0.1/ntset/(e² + median(e)²)`.

### What is missing
1. Static mask comes from one reference atmosphere; does not follow per-epoch airmass /
   water vapour although `atm_par` is fitted per spectrum.
2. No term for model imperfection. Residual after division scales with optical depth
   `τ = -ln T`: `δf ≈ f_corr · δτ`. Deep lines are badly corrected but get photon-noise-only
   weight → that is why a mask is needed at all.

## Proposal

Error model with relative optical-depth error η:

```
σ² = (e/T)² + (η · f_corr · ln T)²     # f_corr = f/T
mask: T < T_min → flag.atm             # per epoch, replaces static 0.25-limit file
```

Patch (same 3 lines at template ~1450, coadd ~1647, RV loop ~2203, after the division):

```python
sp.f /= yatm
sp.e /= yatm
sp.e = np.sqrt(sp.e**2 + (atmerr * sp.f * np.log(yatm))**2)   # telluric model error
sp.bpmap[~(yatm > atmtmin)] |= flag.atm                          # also catches NaN
```

- New options: `-atmerr` (η, start ~0.05), `-atmtmin` (start 0.25 for comparison).
- With `atmspec` set, default `-atmmask ''`; otherwise the static hard cut overrides the weights.
- Coadd: without static mask `tellind` only holds model-NaN pixels (already excluded via
  `emod > 0`) → coadd uses the new σ directly. The commented "high S/N" branch at ~1748
  (`fcont*log(...)`) is essentially this formula already.
- RV error stays consistent: scaled by `rchi`, which now contains the model term.
- Calibrating η: bin RV-fit residuals by τ, choose η such that reduced χ² ≈ 1 in every bin.
  One η per instrument, possibly separate for O₂ / H₂O.
- Known ceiling: model error is correlated within a line, χ² treats pixels independently →
  fine for weighting, slightly optimistic RV errors in telluric-dominated orders.

## Side findings

- `serval.py` ~2213: `ok_o = ... & (spt.f[o] / sp.e[o] > 3)` — likely typo for `sp.f`
  (only affects the `looka` diagnostic plot).
- Sky mask is commented out in all `fit_atm_par` `ok` selections; in the NIR strong OH
  emission may bias the H₂O coefficient → A/B check worth doing.
- `atm.py`: `hdu[0].header['C1_REF']` raises `KeyError` if the keyword is absent, so the new
  `None` check only helps when the key exists with empty value. Use `.get('C1_REF')`.
