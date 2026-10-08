# DECISIONS: WAA_estm

What each step decided, with the numbers behind it, so that the step reports can be deleted by /tidy and rebuilt on
demand (rerun the script named in each heading). Rules of the frozen tail set: README.md, "Tail set v1". Settings:
section `waa` of `research/config.yaml`. TRAIN events only for every number unless stated.

## Step 0, inventory (`step0_inventory.py`, 2026-10-08)

- A(t) is **not stored**: it is `tl_mean` (release) - `baseline_db` (`baselines.csv`), exactly `baseline.a_of_t`.
- `baselines.csv` covers train + val events on the 283 train links (566 sub-links); test has no baseline, so there are
  no test tails. 7375 significant train pairs have a wet bin and `use_in_waa_main`.
- link_state is per link, 5 min; the ±30 min margin is on the dry side only (unknown bins, not wet). So t0
  (`last_wet_end`) is normally an unknown margin bin and the first dry bin comes about 30 min later.
- QC flags that block a sample are already excluded from `tl_mean`; `tsl_step` and `flatline` are information only.
- Decision: build the tails on 1 min `tl_mean` with the frozen baseline v2; no new baseline in this module.

## Step 1, all drying tails (`step1_extract_tails.py`, 2026-10-08)

- Tail = (event, sub-link) from t0 to the earliest of next wet bin / NaN run ≥ 15 min / t0 + 12 h. **No selection.**
- 20 355 tails (19 659 non-empty) on 566 sub-links, 69 events; train 15 370, val 4 985. A(t) equals `a_of_t`
  exactly on 3 random pairs.
- Return to the dry level (A < max(0.3 dB, 2 offset_std)): A_start > 2 dB tails 25 % back at 10 min, 59 % at 30 min,
  79 % at 3 h; A_start < 0.5 dB tails 94 % at 10 min. Strong tails dry slowly: WAA-like memory exists.
- 670 tails start on a wet bin (`last_wet_end` is searched only inside the network event); 566 are clipped at the
  season end. Both kept as rows, excluded later.
- Dry-weather rises ≥ 0.3 dB later in the tail: 56 % of tails, partly RSL quantization (upper bound).
- Decision: `tails_1min.parquet` (47 MB) is git-ignored and rebuilt by this script.

## Step 2, selection and fast-drop test (`step2_select_tails.py`, 2026-10-08)

- Rules fixed before the run (`waa.selection`); no tail kept or dropped by its shape. VAL not looked at.
- Excluded 4023 of 15 370 (mostly `use_in_waa_main` False: 3421). Truncation: dry-weather rise (≥ 0.6 dB for ≥ 30
  min) cut 1249 tails (5.5 % of tail time), a TSL step 29.
- Main set (A_start ≥ 1 dB, ≥ 60 min usable): 959 tails, 373 sub-links, 33 events; control (< 0.5 dB) 8838.
- Fast drop (A_start - A at the first dry bin): median 1.11 dB = 83 % of A_start. It does **not** grow with L
  (Spearman ≥ 25 GHz -0.10 [-0.20, 0.00]; < 25 GHz 0.14 [-0.09, 0.31]), but the path-rain reference wet_A_max barely
  grows with L either (0.13), so the L test has little power. Sub-links 1 vs 2 agree (rho 0.87).
- Open after step 2: light rain after the first dry bin stays in 466 of 959 main tails -> fixed in step 2b.
- Kept as key outputs: `step2_fast_drop_vs_length.png`, `step2_correlations.csv`.

## Step 2b, fixes after the step 2 gallery (`step2b_fix_tails.py`, 2026-10-08)

- Start moved to **t0_trace** = end of the last 5 min bin with any rain trace (radar R_max > 0 on the path, or a gauge
  within 1 km > 0) between t0 and the first dry bin.
- 381 of the 959 step 2 main tails fall below A_start = 1 dB at t0_trace: much of the step 2 "fast drop" happened
  while the radar or a gauge still showed rain. So the fast drop is not evidence for WAA.
- Two ends: lenient (step 2 ends, rise cut moved back to the last quiet minute, excess < 0.15 dB) and strict (also the
  first rain trace after t0_trace). Same 479 tails in both; 220 strict tails end at a rain trace.
- Noise exclusion: sub-links with dry-weather `offset_std` > 0.3 dB (one RSL step): 39 sub-links, 1077 tails.
- `level_shift` (|A| at 6-12 h > max(0.3, 2 offset_std)) and `undershoot` (A < -0.3 dB within 3 h) are flags, not
  exclusions; undershoot looks like a lasting level near -0.3 dB (baseline offset for that event?), not a dip.

## Step 2c, frozen tail set v1 (`step2c_freeze.py`, 2026-10-08; its report is a key output)

- Radar NaN on the path counts as a rain trace (unknown); rise cut v2 on all minutes from t0_trace (10 min median ≥ 0.6
  dB above its running minimum for ≥ 10 min).
- **Frozen set v1**: main 429 train / 123 val tails; alternative (A_start ≥ 0.6 dB) 1155 / 458.
- 3 events give 246 of the 429 train tails: resample or weight by event in the model step.
- Open: night-time rises without any rain trace (E20150701_0030): drizzle under the radar beam or condensation? RH is
  missing for those hours, so this is not decided.
