# Methods: cmlpipe processing of the OpenMRG dataset

Final definitions of every product of the pipeline (cmlpipe 0.2.0, configuration `configs/openmrg.yaml`). Every
threshold below is a key of that file; the code is in `src/cmlpipe/`. Variable-level documentation:
`docs/data_dictionary.md` (generated with the release).

**Data.** OpenMRG: Andersson, J. C. M., Olsson, J., van de Beek, R. (C. Z.), and Hansryd, J. (2022): OpenMRG: Open data
from Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden, Earth Syst. Sci. Data, 14,
5411-5426, https://doi.org/10.5194/essd-14-5411-2022. Data v1.1: https://doi.org/10.5281/zenodo.7107689 (CC BY-SA 4.0).
Period 1 June - 31 August 2015 (JJA), 364 links / 728 sub-links, 10 City gauges (1 min), SMHI gauge GbgA (15 min),
hourly SMHI weather, NORDRAD radar composite (5 min, 2 km).

## 1. Conventions

- **Time**: UTC, tz-naive. Every value is labelled by the START of its interval, bin = [t, t + dt) ("interval_start").
  Sources labelled at the interval end are shifted back once at ingest (City gauges -60 s, SMHI precipitation -900 s).
  Radar stamps and SMHI weather are instants and keep their stamp; a radar stamp t is assigned to the bin
  [t, t + 5 min) (gauge vs radar correlation peaks at lag 0). CML samples are labelled by the start of the 10 s SNMP
  window.
- **Keys**: `cml_id` (link, raw "Link"), sub-link 1 transmits site_0 -> site_1 (raw direction A), sub-link 2 the
  reverse; in the release `sublink_id` = "<cml_id>_<1|2>". `gauge_id`, `event_id` = "E" + start ("E20150728_1455"),
  `run_id` = "<event_id>_Rnn".
- **Units**: dBm (tsl, rsl), dB (tl = tsl - rsl), mm, mm/h, metres in EPSG:3006 for every length and distance, GHz.
- **Levels**: L1 = raw values in a common schema (sentinel 1e10 -> NaN, regular 10 s grid, tl added), int16 x 0.1 dB
  (lossless: all raw values are multiples of 0.1 dB). L2 = QC flags (never edits L1). L3 = aggregates.
- **Reproducibility**: every output carries `stage`, `stage_hash` (sha256 of only the config keys that stage uses; no
  dataset path is hashed), `config_hash`, `code_version` (git commit) and `created_utc`. A stage skips outputs with the
  same stage hash and the same inputs; `--force` rebuilds. Run log: `run_manifest.json`.

## 2. QC flags (L2)

One uint16 bitmask per 10 s sample (`cmlpipe.qc.flags.QCFlag`). FLAG, DON'T FIX: nothing is removed or filled. Flags use
only the CML signals (never radar, gauges or weather). Days are processed with 12 h of context on each side.

| bit | flag | rule (config `qc.flags`) | why |
|---|---|---|---|
| 0 | missing_rsl | rsl is NaN | missing data |
| 1 | missing_tsl | tsl is NaN | missing data |
| 2 | out_of_range | rsl outside [-100, -20] dBm or tsl outside [-15, 30] dBm | safety net beyond the published extremes |
| 3 | rx_floor | rsl <= floor + 1 dB for >= 10 min, where the sub-link floor (lowest level with >= 30 samples within 1 dB) is >= 20 dB below its median and <= -87 dBm | receiver clipping: the attenuation is unknown, only a lower bound |
| 4 | spike | \|rsl - centred 5 min rolling median\| > max(6 scaled MAD, 5 dB), run of <= 3 samples, the samples before and after agree within 1.5 dB and lie within 1.5 dB of the rolling median | isolated outliers, not the sharp bottom of a fade |
| 5 | tsl_step | \|tsl - previous valid tsl\| >= 1 dB | transmit-power change (informational) |
| 6 | flatline | rsl exactly constant over >= 12 h of valid samples | frozen reporting (informational) |
| 7 | network_gap | sample NaN while >= 5 % of the alive sub-links (>= 50 % valid that day) are NaN | collector outage, not a link problem |

**Blocking sets.** L2 summaries (`qc.blocking_flags`): missing_rsl, missing_tsl, out_of_range, rx_floor, spike,
network_gap. Aggregation (`aggregate.blocking_flags`): rsl: missing_rsl, out_of_range, rx_floor, spike, network_gap;
tsl: missing_tsl, out_of_range, network_gap; tl: the union. tsl_step and flatline never block.

## 3. Quality status per sub-link

`tables/inventory/quality_status.csv` (config `inventory.status`), first rule that applies:

- **unusable**: JJA tl completeness < 50 %, or > 20 % of the samples present (valid rsl and tsl) blocking-flagged, or
  rx_floor > 20 %.
- **degraded**: a month with tl completeness < 80 %; rx_floor > 5 %; spike > 0.1 %, flatline > 10 % or tsl_step > 1 % of
  the samples; std of the dry-period 5 min tl means above the network p95; p90 - p10 of the dry-period tl medians above
  the network p95. The reason string lists every rule that fired.
- **good**: otherwise. (The step 3 status, `sublink_qc_status.csv`, is kept for reference.)

**TSL behaviour** (f = share of valid samples with a TSL change >= 1 dB; d = changes per day of data): no_tsl;
toggling f > 1 %; constant: >= 99 % of samples at one value and d <= 0.1; rare_steps d <= 2; frequent_steps otherwise.
`tsl_follows_fades_descriptive` (>= 5 steps in events and an in-event step rate >= 3 x the dry rate) is descriptive only.

## 4. Geometry

All lengths and distances in EPSG:3006 (SWEREF 99 TM). Link ends within 10 m are one site (418 sites). Radar path
weights: the straight path is intersected with the radar grid in the radar's own polar stereographic projection
(`+proj=stere +lat_ts=60 +ellps=bessel +lon_0=14 +lat_0=90`, the only one where cells are exact 2 km squares); the path
is cut at every cell edge and each piece is assigned to the cell of its midpoint (no double counting on edges or
corners); a weight is the fraction of the path length in the cell (fractions are distortion-free). Gauge-to-path
distance is point-to-segment, position 0 = site_0, 1 = site_1. Link / gauge pairs are kept up to 5 km; neighbours up
to 5 km between paths.

## 5. Aggregation (L3)

Resolutions 1 min, 5 min, 15 min, 1 h; every one divides a day, so no bin crosses midnight. A 10 s sample goes to the
bin containing its label.
- **CML**: per bin and variable (tl, rsl, tsl): mean, median, min, max, std (ddof 0; no rsl / tsl std at 1 min) of the
  samples without a blocking flag; n_valid_<v> always kept; all statistics of v are NaN when n_valid < ceil(0.8 x
  samples per bin). flags_or = bitwise OR of all flags in the bin; n_tsl_step = tsl_step samples in the bin. float32.
- **Gauges**: SUMS of the native amounts; NaN unless every native step of the bin is valid (strict). rate = amount /
  dt. City 1 min - 1 h, plus the peak 1 min rate; SMHI 15 min and 1 h only (never disaggregated); qc_worst for SMHI.
- **Radar**: dBZ < 0 -> rain rate 0 (the "no echo" level would otherwise accumulate phantom rain); Z = 200 R^1.5 (file
  attributes). Along each path: path-weighted mean (R_mean), weighted std, max cell (R_max), wet_fraction (share of path
  weight > 0.1 mm/h), valid_weight_fraction; a stamp is valid for a link when >= 50 % of its path weight has data. Bins:
  mean of the valid stamps (max for R_max), NaN when fewer than 60 % of the stamps are valid. At gauges: the cell
  containing the gauge and the 3 x 3 mean. 5 min (native), 15 min, 1 h; never at 1 min.

## 6. Events, dry periods, significance, runs

From the REFERENCES only (radar + 10 City gauges), never from CML; SMHI GbgA is the independent check (100 % of its
rain falls inside events).
- **Domain**: the 720 radar cells whose centre is within 3 km of a link path.
- **State per 5 min bin**: radar valid if >= 50 % of the domain cells have data. Wet if (radar valid and >= 1 % of the
  valid domain cells have a floored rate >= 0.1 mm/h) or >= 1 City gauge has >= 0.1 mm. Dry if not wet, radar valid
  and no wet domain cell. Unknown otherwise (radar missing, or 1-7 wet cells); unknown never counts as dry, and a gauge
  with a missing bin is "gauge unknown".
- **Events**: wet bins merged when separated by < 2 h of non-wet time; events < 15 min dropped; an event within 2 h of
  the period start / end is flagged. Intensity class: terciles of the maximum 5 min City rate over events.
- **Dry periods**: dry bins more than 30 min from any wet bin, in runs >= 6 h.
- **Significance** (`event_significance`, physics-based, not tuned on CML): significant if the mean City gauge total
  >= 1 mm, or the domain-mean radar accumulation >= 1 mm, or the max 5 min City rate >= 5 mm/h; heavy if significant and
  (total >= 10 mm or rate >= 30 mm/h); otherwise none_drizzle.
- **Wet runs**: the wet runs of the state series inside each event (before merging). A run is significant by the same
  rule with both totals lowered to 0.5 mm; the wet core spans the first to the last significant run.

## 7. Per-link reference state (5 min)

From the references only (`link_state`): **wet** if R_max along the path >= 0.1 mm/h or a gauge within 1 km of the path
has >= 0.1 mm in the bin; **dry** if not wet, R_max == 0 in the bin and every bin within +-30 min (no radar gap), and
every gauge within 1 km is 0 and valid over the same window; **unknown** otherwise. Stored as int8 (0 dry, 1 wet,
-1 unknown) with a reason bitmask (1 radar_wet, 2 gauge_wet, 4 radar_missing in the window, 8 rain in the margin).
City gauges at 5 min; the SMHI 15 min amount labels its three 5 min bins. Over JJA: 77 % dry, 9 % wet, 14 % unknown
link-bins; inside significant events 30 % wet; outside events 0 % wet; when the network state is dry, 97 % of the
link-bins are dry.

## 8. Split (frozen)

Events: train / val / test = 60 / 20 / 20 %, stratified by intensity class (seeded RNG, seed 42, over the sorted ids).
Links: train / test = 80 / 20 %, stratified by length class (< 1, 1-3, > 3 km) x frequency band (< 25, >= 25 GHz of
sub-link 1); parallel links (same pair of sites) share a split; links with both sub-links unusable at step 3 are
excluded. Evaluation: fit on train events x train links (tune on val events); unseen events = test events x train
links; unseen links = train events x test links; both = test x test. The link split is NOT a spatial hold-out
(neighbouring links see the same rain and share sites). The split files never change unless `--force-split`.

## 9. Known limitations

- The radar is below the gauge JJA totals at 10 of 11 gauge pixels (by 3-24 %) and missed the first Torslanda peak; it
  is gauge-adjusted, so it is not independent of the gauges.
- 87 of 364 links have a gauge within 1 km of the path; 62 links lie inside a single radar cell.
- Length and frequency are confounded (Pearson r(log10 length, frequency) = -0.57).
- 54 of 89 events are drizzle only; 79 events are merged from several wet runs, so network event bounds are coarse for
  per-link timing (use the wet runs or the per-link state).
- The status rules based on network quantiles mark about 5 % of the sub-links by construction.
- Not modelled yet: baseline, the attenuation A(t) = TSL - RSL - baseline, wet / dry classification of the CML,
  wet-antenna attenuation, the path law and rain retrieval.
