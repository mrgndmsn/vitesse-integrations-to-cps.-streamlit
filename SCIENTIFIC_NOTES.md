# Vitesse · Time & Pixels

Local Streamlit app for continuous Vitesse integrations, reviewed spot coordinates, configurable time bins, and NICE v7-layout CSV packages.

## Start

On this Mac, run `./run.command` from this directory, or double-click `run.command`. Visit http://localhost:8501. The launcher uses the task's installed dependencies when available; otherwise it creates a local virtual environment and installs requirements. The pinned, tested dependencies require Python 3.12 or newer. Set `VITESSE_PYTHON` to the appropriate Python executable if `python3` is older.

Portable setup:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

The app listens only on this computer (127.0.0.1). Usage telemetry is disabled. Original inputs are read-only. Folder/ZIP uploads and generated exports are staged in an operating-system temporary directory per session. Large datasets work best with the Local folder option, avoiding browser upload copies. Temporary staging folders are named `vitesse-app-*`; remove them after the app is stopped if you no longer need them.

## Workflow

One page contains the preview, global trigger offset, output integration and Export .vit button. Load a raw folder or ZIP in the sidebar. Choose an acquisition to inspect, then adjust **Trigger time offset (s)**: positive shifts laser triggers later, negative earlier. This offset is added to each acquisition's saved transit correction and applied both to the preview and every exported spot window. Signal timestamps are not moved.

Choose the output integration time. Enable **Bin counts into an X/Y pixel grid** to set independent pixel width and height (default 5 × 5 µm), with optional grid origin. Export all loaded acquisitions by default, or only the previewed acquisition. Click Export .vit and download the completed package. The supplied v7 template is built in. Mapping and diagnostics are under a collapsed Details & advanced settings area rather than separate tabs.

Stationary spot names and XYZ come from laser.info. Trigger trains are associated in recorded order when their count matches the spot count. When counts differ, an association is accepted only if ordered shot counts (or numeric shot-count suffixes in names when metadata counts are zero) give one unique match. This fallback is explicitly reported as inferred. Missing/ambiguous spots stay in unassigned data, without guessed coordinates or discarded counts. The supplied dataset has absent single-shot triggers; this version cannot automatically name those signal intervals. Trigger windows run from the first recorded pulse through one repetition period after the last; absent first pulses and washout tails are not reconstructed.

No manual interval/coordinate table is required in the main workflow. Source files remain read-only. Isotope mapping and CPS remain provisional as described below.

## Exact v7 layout

`v7_template.csv` preserves the 14 metadata lines and the original 275-column header from the supplied `line_0_20260702161619286342.csv`. A replacement v7 template can be uploaded. Exports preserve metadata key spelling, order, whitespace prefixes, time/X/Y column order, isotope column labels and isotope column order. No extra columns or metadata lines are added to the iolite CSVs.

Each named spot interval becomes a separate CSV. `Laser line name` carries the interval name into both iolite's data Samples and laser-log Samples, linked to its start/end time. Background intervals are exported separately with a `Background -` or `Unassigned background -` prefix. Names and intervals are derived from the raw trigger association described above; unresolved times are explicitly unassigned. In iolite, open Data → Samples and use the named samples for Create selection from sample or Auto selections. Filter to the desired spot names to exclude background entries, and use either data samples or log samples to avoid duplicate selections.

`Timestamp` is converted from `run.info` to the chosen timezone and written without an offset, matching v7. `Cycle time (ms)` remains elapsed time from the acquisition origin. Starting X/Y/Z comes from the first output row, or `nan` if unassigned. Coordinates remain fixed within mapped intervals. Each CSV ends with an all-NaN signal boundary row at the interval's actual end: this preserves sample duration even when a spot is summed into one measurement, without duplicating counts. The direction line is retained for compatibility. Spot metadata uses the matching recorded spot name where available, falling back to the first candidate; review metadata for manually renamed spots.

The `.vit` archive contains only root-level `line_*.csv` files. It does not contain manifests, count columns, additional folders or diagnostic CSVs. Diagnostics are in a separate ZIP.

## Spatial pixels and counts

The spatial grid uses `floor((x-origin_x)/width)` and `floor((y-origin_y)/height)`. Cells include their lower edge and exclude their upper edge. It combines all linked signal contributions in the same cell across the selected acquisitions, including different spots and repeat visits. Rows use the grid-cell center; source coordinates and contributing spot names are retained in `pixel_contributions.csv`. Grid assignment uses the recorded point position, not fractional overlap of the laser footprint. Stationary holes cannot yield finer spatial detail than their recorded positions. The current raw reader supports stationary spots; moving raster coordinates are not yet decoded.

Counts are stored signal divided by AverageSingleIonArea. Spatial bins sum counts and actual contributing exposure, excluding intervening gaps, then write counts/exposure into the required v7 **CPS** columns. Exact counts for every stored channel are in `pixels_counts.csv`. Unmapped isotope columns remain NaN. These estimates do not reproduce NuQuant quantification.

Spatial mode produces three downloads:

- `vitesse_spatial.vit`: one CSV per occupied pixel, with one measurement row plus a NaN endpoint marker. Its time column is a **synthetic cumulative exposure clock**, not original measurement time, because multiple acquisition times may contribute to one pixel. Shared metadata are retained; conflicting spot size/rate/fluence values are NaN. Multiple contributing names are joined.
- `vitesse_export.vit`: the complete time-resolved companion, including all background and unassigned data at the chosen temporal integration. Original acquisition times are preserved. Import spatial and timeline archives into **separate iolite sessions** to avoid counting the same signals twice.
- `vitesse_spatial_details.zip`: exact counts, exposures, source acquisition/time/spot/XYZ contributions, grid dimensions/origin and summary-clock explanation. Timeline diagnostics remain available separately under Export details.

Unknown-coordinate time is not placed into a fabricated spatial pixel. It remains in the complete timeline, including signal from single-shot spots with missing triggers. Spatial export stops if no linked positions exist. Grid resolution does not alter the time-resolved companion or split an individual stationary hole into artificial pixels.

## Validation and current limits

- Automated tests cover exact v7 metadata-key/header ordering, full column counts, ZIP structure, background retention, count conservation, partial bins, interval overlaps, missing chunks, safe upload extraction and separation of disjoint pixels.
- Tested interactively through Streamlit's AppTest with the first supplied acquisition: loading 90,651 records, changing bin interval, assigning coordinates, pixel summation and archive download preparation.
- Reader format and basic scaling were checked against [pewlib](https://github.com/djdt/pewlib/blob/master/src/pewlib/io/nu.py). The first acquisition's decoding agrees with that implementation.
- **Named-spot import verified in iolite 4.10.12:** two test intervals imported separately as `915001-1` and `915001-10`, alongside three background intervals. Both data and laser-log Samples had matching names and boundaries (0.800 s and 1.499 s durations, within millisecond metadata precision). These were test assignments, not validated laser timing identities.
- **Earlier single-acquisition import verified in iolite 4.10.12** in fresh test sessions. Time-resolved export: 4,773 CSV rows loaded as 4,774 iolite points, with one trailing NaN separator. Hf values, relative times (within 0.5 microseconds), X and Y matched the CSV. Pixel-sum export: 4,694 CSV rows loaded as 4,695 points, also with one trailing NaN separator. All 272 isotope columns, relative times, X and Y matched, including missing values. Test reports are in `import_validation.json`. This verifies file transport and time-series values, not isotope quantification or all imaging/reduction workflows.
- **Your iolite default was left unchanged as requested.** Preferences → Import → Vitesse → “Treat channels as already background subtracted” is currently enabled. The tests confirmed that iolite sets `BackgroundSubtracted=True` even though the exported background remains in the actual signal. This can bypass iolite's baseline-selection requirement. The unchanged v7 CSV format does not convey that flag. Review it for this workflow before doing background-dependent data reduction; the app does not change iolite preferences.
- Iolite also warns that these files do not specify channel dwell-time metadata. The CSV preserves actual time coordinates, while exact exposure durations live in the diagnostics ZIP. A single global dwell setting does not represent variable-duration pixel rows. Import success does not validate dwell-dependent uncertainty calculations.
- Use recorded X/Y for imaging; rows with unassigned coordinates cannot locate a pixel. They remain available for time-series/background work rather than being assigned fabricated coordinates.
- **This is a working extraction/export prototype, not a validated replacement for NuQuant quantification.** The original dataset stores 196 channels but its NuQuant CSV has 272 isotope columns. Channel mapping suggestions cannot establish isotope identity or correct interferences. The preliminary silicon comparison differs from NuQuant by about 27%; Hf/U examples were closer. No extra background subtraction, isotope deconvolution, dead-time, detector, or auto-blanking correction is applied. Stored signals may already include instrument-side processing.
- The supported reader accepts continuous single-cycle/single-segment files with fixed channel centers, including gzip-compressed `.integ` files. Missing files, counter gaps, changed centers or unsupported segments cause a clear error instead of a silent partial export. Batch export requires consistent mass channels across acquisitions.
- Every stored integration is retained, but all available mass channels are only present in the diagnostics. The v7 archive follows the reviewed mapping to the fixed isotope columns. Original files remain the authoritative unmodified source.

Run core tests with `python -m pip install pytest` followed by `python -m pytest test_core.py -q`.
