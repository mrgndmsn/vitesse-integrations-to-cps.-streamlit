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

1. Load a results folder, upload the original directory, or upload its ZIP. Preserve subdirectories so identically named files from different acquisitions remain separate. The reader needs `run.info`, `integrated.index` and `.integ` files; include `.pulse`, `pulse.index`, `laser.info` and `TriggerCorrections.dat` for timing and coordinate metadata. A `.raw` waveform file alone is not supported.
2. Inspect a selected acquisition. Choose channels, zoom the trace, and optionally overlay native integrations. Change the target integration interval in the sidebar. The preview and exports use the same binning implementation.
3. Add signal/background time intervals and their fixed X/Y/Z. Copy positions from recorded spot metadata or enter them manually. Save edits with Apply interval edits. Tables can also be uploaded/downloaded as CSV. The app never guesses a spot identity from an incomplete pulse train.
4. Choose regular time bins or a single pixel per mapped signal interval. Background and unassigned intervals remain time resolved. No time interval is discarded or duplicated. Intervals snap to native integration midpoints; diagnostics expose the actual resulting boundaries.
5. Review the channel mapping to your template isotope columns. Suggestions are based on nominal mass proximity only. Apply changes. Unmapped columns remain `nan`, not zero. Multiple isotope labels can reference one channel, but this does not separate isobars or reproduce NuQuant's isotope processing.
6. Select one or more acquisitions, review the export choices, and build the package. Download `.vit` directly, the identical bytes named `.zip`, and the separate diagnostics ZIP.

## Exact v7 layout

`v7_template.csv` preserves the 14 metadata lines and the original 275-column header from the supplied `line_0_20260702161619286342.csv`. A replacement v7 template can be uploaded. Exports preserve metadata key spelling, order, whitespace prefixes, time/X/Y column order, isotope column labels and isotope column order. No extra columns or metadata lines are added to the iolite CSVs.

Each named spot interval becomes a separate CSV. `Laser line name` carries the interval name into both iolite's data Samples and laser-log Samples, linked to its start/end time. Background intervals are exported separately with a `Background -` or `Unassigned background -` prefix. Assign a name and time interval before exporting; the app does not guess spot timings. In iolite, open Data → Samples and use the named samples for Create selection from sample or Auto selections. Filter to the desired spot names to exclude background entries, and use either data samples or log samples to avoid duplicate selections.

`Timestamp` is converted from `run.info` to the chosen timezone and written without an offset, matching v7. `Cycle time (ms)` remains elapsed time from the acquisition origin. Starting X/Y/Z comes from the first output row, or `nan` if unassigned. Coordinates remain fixed within mapped intervals. Each CSV ends with an all-NaN signal boundary row at the interval's actual end: this preserves sample duration even when a spot is summed into one measurement, without duplicating counts. The direction line is retained for compatibility. Spot metadata uses the matching recorded spot name where available, falling back to the first candidate; review metadata for manually renamed spots.

The `.vit` archive contains only root-level `line_*.csv` files. It does not contain manifests, count columns, additional folders or diagnostic CSVs. Diagnostics are in a separate ZIP.

## Counts versus CPS

The reader uses the stored signal divided by `AverageSingleIonArea` to estimate counts. It divides by measured exposure to get estimated CPS. Binning sums counts and divides by total exposure; it never sums CPS values directly. Whole native integrations are used, so requesting 10 ms on the supplied dataset produces 10.1194 ms bins (19 × 0.5326 ms), plus shorter bins at boundaries. Pixel mode sums a complete mapped signal interval into one row while leaving other time intervals intact. Identical X/Y in two disjoint intervals does not cause the intervening background to disappear.

Each diagnostics `group_N_counts.csv` records interval name/type, actual start/end, exposure, native integration count, X/Y/Z and summed mass-channel counts. The manifest records requested/effective binning, missing coordinates, source quality flags and scientific limitations. Channel mapping and interval tables are included there too.

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
