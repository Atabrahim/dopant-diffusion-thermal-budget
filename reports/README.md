# Synthetic process reports

Run `python examples/reproduce.py` to regenerate the tracked files.

Each case contains a JSON report with unit-labelled inputs, thermal budget, dose-balance diagnostics and every threshold crossing, plus a CSV of the initial and final chemical concentration in cm^-3 against depth in nm. `process_comparison.json` groups the two declared anneal cases and six sensitivity rows; `sensitivity.csv` provides the same summary in tabular form.

These files are deterministic numerical simulations from the declared model. They are neither experimental measurements nor a calibrated prediction for a specific fabrication line. The deepest threshold crossing is a descriptive model diagnostic, not an electrical junction depth.
