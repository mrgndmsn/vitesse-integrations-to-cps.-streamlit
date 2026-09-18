# Vitesse — GitHub / Streamlit deployment

This is the updated one-page app with trigger offset and spatial grid controls.

## Deploy

Put the contents of this folder at your repository root, including `.streamlit/config.toml`. In Streamlit Community Cloud, select the repository and branch and set the entrypoint to `app.py`. Choose Python **3.12** in Advanced settings. `requirements.txt` belongs next to `app.py`. Do not upload the ZIP itself as the app, and do not use `run.command` as the entrypoint.

Required runtime files: app.py, core.py, spatial.py, v7_template.csv and requirements.txt. Keep these together. The package versions match the tested local runtime; dependency installation on the deployed host still needs verification in its logs.

## Data

Use Upload ZIP, preserving Image001/00001/etc subdirectories and their run.info, integrated.index, .integ, pulse.index, .pulse, laser.info and TriggerCorrections.dat files. Local folder refers to the server filesystem, not your Mac. The configured upload cap is 1 GB per file; this is not a promise that the host has enough memory to process that amount. Large runs are better processed locally. Uploaded raw data are sent to the hosting server, processed in temporary directories, and are not automatically committed to GitHub. Temporary storage may last until server cleanup/restart.

## Outputs and limits

The complete timeline .vit keeps background and original timing. Optional spatial .vit aggregates by X/Y cell and uses a synthetic exposure clock; import these into separate iolite sessions. Each keeps the supplied v7 layout. Exact counts and source contributions are in diagnostics. CPS and nominal-mass assignments are provisional and do not reproduce NuQuant corrections. The current reader supports stationary spot positions, not moving raster coordinates. Missing triggers cannot be assigned fabricated coordinates. No iolite importer files or preferences are changed.

## Troubleshooting

Check Manage app → logs. A dependency installation error requires checking Python version and the reported package. A ModuleNotFoundError for spatial or core means a runtime file is missing. A file-size rejection means `.streamlit/config.toml` was not included or the file exceeds the configured cap. An out-of-memory failure requires a smaller acquisition or a larger host/local execution. Share the exact log to diagnose other errors.

https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/file-organization
