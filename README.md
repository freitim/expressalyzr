# expressalyzr

`expressalyzr` is an experiment-level workflow wrapper for flow cytometry
analysis in R. It combines the Bioconductor flow-cytometry stack
(`flowCore`, `flowWorkspace`, `openCyto`, `ggcyto`) with a thin orchestration
layer, declarative configurations, and self-contained experiment directories.

## Install

```r
# install devtools once if needed
install.packages("devtools")

devtools::install_github("Ariwor/expressalyzr")
```

## What An Experiment Folder Looks Like

The only thing you must have at the start is a folder with FCS files.

```text
experiment_001/
  01-Well-A1.fcs
  01-Well-A2.fcs
  01-Well-A3.fcs
```

It is recommended that you also include a design CSV file (e.g.`Design_experiment_001.csv`) where you include metadata and relevant info about your experiment.

```text
experiment_001/
  01-Well-A1.fcs
  01-Well-A2.fcs
  01-Well-A3.fcs
  Design_experiment_001.csv
```

## First Run

For a new experiment, you can start with the simplest guided-review run:

```r
library(expressalyzr)

results <- run_pipeline("experiment_001")
```

What happens:

1. FCS files are moved into `data/` if needed, then loaded.
2. Configuration files (`config.yml` and `gt_samples.csv`) are copied from a
   starter profile if they do not already exist. `config.yml` describes
   experiment-level settings such as beads, controls, manual prompts, and
   design-file matching. `gt_samples.csv` is the OpenCyto gating template.
3. The config opens for review.
4. The visual gate designer opens.
5. Gates, compensation, positivity calls, and optional QC/report sections run.
6. Results are written back into the experiment folder.

`results` is an event-level `data.table`. The same table is also written to `experiment_001/experiment_001.csv`.

After the first run, you can turn individual pieces on or off depending on
whether you are still reviewing the experiment or running a finished analysis
script.

## Interactive Vs Non-Interactive Runs

`expressalyzr` supports a range of interactive modes, allowing you to control
how much of the workflow you manage directly.

`interactive = TRUE` means `expressalyzr` is allowed to pause
and ask the user to do something:

- open `config.yml` for review
- launch the visual gate designer
- ask for manual compensation and cutoff adjustments
- show manual review plots during those prompts

Example usage:

```r
results <- run_pipeline(
  data_path = "experiment_001",
  interactive = TRUE,
  report = TRUE
)
```

`interactive = FALSE` means unattended/script mode. `expressalyzr` will not open
editors, launch browser apps, call console prompts, or create required setup
files silently. Use it after the experiment folder is already configured.

Example usage:

```r
results <- run_pipeline(
  data_path = "experiment_001",
  interactive = FALSE,
  gate_designer = FALSE,
  report = TRUE
)
```

You can also control other settings when calling `run_pipeline()`:

- `view_config = TRUE`: opens `config.yml` for review in interactive runs.
- `interactive = TRUE`: fully interactive/review mode is allowed
- `gate_designer = TRUE`: the visual gate designer launches before gating
- `gating_templates = NULL`: use the standard single-template
  `gt_samples.csv` workflow.
- `gated_populations = NULL`: only needed when template labels need to be
  supplied separately.
- `report = TRUE`: write the HTML report without changing the returned
  `data.table`.
- `include_gating_plots = TRUE`: include full per-well gate plots in the HTML
  report for the default guided-review workflow.
- `summaries = FALSE`: do not write extra population/marker summary CSVs unless
  requested.
- `save_gating_set = FALSE`: do not save `flowWorkspace` GatingSet objects
  unless requested.
- `config_profile = "default"`: when a new config is created, use the
  recommended starter profile.
- `gating_output = NULL`: legacy compatibility slot; normally leave it unset.

To create a different starter config, use:

```r
run_pipeline(
  "experiment_001",
  view_config = TRUE,
  interactive = TRUE,
  config_profile = "default"
)
```

Available profiles:

- `default`: recommended for new experiments; includes manual compensation and
  cutoff review during interactive runs.
- `legacy`: closer to older expressalyzr defaults.
- `advanced`: exposes manual/debug settings.

## Running From The Terminal

`expressalyzr` does not require a separate command-line program. You can run the
same R workflow from a terminal with `Rscript`.

For a first guided review:

```sh
Rscript -e 'expressalyzr::run_pipeline("experiment_001")'
```

For a finished folder that should run unattended:

```sh
Rscript -e 'expressalyzr::run_pipeline("experiment_001", interactive = FALSE, gate_designer = FALSE)'
```

## Outputs

`run_pipeline()` writes outputs inside the experiment folder so the analysis can
be reviewed later and shared with the original FCS files.

### Main Event Table

The main CSV is always written:

```text
experiment_001/experiment_001.csv
```

This table contains one row per gated event/cell in the final extracted population. Typical columns include:

- sample/file identifiers such as `File`
- metadata merged from the design file, if one is configured
- processed fluorescence/scatter channel values such as `FSC-A`, `SSC-A`, and
  `FL*` channels
- gating labels such as `gating_template` and `gated_population` when using
  multi-template/co-culture mode
- positivity helper columns such as `positive`, `no_negative`, and channel
  positivity calls when positivity thresholds are available

In R, the same table is available as the return value as a `data.table`.

### HTML Review Report

With the default `report = TRUE`, the pipeline writes:

```text
experiment_001/experiment_001_report.html
```

You can open it in a browser to check:

- sample/event-count QC
- compensation status and spillover heatmap
- gate hierarchy and per-well gate plots
- channel distributions
- population and marker summaries
- files written by the run

### Optional Summary CSVs

For plate summaries, figures, or statistical analysis, add `summaries = TRUE`:

```r
results <- run_pipeline("experiment_001", summaries = TRUE)
```

This writes:

```text
experiment_001/experiment_001_populations.csv
experiment_001/experiment_001_markers.csv
```

`*_populations.csv` has one row per sample and gated population, with event
counts and percentages relative to the parent and total population.

`*_markers.csv` has one row per sample, population, and channel, with median,
mean, spread, min/max, and percent-positive summaries.

### Other Optional Outputs

- `save_gating_set = TRUE` saves reusable `flowWorkspace` GatingSet output in
  `gs/` or `gs_<template>/`.
- `include_gating_plots = FALSE` keeps routine reports smaller by omitting full
  per-well gate plots.
- When compensation controls are matched from the design file, the pipeline
  writes `<experiment>_compensation_controls.csv` so you can review which wells
  were used for each channel.

## Gating

Most experiments use one gating template:

```text
gt_samples.csv
```

The visual designer can edit this file before the pipeline gates all samples:

```r
run_pipeline(
  "experiment_001",
  view_config = FALSE,
  interactive = TRUE,
  gate_designer = TRUE
)
```

In the designer, you can:

- Select a sample/well to preview.
- Draw or adjust boundary and polygon gates.
- Edit gate parameters in the table.
- Switch axes and axis scales.
- Filter the sample dropdown to control wells.
- Save the edited template back to `gt_samples.csv`.

Automated openCyto gates such as `density_gate`, `mixture_gate`, and
`singletGate` remain template rows and can be edited.

## Multiple Cell Populations

For co-culture experiments, different populations may have different scatter
profiles and therefore need different gates. Use named gating templates:

```r
result <- run_pipeline(
  "experiment_001",
  view_config = FALSE,
  interactive = TRUE,
  gate_designer = TRUE,
  gating_templates = c(
    HEK = "gt_HEK.csv",
    Tcells = "gt_Tcells.csv"
  ),
  report = TRUE
)
```

If `gt_HEK.csv` or `gt_Tcells.csv` do not exist yet, an interactive run creates
starter copies from the standard `gt_samples.csv` template and opens them in the
combined gate designer. Edit and save each population template there. During the
run, each template gates the same FCS files independently, and the extracted
events are combined into one returned `data.table`.

The output includes:

- `gating_template`: which template was used.
- `gated_population`: the biological label, such as `HEK` or `Tcells`.

Population labels also let compensation and background thresholds stay separate
for co-culture runs when your controls are population-specific. Add a
`gated_population` column to the control rows in the design CSV and use the same
labels as the named templates. Only FCS files matching `controls_pattern` are
treated as controls; the `gated_population` column tells `expressalyzr` which
population should use each matched control row.

You do not need `gating_templates` for normal single-population experiments.

## Compensation Controls

For new experiments, the clearest workflow is to label single-color controls in
the design CSV. The column is usually named `Channel`.

Example design rows:

```csv
Well,Sample,Channel
A1,empty_vector,unstained
A2,mGL,FITC
A3,miRFP670,APC
```

The config tells `expressalyzr` which FCS files are controls and which design
column contains the control label:

```yaml
default:
  controls_pattern: "Controls"
  controls_design_column: "Channel"
  controls_population_column: "gated_population"
  channel_map: "cytoflex"
```

`expressalyzr` automatically matches labels such as `FITC` or `APC` to the detector channels in
the FCS metadata, orders the controls by machine channel order, and writes a
review table:

```text
experiment_001/experiment_001_compensation_controls.csv
```

For multi-population/co-culture experiments, controls can be scoped to one
population:

```csv
Well,Sample,Channel,gated_population
A1,HEK_unstained,unstained,HEK
A2,HEK_FITC,FITC,HEK
B1,Tcell_unstained,unstained,Tcells
B2,Tcell_FITC,FITC,Tcells
```

When population-specific control rows exist, `expressalyzr` uses them for that
population's compensation matrix and background/positivity thresholds. Blank
`gated_population` rows remain global fallback controls for older or shared
control designs.

Older experiments and scripts can still use `controls_index`.
