# QC5 effective-gain analysis

`QC5graph.C` is a CERN ROOT macro for QC5 measurements of ME0 GE1/1 detector sectors. For each configured sector, it calculates an effective gas gain from the measured current and a plateau rate, then produces a dual-axis plot of gain and trigger rate against divider current.

## Requirements

- CERN [ROOT](https://root.cern/) 6 or a compatible release
- A CSV file for each sector: `sector3_2.csv`, `sector4_2.csv`, and `sector5_2.csv`

## Run

Place the macro and the three input CSV files in the same working directory, then run:

```sh
root -l -q 'QC5graph.C()'
```

The macro reports the selected plateau rows and rate to the terminal. Change the three `DrawSector(...)` calls at the bottom of `QC5graph.C` to analyze different files or labels.

## Input CSV format

The first row is treated as a header and skipped. The macro expects comma-separated numeric values, with at least ten columns. It reads these **zero-based** column positions:

| Column position | Value used |
| ---: | --- |
| 2 (third column) | Divider current, `Imon` [µA] |
| 4 (fifth column) | Background/scaler-off width, `XOFF_WDT` |
| 6 (seventh column) | Scaler-on width, `XON_WDT` |
| 9 (tenth column) | Measured current, `IR` [A] |

For each valid data row, the rate is calculated as `(XON_WDT - XOFF_WDT) / 60`. Rows with a non-positive rate, fewer than ten columns, or non-numeric values in the required fields are skipped with a warning.

## Outputs

For each input sector, ROOT writes:

- `Sector 3.png` and `Sector 3.pdf`
- `Sector 4.png` and `Sector 4.pdf`
- `Sector 5.png` and `Sector 5.pdf`

Files are created in the directory from which ROOT is run. The output filenames follow the sector labels supplied to `DrawSector`.

## Analysis notes and caveats

- The plateau is the consecutive three-row window with the smallest total change in measured rate. With fewer than three valid rows, all available valid rows are used.
- Effective gain is calculated as `abs(IR) / (e × N_primary × plateau_rate)`, using `e = 1.6 × 10⁻¹⁹ C` and `N_primary = 346`.
- The exponential fit is `G = A exp(α I)`, where `I` is divider current in µA. The fit and smoothing splines are visual aids; assess their physics validity and uncertainties before using them in a final result.
- The CSV reader is intentionally simple: it does not support quoted fields containing commas. Confirm the fixed column positions and the 60-second normalization match the data-acquisition export before analysis.

