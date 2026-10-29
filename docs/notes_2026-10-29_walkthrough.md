# Walkthrough notes, 2026-10-29

Questions that came up while reading the results page, with the answers. Items marked **(new)** were checked against the code or data during the walkthrough and are not on the page.

## Measures and methods

**Ceiling.** Half the mean squared difference between side-by-side wells is treated as unpredictable scatter. The ceiling is one minus that scatter divided by the plate's variance, and the dashed line is its median over a panel's surfaces (`variance_decomposition` in `src/experiments/metrics.jl`). It runs low, because some real smooth change between neighbours gets counted as scatter. It is a rough guide, not a limit.

**Row + column fit.** There is one offset per row and per column (20 unknowns on a 96-well plate). Each control gives one equation: its error = its row offset + its column offset. These are solved by least squares with the pseudo-inverse, which picks the smallest-norm solution when there are fewer controls than unknowns, so a row or column with no control gets 0 (`predict_rowcol` in `src/experiments/baselines.jl`). Around 16–20 controls there are about as many equations as unknowns. The model then fits every control exactly, noise included, and its predictions blow up.

**Inverse-distance vs GP.** Inverse-distance isn't a GP kernel, but both predict a weighted average of the control errors. The differences:
- **Weights:** inverse-distance always uses 1/distance². The GP learns from the controls how far the correlation reaches and how noisy they are.
- **Noise:** inverse-distance passes exactly through every control, so it copies control noise to the neighbours. The GP's noise term smooths over it.
- **Clustered controls:** the GP accounts for controls being near each other; inverse-distance double-counts them.
- **Far from controls:** the GP goes back to 0 (the control median); inverse-distance keeps averaging all the controls.
- **Uncertainty:** only the GP gives a ± range.

With few controls the GP can't learn its settings well (most fits hit the 2-well floor), which is why inverse-distance wins at low counts.

**`subtract_and_normalize`.** It has been in the package since the first commit (`src/correction/correction.jl`): subtract the negative-control mean, then divide by the positive-minus-negative gap. The all-control sweep never ran it, because it needs both control types on one plate. It appears only in the poster study, as a reference row. It scores exactly 0, which shows that normalization alone does nothing about the spatial error and that replicate skill ignores shifts and rescaling. Its numbers match the uncorrected row, so it can be dropped from the paper if that point isn't needed. (An older note said it "was never run"; that was true before the poster study.)

## Placement figure (`fig_exp1_placement`)

- **Dots:** 30 per setting in every panel, one per layout, so 240 per panel (verified for all 96 setting × panel combinations). Each dot is the layout's GP skill averaged over the plate's surfaces. At a given control count, all four plate rows use the same layouts.
- **Black bar:** the mean of all 30, including dots below the clipped axis (the corner note counts those).
- **Settings:** random; LHS (each row and column gets its fair share, distances ignored); minisum* (minimizes the summed Manhattan distance from every well to its nearest control, named "minimax" in PlateArrays); hybrid λ = λ × spread + (1 − λ) × LHS, each rescaled to 0–1.
- **Both classes placed:** each layout places n positives and n negatives on the 96-well plate, and the GP fits one class. At n = 40 that fills 80 of 96 wells, so the layouts can't differ much.
- **(new) Two columns are duplicates.** Hybrid λ = 0 optimizes the same objective as LHS, and λ = 1 the same as minisum*. Only the random seeds differ. At n = 12 their median layout scores match exactly (LHS score 4 vs 4, spread score 254 vs 254, from `results/layouts.csv`). The gap between the two columns of a pair therefore shows how far a mean moves by chance: up to about 0.04–0.07 (negative n = 12: −0.02 vs −0.06; gas relative n = 12: 0.12 vs 0.19), and 0.51 vs 0.38 when a dead-well outlier lands in one batch. Most gaps between other settings are smaller than that.
- **(new) Spreading out, not row/column balance, narrows the scatter.** Standard deviation of skill across layouts at n = 12:

  | Plate | random | LHS | minisum* |
  |---|---|---|---|
  | negative | 0.13 | 0.10 | 0.04 |
  | positive | 0.07 | 0.17 | 0.05 |
  | gas relative | 0.19 | 0.26 | 0.06 |
  | gas additive | 0.13 | 0.10 | 0.08 |

- **(new) PlateArrays `:labmda` typo is harmless.** `exchange` reads λ under a misspelled key and defaults to 0.5, but the user's λ still reaches `hybrid` through the splatted keyword arguments, which override the default. The layout scores follow λ as they should.

## Uncertainty vs error figure (`fig_exp2_sd_vs_error`)

- **What a dot is:** one layout (240 per panel). Across: the GP's mean predictive sd over held-out wells, noise included. Up: actual RMSE at those wells. ρ is their rank correlation.
- **Answer:** no, the GP doesn't know which layouts are good. On the negative and positive plates the dots fall (ρ −0.38 to −0.62), and on gas additive there's no pattern (about −0.1).
- **(new) Mechanism, checked in `results/sweep_*.csv`.** At 12 controls, the predictive sd tracks the fitted noise (rank correlation 0.91 negative, 0.83 positive), and the fitted noise correlates negatively with error (−0.62, −0.55). Fits that call little of the variation noise get narrow ranges and chase the wiggles in the controls (75% and 63% sit at the 2-well floor), so they predict worse. The page's range for the sd–noise link (0.78–0.90) differs slightly from this recomputation per layout (0.83–0.91); the story is the same.
- **(new) At 40 controls the noise explanation weakens:** the noise–error correlation drops to about −0.1, yet ρ stays negative (−0.38, −0.55). There is no full explanation for this yet.
- **Gas relative (−0.73):** the page attributes it to the dead wells. The mechanism wasn't checked.

## Centre vs edge figure (`fig_e4_rings`)

- **What it shows:** where the error is, not how well each method does at each place. Up is each middle well minus its glucose level's median, averaged by ring (1 = edge, 8 = centre), in control-window units. The grey uncorrected line is the artifact itself.
- **Why the GP beats inverse-distance most at the centre:** that's where the error is, and the error there is a local peak. Inverse-distance averages all the controls, and the many off-centre ones flatten the peak. The GP's fitted lengthscale (2–6 wells) keeps its predictions local. At the edge there is little to remove, so the methods converge, though the GP is still slightly better (ring 1 at read 90: −0.015 vs −0.05).
- **(new) Controls reach the centre.** The poster layout puts 8 positives and 8 negatives in rings 5–8 (none of the 24 positives are in ring 8), so the GP sees the bump directly.
- **(new) The centre reads high; the edge doesn't read low.** No-glucose controls at read 90 average 0.095 OD in ring 1 (blank-like) and 0.232 OD in rings 5 and in. The rings figure alone can't tell which, because two-thirds of the scored wells are in rings 1–4.
- Nearest control and row + column aren't on this figure. The before/after table shows them leaving more of the centre (centre excess 0.17 and 0.19, against 0.06 for the GP).
- The physical cause of the central block is not established. It shows up on all four 2021 plates.
