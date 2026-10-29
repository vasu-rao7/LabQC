# LabQC 384-well plate

Card submitted to the lab automation team on 2026-10-12. Attachments: glucose_qc384_well_map.csv, glucose_qc384_plate_grid.csv and glucose_qc384_glucose_levels.csv from this folder. The card notes below are as submitted.

---

(Written with Claude)

Design

This is a single 384-well plate to measure plate noise and spatial error (QC). The one factor is initial glucose concentration. I'm treating it as a 96-well experiment where each "well" is a 2x2 block of 384 wells in the same spot as that well on a 96-well plate. For example, 96-well A1 is 384-well A1, A2, B1 and B2.

Every block has the same layout:

- Top left: positive control (10 g/L glucose)
- Bottom left: negative control (no glucose)
- Top right and bottom right: the block's condition, as 2 replicate wells

Organism. Streptococcus mutans UA159 (wild type). All 384 wells are inoculated, including the negative controls.

The plate will contain 96 blocks:

- 12 positive control blocks (both right-hand wells at 10 g/L)
- 12 negative control blocks (both right-hand wells with no glucose)
- 72 test blocks: 24 glucose levels x 3 blocks each

That works out to 120 positive control wells, 120 negative control wells and 144 test wells (6 per glucose level). There are no empty wells.

Factors

Glucose. 26 levels on a linear scale from 0 to 10 g/L in 0.4 g/L steps: 0 g/L (negative control), 0.4, 0.8 ... 9.6 g/L (the 24 test levels) and 10 g/L (positive control). 10 g/L (1% w/v) is the standard CDM level, so the positive control is just standard CDM.

Everything else in the medium (amino acids, vitamins, minerals, bases, buffer) is at standard CDM levels in every well. No stressors. Since only glucose changes, the test levels can be made by mixing standard CDM with glucose-free CDM (for example, 4 g/L = 40% standard CDM + 60% glucose-free CDM), but use whatever works best for you.

Layout

Please use the layout in the attached files as given. Where each block sits is part of the QC design. The 24 control blocks are placed with a Latin hypercube design on the 96-well block grid (every block column has 1 positive and 1 negative control block). The 72 test blocks are randomized, with the 3 blocks for each glucose level in different block rows and columns.

Methods

Inoculate the plate and take OD600 kinetic reads every ~15-20 minutes for 24 hours. Aerobic (5% CO2), standard 37 degrees C, shake the plate before each read.

- glucose_qc384_well_map.csv -- 384 rows, one per well. Each row gives the well's 96-well block, its position in the block (TL/TR/BL/BR), the block type, the well type and the glucose concentration (g/L).
- glucose_qc384_plate_grid.csv -- the same layout as a 16x24 grid of glucose concentrations (g/L), for a quick visual check.
- glucose_qc384_glucose_levels.csv -- 26 rows. Each glucose level and the number of wells that get it.

NOTE: Please send back the OD600 data for all 384 wells, including every control well. The per-block controls are the main point of this plate.
