# Versioning

Last updated 2026-10-29.

## Remotes

| Remote | URL | Access |
|---|---|---|
| `origin` | `git@github.com:vasu-rao7/LabQC.git` (Vasu's fork, public) | push and pull |
| `upstream` | `https://github.com/jensenlab/LabQC.git` (Ben David's original, public) | pull only |

The vasu-rao7 account has no write access to `jensenlab/LabQC`. Every commit there is Ben's, and pushes over HTTPS and SSH both return 403. Work goes to the fork. To bring something into Ben's repo, push it to the fork and open a pull request against `jensenlab/LabQC`.

Both machines push and pull through `origin` with the same SSH key (`ssh -T git@github.com` answers "Hi vasu-rao7!").

## Machines

| Machine | Path | Usual branch |
|---|---|---|
| goblin (lab server, over SSH) | `~/LabQC` | `sweep-rerun-poster-e4` |
| Mac | `~/Documents/GitHub/LabQC` | `main` |

## Branches on the fork

- **`main`**: Ben's history plus `0b71894` "Add glucose QC 384-well plate design" (`designs/glucose_qc384/`, made on the Mac).
- **`sweep-rerun-poster-e4`**: all the GP-correction experiment work. It contains `main` (merged in as `9bf7a7e`), so it is a superset of `main`. It has not been merged back into `main`.

| Commit | What it does |
|---|---|
| `6755b12` | Experiment harness and sweep code as run on 2026-09-04 |
| `9bf7a7e` | Merge of `origin/main` (the glucose QC 384 design) |
| `7afd413` | `.gitignore` for the raw sweep tables (see below) |
| `2d4e520` | Manifest sync (Serialization) |
| `62da341` | 384-well kinetic reader, elapsed-time text; poster E4 input files |
| `72d4671` | Per-axis (ARD) kernel option for `plateerrorGP`; golden test for the default kernel |
| `a5c6d85` | Sweep fixes: baselines on the right scale, ceiling on the right scale, predictive-sd coverage, screened scores, run records |
| `e9a78cf` | 2026-10-09 sweep outputs (tables, layouts, run record) and the 2026-09-04 archive |
| `84df312` | Poster plate E4 study: code, tests, results |
| `449fa3c` | Sweep figures redrawn: gas additive panels, linear control-count axes |
| `62042b4` | Full-plate run (3–96 controls): code, summary, figures |
| later | These notes |

The test suite (`julia --project=test test/runtests.jl`) passed at `62042b4` with 138/138.

## Common tasks

```bash
# get Ben's new commits
git fetch upstream && git merge upstream/main

# Mac: get goblin's experiment work
git fetch origin && git switch sweep-rerun-poster-e4   # creates a tracking branch on first use

# goblin: get anything new pushed from the Mac to main
git fetch origin && git merge origin/main

# after committing
git push            # both branches already track origin
```

## What git leaves out

**Raw per-cell sweep tables.** Ignored in `.gitignore` because each is 43–260 MB and GitHub rejects files over 100 MB. They exist only on goblin:

| File | Size |
|---|---|
| `results/sweep_negative.csv` | 187 MB |
| `results/sweep_positive.csv` | 259 MB |
| `results/sweep_gas.csv` | 50 MB |
| `results/full_n/full_n.csv` | 105 MB |
| `results/archive_2026-09-04/sweep_{negative,positive,gas}.csv` | 162, 223, 43 MB |

The tracked summary tables (`summary_*.csv`, `per_layout_gp.csv`, `summary_full_n.csv`) are built from them, so every figure can be redrawn without them. Regenerating them means rerunning `experiments/run_sweep.jl` and `experiments/run_full_n.jl` (the full-plate run took 26 min on 112 threads). If they need to be shared or backed up, use Git LFS or copy them off goblin.

**Loose images in `raw_data/`.** Ignored since 2026-10-29 (`raw_data/*.png`, `raw_data/**/*.jp*g`); they exist only on goblin. Images git already tracked, such as the 2026-06-04 heatmaps, are unaffected.

- `raw_data/092424_neg_control_plate{1,2}_{endpoint,heatmaps,maxOD}.png` (6 files, 2026-08-26)
- `raw_data/2026_08_21_VR_PosControl_{endpoint_t8.5h,heatmaps,maxOD}.png` (3 files, 2026-08-26)
- `raw_data/2026-06-04/IMG_3089.jpeg` (2026-09-04)

**Old local copies:** `Glucose_Varied/` and `PlateErrorGP/` are ignored on purpose (reference copies of the 2021 project and package).

`PlateErrorGP/` is its own git repo, a clone of `jensenlab/PlateErrorGP` (27 of its 32 commits are Vasu's). It has uncommitted changes that predate its copy to goblin on 2026-10-06: an edited `GP_Notebook.ipynb`, one extra `.gitignore` line, and two untracked data files that are identical to the copies in `raw_data/poster_e4/`. **Don't discard the notebook edits.** The poster reproduction in `experiments/run_poster_e4.jl` cites "cell 56" and "cell 62" of this edited version; in the committed version those cells are 55 and 61. Both versions contain the printed values. VS Code is set to hide this repo from Source Control (`git.ignoredRepositories` in `~/.vscode-server/data/Machine/settings.json`).

## Run records

Each results folder from 2026-10-09 has a `RUN_INFO.txt` and a `RUN_DIFF_*.patch`, written by `write_run_info`. Together they record the commit and the uncommitted code that produced it, so a result can be tied to code even if it ran before a commit:

- `results/`: sweep, started 2026-10-09 01:22 on `6755b12` plus `RUN_DIFF_20261009T012242.patch`
- `results/poster_e4/`: started 01:47 on `6755b12` plus `RUN_DIFF_20261009T014730.patch`
- `results/full_n/`: no RUN_INFO (that script doesn't call `write_run_info`). The code it ran is the code in `62042b4`. Run log: `logs/full_n_2026-10-09.log`.
