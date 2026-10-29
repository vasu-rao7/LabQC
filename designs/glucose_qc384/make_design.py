#!/usr/bin/env python3
"""Generate the 384-well glucose QC plate layout described in request.md.

The plate is run as a 96-well experiment: each 96-well position is a 2x2 block of
384 wells (96-well A1 = 384-well A1, A2, B1, B2). Every block has a positive control
(top left), a negative control (bottom left) and two replicate wells of the block's
condition (top right, bottom right).

Blocks: 12 positive control blocks, 12 negative control blocks and 72 test blocks
(24 glucose levels x 3 blocks).

Control blocks: Latin hypercube on the 8x12 block grid, scored with the same `LHS`
objective as PlateArrays.jl (src/solvers/objectives.jl), with no two control blocks
sharing an edge. Ties are broken with its `minimax` objective (total distance from
every block to its nearest positive and nearest negative control), so the controls are
also spread out.

Test blocks: random, except that the 3 blocks of a level sit in different block rows
and columns, never touch (including diagonally) and are not all on the outer ring.
Among layouts that pass, keep the one whose glucose levels are most even across block
rows, block columns and rings (distance from the plate edge).

Usage: python3 make_design.py   (writes the CSVs next to this script)
"""

import csv
import random
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np

SEED = 20261012
OUT_DIR = Path(__file__).resolve().parent
PREFIX = "glucose_qc384"

MAX_GLUCOSE = 10.0  # g/L; standard CDM (5 g glucose per 500 mL)
N_LEVELS = 24       # test levels, evenly spaced strictly between 0 and MAX_GLUCOSE
N_REPS = 3          # test blocks per level
N_POS = N_NEG = 12  # control blocks of each type
R, C = 8, 12        # block grid (96-well plate)

NEG_LEVEL, POS_LEVEL = 0, N_LEVELS + 1
LETTERS = "ABCDEFGHIJKLMNOP"

# (block row offset, block column offset) -> position name and what the well holds
WELLS_IN_BLOCK = {
    (0, 0): ("TL", "positive_control"),
    (0, 1): ("TR", None),  # block condition
    (1, 0): ("BL", "negative_control"),
    (1, 1): ("BR", None),  # block condition (replicate of TR)
}

BLOCK_ROW = np.repeat(np.arange(R), C)
BLOCK_COL = np.tile(np.arange(C), R)
COL_DIST = np.abs(BLOCK_COL[:, None] - np.arange(C)[None, :])


def glucose(level):
    """g/L for a level index (0 = none, POS_LEVEL = max)."""
    return round(MAX_GLUCOSE * level / POS_LEVEL, 4)


def ring(cell):
    r, c = cell
    return min(r, c, R - 1 - r, C - 1 - c) + 1


def lhs(pos, neg):
    """PlateArrays.jl `LHS`: squared deviation of control counts per row and per column
    from an even spread. 4.0 is the minimum for 12 + 12 controls on an 8x12 grid."""
    score = 0.0
    for cells in (pos, neg):
        n = len(cells)
        rows = np.bincount([r for r, _ in cells], minlength=R)
        cols = np.bincount([c for _, c in cells], minlength=C)
        score += ((rows - n / R) ** 2).sum() + ((cols - n / C) ** 2).sum()
    return score


def minimax(pos, neg):
    """PlateArrays.jl `minimax`: summed Manhattan distance from every block to its nearest
    positive and nearest negative control. pos[c] / neg[c] is the control row in column c."""
    dp = np.abs(BLOCK_ROW[:, None] - np.asarray(pos)[None, :]) + COL_DIST
    dn = np.abs(BLOCK_ROW[:, None] - np.asarray(neg)[None, :]) + COL_DIST
    return int(dp.min(axis=1).sum() + dn.min(axis=1).sum())


def touching(pos, neg):
    """Pairs of control blocks that share an edge."""
    n = sum(abs(p - q) == 1 for p, q in zip(pos, neg))
    return n + sum(len({pos[c], neg[c]} & {pos[c + 1], neg[c + 1]}) for c in range(C - 1))


def place_controls(rng, restarts=300, iterations=3000):
    """One positive and one negative control block per column; per row, 1-2 of each type
    and 3 in total; no two control blocks share an edge. Every such layout hits the LHS
    minimum; keep the lowest minimax."""
    assert N_POS == N_NEG == C and R < N_POS <= 2 * R
    best = None
    for _ in range(restarts):
        doubles = set(rng.sample(range(R), N_POS - R))  # rows with 2 positives and 1 negative
        pos = [r for r in range(R) for _ in range(2 if r in doubles else 1)]
        neg = [r for r in range(R) for _ in range(1 if r in doubles else 2)]
        rng.shuffle(pos)
        rng.shuffle(neg)

        def cost():
            # Each row holds 3 controls in total, so keeping positives at 1-2 per row
            # keeps negatives at 1-2 per row as well.
            bad_rows = sum(Counter(pos)[r] not in (1, 2) for r in range(R))
            clashes = sum(p == n for p, n in zip(pos, neg))
            return 1000 * (bad_rows + clashes + touching(pos, neg)) + minimax(pos, neg)

        def swap(move, a, b):  # every move is its own inverse
            if move < 0.4:
                pos[a], pos[b] = pos[b], pos[a]
            elif move < 0.8:
                neg[a], neg[b] = neg[b], neg[a]
            else:  # swap types within a column; row totals stay at 3
                pos[a], neg[a] = neg[a], pos[a]

        current = cost()
        for _ in range(iterations):
            a, b = rng.sample(range(C), 2)
            move = rng.random()
            swap(move, a, b)
            new = cost()
            if new <= current:
                current = new
            else:
                swap(move, a, b)
        if current < 1000 and (best is None or current < best[0]):
            best = (current, pos[:], neg[:])
    _, pos, neg = best
    return [(r, c) for c, r in enumerate(pos)], [(r, c) for c, r in enumerate(neg)]


def assign_levels(rng, cells, restarts=20, iterations=10000):
    """Give each test level N_REPS blocks in different rows and columns, not touching and
    not all on the outer ring; among those, keep the most even spread of levels."""
    rows = np.array([r for r, _ in cells])
    cols = np.array([c for _, c in cells])
    rings = np.array([ring(x) for x in cells])
    groups = [[f == g for g in np.unique(f)] for f in (rows, cols, rings)]

    def violations(levels):
        by_level = {}
        for cell, level in zip(cells, levels):
            by_level.setdefault(level, []).append(cell)
        v = 0
        for pts in by_level.values():
            v += len(pts) - len({r for r, _ in pts})
            v += len(pts) - len({c for _, c in pts})
            v += sum(max(abs(r1 - r2), abs(c1 - c2)) < 2 for (r1, c1), (r2, c2) in combinations(pts, 2))
            v += all(ring(p) == 1 for p in pts)
        return v

    def imbalance(levels):
        x = np.asarray(levels, dtype=float)
        mu = x.mean()
        return sum(m.sum() * (x[m].mean() - mu) ** 2 for factor in groups for m in factor)

    def cost(levels):
        return 1000 * violations(levels) + imbalance(levels)

    best = None
    for _ in range(restarts):
        levels = [k for k in range(1, N_LEVELS + 1) for _ in range(N_REPS)]
        rng.shuffle(levels)
        current = cost(levels)
        for _ in range(iterations):
            a, b = rng.sample(range(len(cells)), 2)
            if levels[a] == levels[b]:
                continue
            levels[a], levels[b] = levels[b], levels[a]
            new = cost(levels)
            if new <= current:
                current = new
            else:
                levels[a], levels[b] = levels[b], levels[a]
        if current < 1000 and (best is None or current < best[0]):
            best = (current, levels[:])
    return dict(zip(cells, best[1]))


def main():
    rng = random.Random(SEED)
    pos_cells, neg_cells = place_controls(rng)
    block_level = {cell: POS_LEVEL for cell in pos_cells}
    block_level.update({cell: NEG_LEVEL for cell in neg_cells})
    test_cells = [(r, c) for r in range(R) for c in range(C) if (r, c) not in block_level]
    block_level.update(assign_levels(rng, test_cells))

    block_type = {NEG_LEVEL: "negative_control", POS_LEVEL: "positive_control"}
    wells = []
    for (r, c), level in block_level.items():
        for (dr, dc), (position, fixed_type) in WELLS_IN_BLOCK.items():
            row, col = 2 * r + dr, 2 * c + dc
            if fixed_type == "positive_control":
                well_level = POS_LEVEL
            elif fixed_type == "negative_control":
                well_level = NEG_LEVEL
            else:
                well_level = level
            wells.append({
                "Well": f"{LETTERS[row]}{col + 1}",
                "Row": LETTERS[row],
                "Column": col + 1,
                "Block": f"{LETTERS[r]}{c + 1}",
                "Position_In_Block": position,
                "Block_Type": block_type.get(level, "test"),
                "Well_Type": block_type.get(well_level, "test"),
                "Glucose_g_per_L": glucose(well_level),
                "_row": row,
                "_col": col,
                "_level": well_level,
            })
    wells.sort(key=lambda w: (w["_row"], w["_col"]))

    columns = ["Well", "Row", "Column", "Block", "Position_In_Block", "Block_Type", "Well_Type", "Glucose_g_per_L"]
    with open(OUT_DIR / f"{PREFIX}_well_map.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(wells)

    grid = {(w["_row"], w["_col"]): w["Glucose_g_per_L"] for w in wells}
    with open(OUT_DIR / f"{PREFIX}_plate_grid.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Row"] + list(range(1, 2 * C + 1)))
        for row in range(2 * R):
            writer.writerow([LETTERS[row]] + [grid[row, col] for col in range(2 * C)])

    well_counts = Counter(w["_level"] for w in wells)
    where = {
        NEG_LEVEL: "BL well of every block + TR/BR wells of the 12 negative control blocks",
        POS_LEVEL: "TL well of every block + TR/BR wells of the 12 positive control blocks",
    }
    with open(OUT_DIR / f"{PREFIX}_glucose_levels.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Glucose_g_per_L", "Glucose_percent_w_v", "Well_Type", "Wells", "Where"])
        for level in range(POS_LEVEL + 1):
            writer.writerow([
                glucose(level),
                round(glucose(level) / 10, 4),
                block_type.get(level, "test"),
                well_counts[level],
                where.get(level, f"TR/BR wells of {N_REPS} test blocks"),
            ])

    # Summary for checking the design
    print(f"LHS score (PlateArrays, min 4.0): {lhs(pos_cells, neg_cells):.2f}")
    print(f"minimax score (PlateArrays): {minimax([r for r, _ in pos_cells], [r for r, _ in neg_cells])}")
    print("Block grid (P = positive control block, N = negative control block, else test level 1-24):")
    print("     " + " ".join(f"{c + 1:>3}" for c in range(C)))
    for r in range(R):
        cells = []
        for c in range(C):
            level = block_level[r, c]
            cells.append({POS_LEVEL: "P", NEG_LEVEL: "N"}.get(level, str(level)))
        print(f"  {LETTERS[r]}  " + " ".join(f"{x:>3}" for x in cells))
    print("Wells per glucose level:", dict(sorted((glucose(k), v) for k, v in well_counts.items())))


if __name__ == "__main__":
    main()
