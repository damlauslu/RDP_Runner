# RDP Runner

A Tkinter-based implementation of the Relationship Diagramming Process (RDP) with an interactive REL matrix editor, step-by-step execution, layout visualization, and report generation.

## Features
- **Responsive 3‑pane layout** (left / center / right) with user-resizable sashes.
- **REL matrix editor** with symmetric updates and color-coded relationship strength.
- **Step-by-step or run‑to‑end** execution of the RDP selection and layout phases.
- **Layout visualization** with scrollable canvas and candidate WPV “ghosts.”
- **TCR table + Sequence (Pi)** display.
- **Report preview (in-app)** and **export to `.txt`**.
- **Guide** and **RDP steps info** popups.
- **Two reset modes:**
  - Reset All (clears everything, matrix -> U)
  - Reset Steps (keeps matrix, clears progress)

## Requirements
- Python 3.9+ (tested with 3.11)
- Tkinter (usually bundled with Python)

## Run
From the project folder:

```bash
python RDP.py
```

## UI Overview
**Left pane**
- Number of Departments (2–20, auto-generates matrix)
- Randomize Relations
- Configure Values (A/E/I/O/U/X)
- Start Calculation
- Next Step / Jump to Layout / Run to End
- Reset All / Reset Steps
- Preview Report / Export Report (.txt)
- Guide / RDP Steps Info
- Log panel

**Center pane**
- REL Matrix (scrollable)
- Layout Visualization (scrollable)

**Right pane**
- TCR Table
- Sequence (Pi)
- Placement Log & WPV

## Core Behavior Notes
- **Calculations only begin after `Start Calculation`.**
- If the matrix is edited or randomized, outputs are cleared and need recalculation.
- `Preview Report` and `Export Report` reflect the current model state.
- Jump to Layout generates selection steps and prepares layout without placing departments.

## Reporting
Reports include:
- Inputs and REL matrix
- TCR (if computed)
- Final sequence (if computed)
- Selection steps (if computed)
- Layout steps (if computed)
- Final layout positions + ASCII grid

## Files
- `RDP.py` — main application
- `README.md` — overview
- `USAGE.md` — step-by-step usage guide

## Troubleshooting
- **Center column not widest on first load:** the app enforces sash positions after geometry is computed. If your OS/theme overrides, resize once to reapply.
- **No TCR/Sequence shown:** press `Start Calculation`.
- **Logs not updating:** use `Start Calculation`, then `Next Step` or `Run to End`.

## License
See `LICENSE`.
