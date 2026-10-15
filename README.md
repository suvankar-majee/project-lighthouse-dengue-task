```
tasks/biology/mathematical-epidemiology/seasonal-dengue-transmission-potential/
├── instruction.md
├── task.toml
├── environment/   Dockerfile (all versions pinned) and data/
├── solution/      solve.sh + seasonal_r0.py (reference solution)
└── tests/         test.sh, test_outputs.py, reference.py (independent recomputation)
```

Root files: `solution_explanation.md`, `VALIDATION.md`, `DECLARATION.md`.

Domain: Biology, specifically mathematical epidemiology. The task asks for thermal-trait-driven dengue
transmission potential under a seasonal climate, which requires multi-laboratory data harmonisation,
climate-record homogenisation, survey-based calibration of a stage-structured vector model, and the
periodic basic reproduction number (Bacaër–Guernaoui / Wang–Zhao).

To validate:

```bash
harbor run -p tasks/biology/mathematical-epidemiology/seasonal-dengue-transmission-potential -a oracle   # mean 1.000
harbor run -p tasks/biology/mathematical-epidemiology/seasonal-dengue-transmission-potential -a nop      # mean 0.000
```
