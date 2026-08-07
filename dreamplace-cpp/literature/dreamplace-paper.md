# DREAMPlace Paper Notes

Source: Y. Lin et al., "DREAMPlace: Deep Learning Toolkit-Enabled GPU
Acceleration for Modern VLSI Placement," IEEE TCAD 40(4), 2021.
DOI: 10.1109/TCAD.2020.3003843.

Key reconstruction requirements:

- Objective: weighted-average wirelength plus electrostatic density energy.
- Initialization: layout center plus Gaussian noise with standard deviation
  0.1% of layout width and height.
- Poisson solve: DCT for charge, inverse cosine/sine transforms for potential
  and fields under Neumann boundary conditions.
- Density weight: initialize from the L1 wirelength/density gradient norm ratio;
  update with the RePlAce HPWL-driven rule (paper Equation 18).
- Solver: Nesterov accelerated gradient with Lipschitz/BB step estimation.
- Legalization: Tetris-like Greedy followed by Abacus on CPU.
- Paper Table II full-flow HPWL: adaptec1 73.22M, adaptec2 82.22M.

The public repository configuration additionally specifies seed 1000, 1000
iterations, target density 1.0, stop overflow 0.07, degree cutoff 100, fillers,
and a 2.5% first-step position-noise perturbation scaled by node dimensions.

