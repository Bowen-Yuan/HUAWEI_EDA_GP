# adaptive_lambda_gp parameter units

- iterations: optimizer iterations for this stage.
- direction.epsilon_bin_scales: HPWL epsilon as multiples of the density bin size; each scale yields one exact active direction, RMS-normalized, then minimum-norm combined on the simplex.
- optimizer: forwarded to the shared ea::make_optimizer; learning_rate is the initial/constant rate.
- step_policy: constant, cosine, or trust. Bin-suffixed fields convert with the smaller density bin dimension. max_backtracks bounds the shared exact backtracking loop (scale = 1, 1/2, ...).
- lambda_policy: funnel_pid (log-domain PID on the overflow funnel) or fixed (ablation). Percent fields are overflow percentages (15.0 means 15%).
- acceptance: exact_funnel only. Percent budgets are relative HPWL costs for recovery steps.
- reset: optimizer state reset triggers; never applied to the lambda controller.
- verbose: print per-iteration key=value telemetry to stdout (console only; never persisted).
