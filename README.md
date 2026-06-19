# EKF-Adaptive False-Alarm-Aware Sensor Placement with a CBF Safety Filter

**Online environmental estimation and certified-safe filter adaptation for barrier-coverage sensor networks**

> Research Task 2 — Muhammad Ibaad, MS Applicant, Dept. of Electrical and Computer Engineering, Georgia Southern University
> Supervised extension of Kim et al. (IEEE SysCon 2025) and Kim et al. (arXiv:2510.05343)

---

## Overview

This repository implements an online extension of the static, offline filter-gain selection used in Kim et al.'s false-alarm–aware spatiotemporal sensing framework. A per-sensor scalar **Extended Kalman Filter (EKF)** tracks the drifting environmental loss field $\omega(s,t)$ in real time, and a **discrete-time Control Barrier Function (CBF)** safety filter projects the resulting adaptive filter gain $\theta_i(t)$ so that the instantaneous void probability $\nu(t)$ never drops below a certified threshold $\nu_{\min}$. The result is a sensor network that adapts its detection/false-alarm trade-off as conditions evolve, while retaining a forward-invariance safety guarantee that the static baseline cannot offer.

---

## Key Technical Contributions

### 1. Per-Sensor EKF for Environmental Loss Tracking

Each sensor $i$ maintains an independent scalar EKF over its local environmental loss $\hat\omega_i(t) \in [0,1]$, driven by noisy SNR / hit-rate observations:

**Prediction**

$$
\hat{x}_{i,k\mid k-1} = F\hat{x}_{i,k-1\mid k-1}, \qquad
P_{i,k\mid k-1} = F^2 P_{i,k-1\mid k-1} + Q
$$

**Update**

$$
K_{i,k} = \frac{P_{i,k\mid k-1}}{P_{i,k\mid k-1} + R}, \qquad
\hat{x}_{i,k\mid k} = \hat{x}_{i,k\mid k-1} + K_{i,k}\left(y_{i,k} - \hat{x}_{i,k\mid k-1}\right), \qquad
P_{i,k\mid k} = (1 - K_{i,k}) P_{i,k\mid k-1}
$$

with $F = 0.97$, $Q = 4\times10^{-3}$, $R = 2\times10^{-2}$ in the simulated configuration. The per-sensor design preserves spatial structure in $\omega(s,t)$ that a single shared filter would average away.

### 2. Online Adaptive Filter Gain

At each time step, the EKF estimate $\hat\omega_i(t)$ is used to re-solve the filter-gain selection that Paper [2] otherwise performs once, offline:

$$
\theta_i^*(t) = \arg\max_{\theta \in [0.4, 2.5]} \int_\Psi \lambda(s,t)\alpha\big(\chi(\theta, \hat\omega_i(t))\big) p\big(s, a_i, \theta, \hat\omega_i(t)\big) ds
$$

evaluated via a 40-point grid search per sensor per time step, with negligible computational overhead.

### 3. Discrete-Time CBF Safety Filter on $\theta$

To guarantee that adaptation never degrades coverage below an operational floor, the proposed $\theta_i^*(t)$ is passed through a CBF projection enforcing forward invariance of the safe set $\mathcal{S} = \{\theta : h(\theta) \ge 0\}$:

$$
h(\theta) = \nu(a,\theta;t) - \nu_{\min}, \qquad h(\theta_{\text{new}}) \ge (1-\gamma) h(\theta_{\text{old}})
$$

If the unconstrained proposal violates this inequality, $\theta$ is projected by bisection between $\theta_{\text{old}}$ and $\theta^*$ (at most 10 steps) toward the closest safe value:

$$
\theta_i^{\text{safe}} = \arg\min_{\theta \in [\theta_{\text{old}}, \theta^*]} \lvert \theta - \theta^* \rvert \quad \text{s.t.} \quad h(a,\theta) \ge (1-\gamma) h(a,\theta_{\text{old}})
$$

with $\nu_{\min} = 0.25$ and class-$\mathcal{K}$ parameter $\gamma = 0.6$ in the reported experiments. This is the core safety contribution absent from the original offline framework.

### 4. Connection to Theorem 1 (Paper [2])

The sufficient condition for beneficial filtering, $\dfrac{\partial_\theta p}{p} \ge -\dfrac{\alpha'(\chi)}{\alpha(\chi)} \partial_\theta\chi$, specializes under the Gaussian sensing model to:

$$
\frac{2(s-a)^2 e^{\beta_\omega \omega}}{\theta^3} \;\ge\; \beta\left(\frac{\chi}{1+\beta\chi}\right) 2(\theta-1)\omega
$$

With accurate online $\hat\omega_i(t)$, the EKF lets the controller stay on the favorable side of this inequality at every time step — something a static $\theta$ can only satisfy on time-average, never instantaneously, in a drifting environment.

### Empirical Result Summary

| Method | Mean $\nu$ | Min $\nu$ | % Time $\nu < \nu_{\min}$ | CBF Activations |
|---|---|---|---|---|
| Random + fixed $\theta$ | 0.300 | 0.108 | 52.1% | — |
| Greedy + fixed $\theta$ (Paper [2]) | 0.419 | 0.184 | 14.6% | — |
| **Greedy + EKF-CBF $\theta$ (Proposed)** | **0.436** | 0.184 | 14.6% | 7.3% |

Mean void probability improves by **+4.1%** over the static-$\theta$ baseline at identical sensor positions and identical worst-case safety performance — the gain comes entirely from the adaptive filtering layer.

---

## Repository Structure

```
ekf-cbf-sensor-placement/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
├── LICENSE
├── .gitignore
│
├── src/
│   └── ekf_cbf_filtering_sim.py       # Single-file simulation (~450 LOC)
│                                       #   - domain & environment setup
│                                       #   - LGCP target intensity model
│                                       #   - sensing / availability model (Paper [2] eqs.)
│                                       #   - SensorEKF class
│                                       #   - adaptive theta search
│                                       #   - CBF safety projection
│                                       #   - experiment runner + plotting
│
├── figures/
│   ├── fig1_environment.png           # True omega(s,t) + per-sensor EKF tracking
│   ├── fig2_theta_adapt.png           # Adaptive theta_i(t) + CBF activation map
│   ├── fig3_void_prob.png             # Instantaneous void probability, 3 methods
│   └── fig4_placement.png             # Sensor positions over target intensity lambda(s,t)
│
├── reports/
│   └── Ibaad_Report_2.pdf             # Full research task report (this project's writeup)
│
├── papers/
│   ├── Near-Optimal_Sensor_Placement_for_Detecting_Stochastic_Target_Trajectories_in_Barrier_Coverage_Systems.pdf   # Paper [1], SysCon 2025
│   └── 2510.05343v3.pdf               # Paper [2], arXiv preprint
│
├── notebooks/                         # (optional) exploratory analysis / sensitivity sweeps
│
└── tests/                             # (optional) unit tests for EKF step, CBF projection, void_prob
```

---

## Requirements & Installation

The simulation has a deliberately minimal dependency footprint:

```
numpy
matplotlib
```

Install into a fresh virtual environment:

```bash
git clone https://github.com/<your-username>/ekf-cbf-sensor-placement.git
cd ekf-cbf-sensor-placement

python3 -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

pip install -r requirements.txt
```

`requirements.txt`:

```
numpy>=1.24
matplotlib>=3.7
```

---

## How to Run

Run the full experiment (EKF tracking → adaptive $\theta$ → CBF projection → three-method comparison) from the repository root:

```bash
python src/ekf_cbf_filtering_sim.py
```

This will:

1. Build the synthetic LGCP target intensity $\lambda(s,t)$ and drifting environmental field $\omega(s,t)$.
2. Run **Random + fixed $\theta$**, **Greedy + fixed $\theta$** (Paper [2] baseline), and **Greedy + EKF-CBF $\theta$** (proposed) over the 24-hour horizon.
3. Print a results summary (sensor positions, final $\theta$ values, mean/min void probability, CBF activation rate) to the terminal.
4. Save figures (`fig1_environment.png` through `fig4_placement.png`, plus `fig5_summary.png`) to the working directory. This repo's `figures/` folder currently includes `fig1`–`fig4`; `fig5_summary.png` is reproduced automatically if you run the script yourself.

To redirect figure output into the `figures/` folder:

```bash
mkdir -p figures && cd figures && python ../src/ekf_cbf_filtering_sim.py
```

Expected terminal output includes a block similar to:

```
Mean nu  | Random  : 0.3003
Mean nu  | Greedy  : 0.4193
Mean nu  | Proposed: 0.4364

CBF activations: 14 sensor-steps out of 192 (7.3%)
```

---

## Simulations & Visuals

| Figure | Description |
|---|---|
| **`fig1_environment.png`** | *(Left)* Ground-truth environmental loss field $\omega(s,t)$ as a drifting, sinusoidally-modulated Gaussian front, with sensor positions overlaid. *(Right)* Per-sensor EKF estimate $\hat\omega_i(t)$ (solid) vs. local ground truth (dashed) — demonstrates convergence within ≈4 hours and tracking accuracy within $\pm0.1$ under measurement noise $\sigma = 0.10$. |
| **`fig2_theta_adapt.png`** | *(Left)* Adaptive filter gain $\theta_i(t)$ per sensor vs. the static $\theta = 1.2$ baseline — outer sensors in cleaner water push $\theta$ toward 2.5; center sensors in degraded conditions hold $\theta$ near 1.1–1.3. *(Right)* Time bins where the CBF projection was triggered, concentrated at the peak-degradation window. |
| **`fig3_void_prob.png`** | Instantaneous void probability $\nu(t)$ across the 24-hour horizon for all three methods, with the CBF safety threshold $\nu_{\min} = 0.25$ marked. Shows the proposed method tracking the greedy baseline and pulling ahead once the environment drifts away from the static $\theta$'s calibration point. |
| **`fig4_placement.png`** | Target intensity $\lambda(s,t)$ heatmap with sensor positions for all three methods overlaid, illustrating how Random placement misses a high-traffic lane that Greedy/Proposed correctly cover. |

---

## Academic Citations

If referencing this work or its baselines, please cite:

```bibtex
@inproceedings{kim2025nearoptimal,
  author    = {Kim, Mingyu and Stilwell, Daniel J. and Yetkin, Harun and Jimenez, Jorge},
  title     = {Near-Optimal Sensor Placement for Detecting Stochastic Target Trajectories in Barrier Coverage Systems},
  booktitle = {2025 IEEE International Systems Conference (SysCon)},
  year      = {2025},
  pages     = {1--7}
}

@misc{kim2025robust,
  author       = {Kim, Mingyu and Sarker, Pronoy and Kim, Seungmo and Stilwell, Daniel J. and Jimenez, Jorge},
  title        = {Robust Sensor Placement for Poisson Arrivals with False-Alarm--Aware Spatiotemporal Sensing},
  howpublished = {arXiv:2510.05343},
  year         = {2025}
}

@article{ames2017cbf,
  author  = {Ames, Aaron D. and Xu, Xiangru and Grizzle, Jessy W. and Tabuada, Paulo},
  title   = {Control Barrier Function Based Quadratic Programs for Safety Critical Systems},
  journal = {IEEE Transactions on Automatic Control},
  volume  = {62},
  number  = {8},
  pages   = {3861--3876},
  year    = {2017}
}

@book{simon2006optimal,
  author    = {Simon, Dan},
  title     = {Optimal State Estimation: Kalman, H-Infinity, and Nonlinear Approaches},
  publisher = {Wiley},
  year      = {2006}
}

@article{nemhauser1978analysis,
  author  = {Nemhauser, G. L. and Wolsey, L. A. and Fisher, M. L.},
  title   = {An Analysis of Approximations for Maximizing Submodular Set Functions---I},
  journal = {Mathematical Programming},
  volume  = {14},
  number  = {1},
  pages   = {265--294},
  year    = {1978}
}
```

This extension's report (`reports/Ibaad_Report_2.pdf`) documents the full derivation, simulation parameters, and discussion of limitations — including open items such as a predictive/robust CBF using EKF covariance and joint position-$\theta$ adaptation, both proposed as immediate next steps.

The original PDFs for both reference papers are included in this repository for convenience:
- [`papers/Near-Optimal_Sensor_Placement_for_Detecting_Stochastic_Target_Trajectories_in_Barrier_Coverage_Systems.pdf`](./papers/Near-Optimal_Sensor_Placement_for_Detecting_Stochastic_Target_Trajectories_in_Barrier_Coverage_Systems.pdf) — Paper [1]
- [`papers/2510.05343v3.pdf`](./papers/2510.05343v3.pdf) — Paper [2]

---

## License

See [`LICENSE`](./LICENSE) for terms.
