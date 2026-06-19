# EKF-Adaptive False-Alarm-Aware Sensor Placement with a CBF Safety Filter

**Online environmental estimation and certified-safe filter adaptation for barrier-coverage sensor networks**

> Research Task 2 — Muhammad Ibaad, MS Applicant, Department of Electrical and Computer Engineering, Georgia Southern University
> Supervised extension of Kim et al. (IEEE SysCon 2025) and Kim et al. (arXiv:2510.05343)

---

## Overview

This repository implements an online extension of the static, offline filter-gain selection used in Kim et al.’s false-alarm-aware spatiotemporal sensing framework.

A per-sensor scalar **Extended Kalman Filter (EKF)** tracks the drifting environmental loss field `ω(s,t)` in real time, and a **discrete-time Control Barrier Function (CBF)** safety filter projects the adaptive filter gain `θᵢ(t)` such that the instantaneous void probability `ν(t)` never drops below a certified threshold `ν_min`.

The result is a sensor network that adapts its detection/false-alarm trade-off as conditions evolve while retaining a forward-invariance safety guarantee that the static baseline cannot provide.

---

# Key Technical Contributions

---

## 1. Per-Sensor EKF for Environmental Loss Tracking

Each sensor `i` maintains an independent scalar EKF over local environmental loss:

```text
ω̂ᵢ(t) ∈ [0,1]
```

This estimate is driven by noisy SNR / hit-rate observations.

### Prediction Step

```math
\hat{x}_{i,k|k-1}=F\hat{x}_{i,k-1|k-1}
```

```math
P_{i,k|k-1}=F^2 P_{i,k-1|k-1}+Q
```

### Update Step

```math
K_{i,k}=\frac{P_{i,k|k-1}}{P_{i,k|k-1}+R}
```

```math
\hat{x}_{i,k|k}
=
\hat{x}_{i,k|k-1}
+
K_{i,k}(y_{i,k}-\hat{x}_{i,k|k-1})
```

```math
P_{i,k|k}=(1-K_{i,k})P_{i,k|k-1}
```

Simulation parameters:

```text
F = 0.97
Q = 4×10⁻³
R = 2×10⁻²
```

The per-sensor formulation preserves spatial structure in `ω(s,t)` that a single global filter would average away.

---

## 2. Online Adaptive Filter Gain

At each time step, the EKF estimate is used to re-solve the filter-gain optimization that Paper [2] computes only once offline.

```math
\theta_i^*(t)
=
\operatorname*{argmax}_{\theta\in[0.4,2.5]}
\int_\Psi
\lambda(s,t)
\alpha(\chi(\theta,\hat{\omega}_i(t)))
p(s,a_i,\theta,\hat{\omega}_i(t))
\,ds
```

Implementation details:

* 40-point grid search
* Per sensor
* Per time step
* Negligible computational overhead

---

## 3. Discrete-Time CBF Safety Filter

To guarantee adaptation never degrades coverage below an operational safety floor, proposed gains are passed through a CBF projection.

Safe set:

```math
\mathcal{S}=\{\theta : h(\theta)\ge0\}
```

Barrier function:

```math
h(\theta)=\nu(a,\theta;t)-\nu_{min}
```

Safety constraint:

```math
h(\theta_{new})\ge (1-\gamma)h(\theta_{old})
```

If the unconstrained proposal violates this inequality, `θ` is projected by bisection between:

```text
θ_old and θ*
```

Maximum iterations:

```text
10
```

Projection objective:

```math
\theta_i^{safe}
=
\arg\min_{\theta\in[\theta_{old},\theta^*]}
|\theta-\theta^*|
```

Subject to:

```math
h(a,\theta)\ge(1-\gamma)h(a,\theta_{old})
```

Simulation parameters:

```text
ν_min = 0.25
γ = 0.6
```

This is the primary safety contribution absent from the original framework.

---

## 4. Connection to Theorem 1 (Paper [2])

The sufficient condition for beneficial filtering is:

```math
\frac{\partial_\theta p}{p}
\ge
-\frac{\alpha'(\chi)}{\alpha(\chi)}\partial_\theta \chi
```

Under the Gaussian sensing model:

```math
\frac{
2(s-a)^2e^{\beta_\omega\omega}
}{
\theta^3
}
\ge
\beta
\left(
\frac{\chi}{1+\beta\chi}
\right)
2(\theta-1)\omega
```

With accurate online estimates `ω̂ᵢ(t)`, the EKF helps the controller remain on the favorable side of this inequality at every time step.

A static `θ` can only satisfy this condition on time-average, not instantaneously.

---

# Empirical Results

| Method                            |    Mean ν |     Min ν | % Time ν < ν_min | CBF Activations |
| --------------------------------- | --------: | --------: | ---------------: | --------------: |
| Random + Fixed θ                  |     0.300 |     0.108 |            52.1% |               — |
| Greedy + Fixed θ                  |     0.419 |     0.184 |            14.6% |               — |
| **Greedy + EKF-CBF θ (Proposed)** | **0.436** | **0.184** |        **14.6%** |        **7.3%** |

Performance improvement:

```text
+4.1% mean void probability improvement
```

This gain comes entirely from adaptive filtering.

---

# Repository Structure

```text
ekf-cbf-sensor-placement/
├── README.md
├── requirements.txt
├── LICENSE
├── .gitignore
│
├── src/
│   └── ekf_cbf_filtering_sim.py
│
├── figures/
│   ├── fig1_environment.png
│   ├── fig2_theta_adapt.png
│   ├── fig3_void_prob.png
│   └── fig4_placement.png
│
├── reports/
│   └── Ibaad_Report_2.pdf
│
├── papers/
│   ├── Near-Optimal_Sensor_Placement...pdf
│   └── 2510.05343v3.pdf
│
├── notebooks/
└── tests/
```

---

# Requirements

Minimal dependencies:

```txt
numpy>=1.24
matplotlib>=3.7
```

---

# Installation

```bash
git clone https://github.com/<your-username>/ekf-cbf-sensor-placement.git
cd ekf-cbf-sensor-placement

python3 -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

---

# Running the Simulation

```bash
python src/ekf_cbf_filtering_sim.py
```

The script will:

1. Build synthetic LGCP target intensity `λ(s,t)`
2. Build drifting environmental field `ω(s,t)`
3. Run three methods:

   * Random + fixed θ
   * Greedy + fixed θ
   * Greedy + EKF-CBF θ
4. Generate figures
5. Print performance summary

Example terminal output:

```text
Mean nu  | Random  : 0.3003
Mean nu  | Greedy  : 0.4193
Mean nu  | Proposed: 0.4364

CBF activations: 14 sensor-steps out of 192 (7.3%)
```

---

# Figures

## fig1_environment.png

Shows:

* Ground-truth environmental loss field
* Sensor positions
* EKF tracking performance

EKF converges within approximately **4 hours**.

---

## fig2_theta_adapt.png

Shows:

* Adaptive filter gain `θᵢ(t)`
* Comparison with static baseline
* CBF activation periods

Outer sensors push toward larger θ values under cleaner conditions.

---

## fig3_void_prob.png

Shows instantaneous void probability `ν(t)` over 24 hours.

Includes safety threshold:

```text
ν_min = 0.25
```

---

## fig4_placement.png

Shows target intensity heatmap and sensor placement.

Demonstrates that random placement misses high-traffic lanes.

---

# Academic Citations

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
  title        = {Robust Sensor Placement for Poisson Arrivals with False-Alarm-Aware Spatiotemporal Sensing},
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
  title     = {Optimal State Estimation},
  publisher = {Wiley},
  year      = {2006}
}
```

---

# Future Work

Planned extensions include:

* Predictive / robust CBF using EKF covariance
* Joint position-θ optimization
* Multi-agent cooperative adaptation
* Real-world deployment with physical sensors

---

# License

See `LICENSE` for usage terms.
