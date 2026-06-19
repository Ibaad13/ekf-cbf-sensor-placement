import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(7)

# =====================================================================
# 1. DOMAIN SETUP 
# =====================================================================
S_MIN, S_MAX = 0.0, 10.0      # km — 1D barrier corridor
T_MIN, T_MAX = 0.0, 24.0      # hours
N_S = 100                     # spatial bins  (100 m each)
N_T = 48                      # temporal bins (30 min each)
ds = (S_MAX - S_MIN) / N_S
dt = (T_MAX - T_MIN) / N_T
s_grid = np.linspace(S_MIN + ds/2, S_MAX - ds/2, N_S)
t_grid = np.linspace(T_MIN + dt/2, T_MAX - dt/2, N_T)

K = 4                          # number of sensors

# Paper [2] hyper-parameters
BETA   = 5.0                   # FA penalty in alpha(chi) = 1/(1+beta*chi)
XI     = 0.2                   # min FA penalty floor
BETA_W = 1.5                   # how strongly omega shrinks sensing range
THETA0 = 1.0                   # nominal filter gain (chi=0 here)
NU_MIN = 0.25                  # CBF safe threshold on void probability
GAMMA_CBF = 0.6                # CBF class-K parameter

# =====================================================================
# 2. TRUE ENVIRONMENT omega(s,t) — 
# =====================================================================
def true_omega(s, t):
    """Time-varying environmental sensing loss field in [0,1].
    Models a turbidity / weather front that drifts across the corridor.
    """
    drift_center = 2.0 + 6.0 * (t / T_MAX)        # front moves left->right
    spatial = np.exp(-((s - drift_center) ** 2) / (2 * 1.5 ** 2))
    daily   = 0.5 + 0.4 * np.sin(2 * np.pi * t / T_MAX - np.pi / 2)
    return np.clip(daily * spatial, 0.0, 1.0)

OMEGA_TRUE = np.array([[true_omega(s, t) for s in s_grid] for t in t_grid])  # (N_T,N_S)
# Prior omega used by the offline greedy placer (temporal mean of the true field).
# This represents what a designer knows at deployment time from historical data.
OMEGA_PRIOR = np.mean(OMEGA_TRUE, axis=0, keepdims=True) * np.ones_like(OMEGA_TRUE)

# =====================================================================
# 3. LGCP TARGET INTENSITY lambda(s,t) — peaks in two lanes
# =====================================================================
def true_lambda(s, t):
    lane1 = 1.4 * np.exp(-((s - 3.0) ** 2) / (2 * 0.7 ** 2))
    lane2 = 1.0 * np.exp(-((s - 7.5) ** 2) / (2 * 0.9 ** 2))
    rush  = 0.5 + 0.5 * np.exp(-((t - 9.0) ** 2) / 8.0) \
                + 0.4 * np.exp(-((t - 17.0) ** 2) / 8.0)
    return rush * (lane1 + lane2)

LAMBDA = np.array([[true_lambda(s, t) for s in s_grid] for t in t_grid])

# =====================================================================
# 4. SENSING MODEL  (eqs. 3-6 of paper [2])
# =====================================================================
def sensing_range(theta, omega):
    """ell(theta,omega) = theta * exp(-beta_omega * omega)."""
    return theta * np.exp(-BETA_W * omega) + 1e-6

def p_raw(s, a, theta, omega):
    """Raw detection probability of sensor at a for a target at s."""
    ell = sensing_range(theta, omega)
    return np.exp(-((s - a) ** 2) / ell)

def chi_fa(theta, omega):
    """False-alarm surrogate chi = omega * ((theta-1)^2 + xi)."""
    return omega * ((theta - 1.0) ** 2 + XI)

def alpha_avail(chi):
    """Availability alpha(chi) = 1 / (1 + beta*chi)."""
    return 1.0 / (1.0 + BETA * chi)

def p_eff(s, a, theta, omega):
    """Effective detection prob: tilde-p = alpha(chi) * p_raw."""
    return alpha_avail(chi_fa(theta, omega)) * p_raw(s, a, theta, omega)

# =====================================================================
# 5. VOID PROBABILITY (Jensen lower bound, eq. 7 of paper [2])
#    nu(a, theta) = exp( - sum_{s,t} lambda(s,t) * prod_i (1 - p_eff_i) )
# =====================================================================
def void_prob(a_list, theta_list, omega_field):
    """a_list: K sensor positions (km).  theta_list: K filter gains.
       omega_field: (N_T, N_S) environment.
       Returns the void probability *per unit time* (averaged over T).
    """
    prod_miss = np.ones_like(LAMBDA)
    for a_k, th_k in zip(a_list, theta_list):
        p_eff_k = np.zeros_like(LAMBDA)
        for ti in range(N_T):
            p_eff_k[ti, :] = p_eff(s_grid, a_k, th_k, omega_field[ti, :])
        prod_miss *= (1.0 - p_eff_k)
    # expected undetected per unit time at each time bin
    U_t = np.sum(LAMBDA * prod_miss, axis=1) * ds       # length N_T
    # average void probability across the horizon
    return float(np.mean(np.exp(-U_t)))


def void_prob_snapshot(a_list, theta_list, omega_snapshot, lam_snapshot):
    """Instantaneous void probability at a single time using snapshot fields.
    omega_snapshot, lam_snapshot are 1-D arrays of length N_S.
    """
    prod_miss = np.ones(N_S)
    for a_k, th_k in zip(a_list, theta_list):
        prod_miss *= (1.0 - p_eff(s_grid, a_k, th_k, omega_snapshot))
    U = np.sum(lam_snapshot * prod_miss) * ds
    return float(np.exp(-U))

# =====================================================================
# 6. EKF FOR omega(s,t) AT EACH SENSOR LOCATION
#    State per sensor:  x_k = omega_k  (scalar)
#    Dynamics:          x_{k+1} = F x_k + w_k,    w_k ~ N(0, Q)
#    Measurement:       y_k = x_k + v_k,          v_k ~ N(0, R)
#       (in practice y_k comes from the sensor's reported SNR or
#        empirical hit-rate; here we feed noisy ground-truth)
# =====================================================================
class SensorEKF:
    def __init__(self, x0=0.3, P0=0.1, F=0.97, Q=4e-3, R=2e-2):
        self.x = x0; self.P = P0
        self.F = F; self.Q = Q; self.R = R
        self.hist_x = [x0]; self.hist_P = [P0]

    def step(self, y):
        # Predict
        x_pred = self.F * self.x
        P_pred = self.F * self.P * self.F + self.Q
        # Update
        K = P_pred / (P_pred + self.R)
        self.x = float(np.clip(x_pred + K * (y - x_pred), 0.0, 1.0))
        self.P = float((1 - K) * P_pred)
        self.hist_x.append(self.x); self.hist_P.append(self.P)
        return self.x, self.P

# =====================================================================
# 7. ADAPTIVE THETA  (closed-form 1-D search over theta for each sensor)
#    Pointwise objective: maximise tilde-p = alpha(chi(theta,omega_hat))
#                                            * exp(-d^2 / ell(theta,omega_hat))
#    Done numerically over a small theta grid; reproduces the trade-off
#    discussed in Theorem 1 of paper [2].
# =====================================================================
THETA_GRID = np.linspace(0.4, 2.5, 40)

def optimal_theta(omega_hat, a_k, lambda_t_slice):
    """Pick theta that maximises integrated effective detection at sensor a_k
    given current omega estimate and current-time intensity slice."""
    best_th, best_val = THETA0, -1.0
    for th in THETA_GRID:
        # contribution to expected detections at this snapshot:
        det = p_eff(s_grid, a_k, th, omega_hat) * lambda_t_slice
        val = np.sum(det) * ds
        if val > best_val:
            best_val = val; best_th = th
    return best_th

# =====================================================================
# 8. CBF SAFETY FILTER on the theta update
#    Safety function h(theta) = nu(a, theta) - NU_MIN >= 0
#    CBF condition (discrete time, class-K = gamma * h):
#         h(theta_new) >= (1 - gamma) * h(theta_old)
#    If the proposed theta_new violates the inequality, we project by
#    line-search between theta_old and theta_new toward theta_old until
#    the condition is met.
# =====================================================================
def cbf_project(a_list, theta_old, theta_new, omega_snapshot, lam_snapshot):
    """Project theta_new toward theta_old until the discrete-time CBF
    inequality  h(x+) >= (1-gamma) h(x)  holds, where
        h(theta) = nu_instant(a, theta) - NU_MIN.
    Returns the safe theta_new vector and a 0/1 flag per sensor."""
    h_old = void_prob_snapshot(a_list, theta_old, omega_snapshot, lam_snapshot) - NU_MIN
    activated = np.zeros(len(theta_old))
    safe_new = list(theta_new)
    for _ in range(10):
        h_new = void_prob_snapshot(a_list, safe_new, omega_snapshot, lam_snapshot) - NU_MIN
        if h_new >= (1 - GAMMA_CBF) * h_old:
            break
        for i in range(len(safe_new)):
            if abs(safe_new[i] - theta_old[i]) > 1e-3:
                safe_new[i] = 0.5 * (safe_new[i] + theta_old[i])
                activated[i] = 1
    return safe_new, activated

# =====================================================================
# 9. SENSOR PLACEMENT METHODS
# =====================================================================
def random_placement(K=K, rng=None, min_sep=0.5):
    rng = rng or np.random.default_rng(3)
    out = []
    while len(out) < K:
        c = rng.uniform(S_MIN + 0.5, S_MAX - 0.5)
        if all(abs(c - p) > min_sep for p in out):
            out.append(c)
    return sorted(out)

def greedy_placement_fixed_theta(theta_fixed=1.2, min_sep=0.5):
    """Paper [2] baseline: greedy maximisation of (time-averaged) void
    probability under *static* filter gain theta_fixed.  A small minimum
    spacing avoids degenerate duplicates.

    Note: The greedy algorithm uses OMEGA_PRIOR (temporal mean of the
    true field) as the environmental prior -- this is what paper [2]
    would have at placement time.  The online EKF then corrects for the
    actual drifting realization during deployment.
    """
    candidates = np.linspace(S_MIN + 0.3, S_MAX - 0.3, 30).tolist()
    placed = []; thetas = []
    for _ in range(K):
        best_v, best_a = -1, None
        for c in candidates:
            if any(abs(c - p) < min_sep for p in placed):
                continue
            # Use prior omega (temporal mean) for offline greedy placement
            v = void_prob(placed + [c], thetas + [theta_fixed], OMEGA_PRIOR)
            if v > best_v:
                best_v, best_a = v, c
        if best_a is None:   # fallback if no candidate left
            best_a = candidates[0]
        placed.append(best_a); thetas.append(theta_fixed)
    return placed, thetas

# =====================================================================
# 10. MAIN ONLINE LOOP — runs all three methods through time
# =====================================================================
def run_experiment():
    # ---------- offline placement (same for all online methods) -------
    a_random  = random_placement()
    theta_random = [1.2] * K               # fixed theta for the random baseline
    a_greedy, theta_greedy = greedy_placement_fixed_theta(1.2)
    a_prop = list(a_greedy)                # proposed uses same positions

    # ---------- EKF bank for the proposed method ----------------------
    # initial guess: 0.3 (mild loss), large initial variance
    ekfs = [SensorEKF(x0=0.3, P0=0.2) for _ in range(K)]
    theta_prop_hist = []     # K x N_T
    cbf_hist        = []     # K x N_T (binary)
    omega_hat_hist  = []     # K x N_T
    nu_random_hist  = []
    nu_greedy_hist  = []
    nu_prop_hist    = []
    theta_prop = list(theta_greedy)

    for ti in range(N_T):
        omega_t = OMEGA_TRUE[ti, :]                # true env at this time
        lam_t   = LAMBDA[ti, :]

        # ---- EKF updates: each sensor samples noisy local omega ----
        omega_hats = []
        for i, ekf in enumerate(ekfs):
            # find nearest spatial bin to sensor i
            idx = int(np.clip(np.round((a_prop[i] - S_MIN) / ds - 0.5), 0, N_S - 1))
            y_meas = omega_t[idx] + np.random.randn() * 0.10  # sensor noise
            x_hat, _ = ekf.step(y_meas)
            omega_hats.append(x_hat)
        omega_hat_hist.append(omega_hats)

        # ---- adaptive theta from EKF estimates ----
        theta_proposed = [optimal_theta(omega_hats[i], a_prop[i], lam_t)
                          for i in range(K)]

        # ---- CBF safety projection ----
        theta_safe, activated = cbf_project(
            a_prop, theta_prop, theta_proposed, omega_t, lam_t)
        theta_prop = theta_safe
        theta_prop_hist.append(list(theta_prop))
        cbf_hist.append(activated.tolist())

        # ---- evaluate instantaneous void probability for each method ----
        nu_random_hist.append(void_prob_snapshot(a_random, theta_random, omega_t, lam_t))
        nu_greedy_hist.append(void_prob_snapshot(a_greedy, theta_greedy, omega_t, lam_t))
        nu_prop_hist.append(  void_prob_snapshot(a_prop,   theta_prop,   omega_t, lam_t))

    return {
        'a_random': a_random, 'theta_random': theta_random,
        'a_greedy': a_greedy, 'theta_greedy': theta_greedy,
        'a_prop':   a_prop,
        'theta_prop_hist': np.array(theta_prop_hist),
        'cbf_hist':        np.array(cbf_hist),
        'omega_hat_hist':  np.array(omega_hat_hist),
        'nu_random_hist':  np.array(nu_random_hist),
        'nu_greedy_hist':  np.array(nu_greedy_hist),
        'nu_prop_hist':    np.array(nu_prop_hist),
        'ekfs': ekfs,
    }

# =====================================================================
# 11. PLOTTING
# =====================================================================
def make_plots(R):
    # ---- Fig 1 : true environment + EKF-tracked omega per sensor ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    im = axes[0].imshow(OMEGA_TRUE, origin='lower', aspect='auto',
                        extent=[S_MIN, S_MAX, T_MIN, T_MAX], cmap='viridis')
    axes[0].set_title("True environment $\\omega(s,t)$")
    axes[0].set_xlabel("s [km]"); axes[0].set_ylabel("t [h]")
    for ak in R['a_prop']:
        axes[0].axvline(ak, color='red', lw=1, alpha=0.7)
    plt.colorbar(im, ax=axes[0], shrink=0.85)

    colors = ['C0', 'C1', 'C2', 'C3']
    for i in range(K):
        # true omega seen at sensor location through time
        idx = int(np.clip(np.round((R['a_prop'][i] - S_MIN) / ds - 0.5),
                          0, N_S - 1))
        true_series = OMEGA_TRUE[:, idx]
        axes[1].plot(t_grid, true_series, '--', color=colors[i], alpha=0.6,
                     label=f"true (s={R['a_prop'][i]:.1f})")
        axes[1].plot(t_grid, R['omega_hat_hist'][:, i], '-',
                     color=colors[i], lw=1.8,
                     label=f"EKF (s={R['a_prop'][i]:.1f})")
    axes[1].set_xlabel("t [h]"); axes[1].set_ylabel("$\\omega$")
    axes[1].set_title("EKF tracking of $\\omega$ at each sensor")
    axes[1].legend(fontsize=7, ncol=2); axes[1].grid(alpha=0.3)
    plt.tight_layout(); plt.savefig("fig1_environment.png", dpi=140); plt.close()

    # ---- Fig 2 : adaptive theta + CBF activation map ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for i in range(K):
        axes[0].plot(t_grid, R['theta_prop_hist'][:, i],
                     color=colors[i], label=f"sensor {i+1}", lw=1.6)
    axes[0].axhline(1.2, ls='--', color='k', alpha=0.6, label="fixed $\\theta$=1.2")
    axes[0].set_xlabel("t [h]"); axes[0].set_ylabel("$\\theta$")
    axes[0].set_title("Adaptive filter gain $\\theta_i(t)$ from EKF")
    axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

    axes[1].imshow(R['cbf_hist'].T, aspect='auto', cmap='Reds',
                   extent=[T_MIN, T_MAX, 0.5, K + 0.5])
    axes[1].set_xlabel("t [h]"); axes[1].set_ylabel("sensor index")
    axes[1].set_title("CBF activations (red = projection used)")
    axes[1].set_yticks(range(1, K + 1))
    plt.tight_layout(); plt.savefig("fig2_theta_adapt.png", dpi=140); plt.close()

    # ---- Fig 3 : void probability over time ----
    plt.figure(figsize=(8, 4))
    plt.plot(t_grid, R['nu_random_hist'], label="Random + fixed $\\theta$",
             color='#d62728', lw=1.5)
    plt.plot(t_grid, R['nu_greedy_hist'], label="Greedy + fixed $\\theta$  (paper [2])",
             color='#1f77b4', lw=1.8)
    plt.plot(t_grid, R['nu_prop_hist'],   label="Greedy + EKF-CBF $\\theta$ (proposed)",
             color='#2ca02c', lw=2.0)
    plt.axhline(NU_MIN, ls=':', color='k', label=f"CBF threshold $\\nu_{{min}}={NU_MIN}$")
    plt.xlabel("t [h]"); plt.ylabel("instantaneous void probability $\\nu$")
    plt.title("Void probability over time (higher is better)")
    plt.legend(fontsize=8); plt.grid(alpha=0.3)
    plt.tight_layout(); plt.savefig("fig3_void_prob.png", dpi=140); plt.close()

    # ---- Fig 4 : sensor placement on intensity map ----
    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(LAMBDA, origin='lower', aspect='auto',
                   extent=[S_MIN, S_MAX, T_MIN, T_MAX], cmap='YlOrRd')
    plt.colorbar(im, ax=ax, label="$\\lambda(s,t)$")
    for a in R['a_random']:  ax.axvline(a, color='red',    ls='--', lw=1.2)
    for a in R['a_greedy']:  ax.axvline(a, color='blue',   ls='-.', lw=1.2)
    for a in R['a_prop']:    ax.axvline(a, color='green',  ls='-',  lw=1.6)
    ax.set_xlabel("s [km]"); ax.set_ylabel("t [h]")
    ax.set_title("Target intensity & sensor positions  "
                 "(red=Random, blue=Greedy-fixed, green=Proposed)")
    plt.tight_layout(); plt.savefig("fig4_placement.png", dpi=140); plt.close()

    # ---- Fig 5 : summary bars ----
    methods = ['Random+fixed', 'Greedy+fixed [2]', 'EKF-CBF (prop.)']
    mean_nu = [R['nu_random_hist'].mean(),
               R['nu_greedy_hist'].mean(),
               R['nu_prop_hist'].mean()]
    min_nu  = [R['nu_random_hist'].min(),
               R['nu_greedy_hist'].min(),
               R['nu_prop_hist'].min()]
    viol    = [np.mean(R['nu_random_hist'] < NU_MIN),
               np.mean(R['nu_greedy_hist'] < NU_MIN),
               np.mean(R['nu_prop_hist']   < NU_MIN)]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    colors_b = ['#d62728', '#1f77b4', '#2ca02c']
    for ax, vals, ttl, ylab in [
        (axes[0], mean_nu, "Mean void probability", r"$\bar\nu$"),
        (axes[1], min_nu,  "Worst-case void probability", r"$\min \nu$"),
        (axes[2], viol,    "Fraction of time below threshold", "violation rate")]:
        ax.bar(methods, vals, color=colors_b, edgecolor='black')
        ax.set_title(ttl); ax.set_ylabel(ylab); ax.grid(alpha=0.3, axis='y')
        ax.set_ylim(0, max(vals) * 1.25 + 0.05)
        for i, v in enumerate(vals):
            ax.text(i, v + 0.01, f"{v:.3f}", ha='center', fontsize=9)
        plt.setp(ax.get_xticklabels(), rotation=10)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.18, wspace=0.3)
    plt.savefig("fig5_summary.png", dpi=140); plt.close()


# =====================================================================
# 12. RUN
# =====================================================================
if __name__ == "__main__":
    print("Running EKF-CBF false-alarm-aware sensor placement simulation ...")
    R = run_experiment()
    make_plots(R)

    print("\n================ RESULTS ================")
    print(f"Sensors (Random)   : {[f'{a:.2f}' for a in R['a_random']]}")
    print(f"Sensors (Greedy)   : {[f'{a:.2f}' for a in R['a_greedy']]}")
    print(f"Sensors (Proposed) : {[f'{a:.2f}' for a in R['a_prop']]}")
    print(f"Final theta (Prop) : {[f'{t:.2f}' for t in R['theta_prop_hist'][-1]]}")
    print()
    print(f"Mean nu  | Random  : {R['nu_random_hist'].mean():.4f}")
    print(f"Mean nu  | Greedy  : {R['nu_greedy_hist'].mean():.4f}")
    print(f"Mean nu  | Proposed: {R['nu_prop_hist'].mean():.4f}")
    print()
    print(f"Min  nu  | Random  : {R['nu_random_hist'].min():.4f}")
    print(f"Min  nu  | Greedy  : {R['nu_greedy_hist'].min():.4f}")
    print(f"Min  nu  | Proposed: {R['nu_prop_hist'].min():.4f}")
    print()
    print(f"Violations < nu_min={NU_MIN}:")
    print(f"  Random   : {np.mean(R['nu_random_hist'] < NU_MIN)*100:.1f}% of time")
    print(f"  Greedy   : {np.mean(R['nu_greedy_hist'] < NU_MIN)*100:.1f}% of time")
    print(f"  Proposed : {np.mean(R['nu_prop_hist']   < NU_MIN)*100:.1f}% of time")
    print(f"\nCBF activations: {int(R['cbf_hist'].sum())} sensor-steps "
          f"out of {N_T*K} ({R['cbf_hist'].mean()*100:.1f}%)")
    print("\nFigures saved: fig1..fig5_*.png")
