#checks if the spearman and binomial tests are relevant using bootstrapping
import numpy as np
from scipy.stats import spearmanr, binomtest
import matplotlib.pyplot as plt


grad = np.load('grad.npy')
baz = np.load('baz.npy')
tol = .5

mask = (np.abs(baz) > tol)

grad_f = grad[mask]
baz_f = baz[mask]

#random resampling
rng = np.random.default_rng(seed=442)
n_perm = 10000

#observed statistics
obs_rho, _ = spearmanr(grad_f, baz_f)

sign_agree  = np.sign(grad_f) == np.sign(baz_f)
obs_k       = np.sum(sign_agree)
obs_n       = len(sign_agree)
obs_pct     = obs_k / obs_n          # 0.62

#permutation distributions
# Null: shuffle baz_f and compute both statistics

perm_rho = np.empty(n_perm)
perm_pct = np.empty(n_perm)

for i in range(n_perm):
    baz_shuf      = rng.permutation(baz_f)
    perm_rho[i]   = spearmanr(grad_f, baz_shuf).statistic
    perm_pct[i]   = np.mean(np.sign(grad_f) == np.sign(baz_shuf))

#permutation p-values (one-tailed: observed > null)
p_perm_rho = np.mean(perm_rho >= obs_rho)
p_perm_pct = np.mean(perm_pct >= obs_pct)

#analytic binomial test for sign agreement
binom_res = binomtest(obs_k, obs_n, p=0.5, alternative='greater')

print("------\n")
print(f"  N observations          : {obs_n}")
print(f"  N permutations          : {n_perm:,}")
print("------\n")
print("  Spearman rho")
print(f"    Observed              : {obs_rho:.4f}")
print(f"    Permutation mean      : {perm_rho.mean():.4f} ")
print(f"    Permutation p-value   : {p_perm_rho:.4f}")
print("------\n")
print("  Sign-agreement rate")
print(f"    Observed              : {obs_pct:.4f}  ({obs_k}/{obs_n})")
print(f"    Permutation mean      : {perm_pct.mean():.4f} ") #  (should be around 0.50)
print(f"    Permutation p-value   : {p_perm_pct:.4f}")
print(f"    Binomial p-value      : {binom_res.pvalue:.4f}")
print("------\n")


fig, axes = plt.subplots(1, 2, figsize=(12, 4))
plt.rcParams.update({'font.size': 13.5})

for ax, perm_vals, obs_val, xlabel, label in zip(
    axes,[perm_rho,perm_pct],[obs_rho,obs_pct],
    ["Spearman ρ", "Sign-agreement r"],
    [f"obs ρ={obs_rho:.2f}", f"obs r={obs_pct:.2f}"],):

    ax.hist(perm_vals, bins=50, color="orchid", alpha=0.85, edgecolor="white")
    ax.axvline(obs_val, color="navy", lw=1.5,ls='--', label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.legend(loc='upper left')#,frameon=False
    ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(f" N={n_perm:,} shuffles, n={obs_n} obs")
fig.tight_layout()
fig.savefig("permutation_null_dist.png", dpi=350,bbox_inches='tight', pad_inches=0.1)
plt.show()
