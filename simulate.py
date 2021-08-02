from utils import *
import numpy as np
import pandas as pd
import pylab as plt


def gen_bursts(
    mu_params,
    sigma_params,
    mu_dist="uniform",
    sigma_dist="uniform",
    N=10000,
    alpha=-1.5,
    E_min_log=37,
    E_max_log=39,
    name="",
    save=True,
):
    """
    Generates a sample of FRBs and returns a dataframe.
    Each row corresponds to an FRB. Spectra is modeled
    using Gaussian, so each row consists of mu_f, sigma_f
    and energy of that FRB.

    mu_params: list containing the parameter of mu_f distribution
    sigma_params: list containing the parameter of sigma_f distribution
    mu_dist: distribution of mu_f
    sigma_dist: distribution of sigma_f
    N: Number of FRBs to generate
    alpha: Energy distribution slope
    E_min_log: Minimum energy cutoff in units of log10(E ergs)
    E_max_log: Maximum energy cutoff in units of log10(E ergs)
    name: Name of the output file
    save: Save the dataframe to a csv
    """
    if mu_dist == "uniform":
        low = mu_params[0]
        high = mu_params[1]
        mu_f = np.random.uniform(low=low, high=high, size=N)  # MHz
    elif mu_dist == "gauss" or mu_dist == "norm":
        loc = mu_params[0]
        scale = mu_params[1]
        mu_f = np.random.normal(loc=loc, scale=scale, size=N)  # MHz
    else:
        raise ValueError

    if sigma_dist == "uniform":
        sig_f = np.random.uniform(
            low=sigma_params[0], high=sigma_params[1], size=N
        )  # MHz
    elif sigma_dist == "gauss" or sigma_dist == "norm":
        loc = sigma_params[0]
        scale = sigma_params[1]
        sig_f = np.random.normal(loc=loc, scale=scale, size=N)  # MHz
    else:
        raise ValueError

    m = sig_f > 0
    N_new = m.sum()

    sig_f = sig_f[m]
    mu_f = mu_f[m]

    Es = []
    for i in range(N_new):
        Es.append(sample_from_powerlaw(E_min_log, E_max_log, alpha))
    Es = np.array(Es)

    vals = {}
    vals["in_mu_f"] = mu_f
    vals["in_sig_f"] = sig_f
    vals["E"] = Es
    df_vals = pd.DataFrame.from_dict(vals)

    mu_params_n = f"{mu_params[0]}_{mu_params[1]}"
    sigma_params_n = f"{sigma_params[0]}_{sigma_params[1]}"

    if len(name):
        name += "_"
    name += f"data_mu_{mu_dist}_{mu_params_n}_sig_{sigma_dist}_{sigma_params_n}_"
    name += f"alpha_{alpha}_Elog_{E_min_log}_{E_max_log}.csv"
    if save:
        df_vals.to_csv(name)
    return df_vals, name


def run_search(
    bursts, fstart, fend, fluence_threshold, in_band_sig=1, ret="all", distance=972
):
    """
    Runs a search on the bursts given in the input dataframe.
    Integrates the spectra within bandwidth bounds to determine
    the fraction of signal within the observing band. Checks
    if that value crosses the fluence threshold to determine
    if the burst is detected.
    Estimates the energy of the burst
    from snr-fluence and approx bandwidth of the burst. Here,
    the bandwidth is the burst signal seen within the observing band.
    It also estimates fit-energy using fit-fluence and actual
    bandwidth of the burst.
    Also, returns the energies of in-band bursts.

    bursts: dataframe with burst information
    fstart, fend (MHz): Observing band
    fluence_threshold: Jy ms
    in_band_sig: sigma to use to determine in-band bursts
    ret: To select the energies to return
    distance: (Mpc) Distance to the FRB

    """

    mu_f = bursts["in_mu_f"]
    sig_f = bursts["in_sig_f"]
    Es = bursts["E"]
    spectra_frac = gauss_integral(fstart, fend, mu_f, sig_f)

    S = energy_to_fluence(Es, sig_f, distance)  # Jy s
    assert len(np.where(sig_f < 0)[0]) == 0

    S_seen = S * spectra_frac  # fluence within the observable band
    bursts["in_S"] = S
    bursts["snr_S"] = S_seen
    fluence_threshold_jys = fluence_threshold / 1000
    detected_df = bursts[bursts["snr_S"] > fluence_threshold_jys]
    # ----------

    original_E = np.array(bursts["E"])

    u = np.array(detected_df["in_mu_f"] + 1 * detected_df["in_sig_f"])
    l = np.array(detected_df["in_mu_f"] - 1 * detected_df["in_sig_f"])

    uu = u >= fend
    ul = u <= fstart
    lu = l >= fend
    ll = l <= fstart
    u[uu & lu] = fend
    l[uu & lu] = fstart

    u[ll & ul] = fend
    l[ll & ul] = fstart

    u[uu & ll] = fend
    l[uu & ll] = fstart

    # u[ul & ll] = u[ul & ll]
    l[ul & ll] = fstart

    u[uu & lu] = fend
    # l[uu & lu] = l[uu & lu]

    freq_width = u - l
    assert len(np.where(freq_width < 0)[0]) == 0
    detected_snr_E = np.array(
        fluence_to_energy(detected_df["snr_S"], freq_width, distance)
    )
    detected_fit_E = np.array(detected_df["E"])

    upper = detected_df["in_mu_f"] + in_band_sig * detected_df["in_sig_f"] < fend
    lower = detected_df["in_mu_f"] - in_band_sig * detected_df["in_sig_f"] > fstart
    detected_in_band_df = detected_df[upper & lower]
    detected_in_band_E = detected_in_band_df["E"]

    if ret == "all":
        Es = {}
        Es["original_E"] = original_E
        Es["detected_snr_E"] = detected_snr_E
        Es["detected_fit_E"] = detected_fit_E
        Es["detected_in_band_E"] = detected_in_band_E
        return detected_df, detected_in_band_df, Es
    elif ret == "snr_E":
        return detected_snr_E
    elif ret == "in_band_E":
        return detected_in_band_E
    else:
        raise ValueError


def analyse_and_plot(bursts, detected, Es):
    "Plot the histogram of energies and print statistics of injected and recovered parameters"

    qs = np.quantile(bursts["in_mu_f"], [0.16, 0.5, 0.84])
    print(f"Injected mu_f is: {qs[1]:.3f}+{(qs[2] - qs[1]):.2f}-{(qs[1] - qs[0]):.2f}")

    qs = np.quantile(detected["in_mu_f"], [0.16, 0.5, 0.84])
    print(f"Recovered mu_f is: {qs[1]:.3f}+{(qs[2] - qs[1]):.2f}-{(qs[1] - qs[0]):.2f}")

    qs = np.quantile(bursts["in_sig_f"], [0.16, 0.5, 0.84])
    print(f"Injected sig_f is: {qs[1]:.2f}+{(qs[2] - qs[1]):.2f}-{(qs[1] - qs[0]):.2f}")

    qs = np.quantile(detected["in_sig_f"], [0.16, 0.5, 0.84])
    print(
        f"Recovered sig_f is: {qs[1]:.2f}+{(qs[2] - qs[1]):.2f}-{(qs[1] - qs[0]):.2f}"
    )

    fig, axes = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(12, 4))
    v, bins, _ = axes[0].hist(bursts["in_mu_f"], alpha=0.5, density=True, label="in")
    axes[0].hist(
        detected["in_mu_f"], alpha=0.5, density=True, label="detected", bins=bins
    )
    axes[0].set_title("mu_f")
    axes[0].legend()

    original_E = Es["original_E"]
    detected_snr_E = Es["detected_snr_E"]
    detected_fit_E = Es["detected_fit_E"]
    detected_in_band_E = Es["detected_in_band_E"]

    v, bins, _ = axes[1].hist(
        bursts["in_sig_f"], alpha=0.5, density=True, bins=30, label="in"
    )
    axes[1].hist(
        detected["in_sig_f"], alpha=0.5, density=True, label="detected", bins=bins
    )
    axes[1].set_title("sig_f")
    axes[1].legend()

    axes[2].scatter(detected["in_S"], detected["snr_S"])
    axes[2].plot(detected["in_S"], detected["in_S"], c="k")
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("in S")
    axes[2].set_ylabel("snr S")

    fig, axes = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(12, 4))
    axes[0].hist(np.log10(original_E), bins=30, density=True)
    axes[0].set_yscale("log")
    axes[0].set_title("original_E")

    axes[1].hist(np.log10(detected_snr_E), bins=30, density=True)
    axes[1].set_yscale("log")
    axes[1].set_title("detected_snr_E")

    axes[2].hist(
        np.log10(detected_fit_E), bins=30, density=True, alpha=0.5, label="fit_E"
    )
    axes[2].hist(
        np.log10(detected_in_band_E),
        bins=30,
        density=True,
        alpha=0.5,
        label="in_band_E",
    )
    axes[2].set_yscale("log")
    axes[2].set_title("fit and in-band E")
    axes[2].legend()
