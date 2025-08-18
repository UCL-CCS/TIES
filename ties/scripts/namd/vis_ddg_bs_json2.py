#!/usr/bin/env python3
"""
Visualise the bootstrapped replicas. We have 1, 3, 5, and some other replica lengths with their ddGs bootstrapped.
"""

import os
from pathlib import Path
from collections import OrderedDict
import glob
import time
import itertools
import random
import sys
import json
import re

import numpy as np
from numpy import genfromtxt
import scipy.stats as st
import matplotlib

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as trans


def combine(transformations):
    """
    Combine the different replicas (4 sets of 5)
    """

    # sort by transformation
    transformations.sort(key=lambda k: str(k).rsplit("/", maxsplit=1)[1])

    # group by transformation
    trans_dgs = {}
    for transformation, replica_sets in itertools.groupby(
        transformations, key=lambda t: str(t).rsplit("_", maxsplit=1)[0]
    ):
        print(f"Next: {transformation}")

        # create keys
        dg_replicas_samples = {}  # repsNum: [sample1, sample2, ]

        for repset in replica_sets:
            # load pickled lig and complex data
            with open(repset) as F:
                lig_all = json.load(F)
                # verify that the length is the same for each "replica number"
                sample_num = len(list(lig_all.values())[0])
                for next_rep_num, samples in lig_all.items():
                    int_next_rep = int(next_rep_num)
                    if sample_num != len(samples):
                        print(
                            f"In file {repset}, replica {next_rep_num} has {len(samples)} samples, "
                            f"but replica 1 has {sample_num}"
                        )

                    # merge the data
                    dg_replicas_samples.setdefault(int_next_rep, [])
                    dg_replicas_samples[int_next_rep].extend(samples)

        print(
            f"in each set there is {len(list(dg_replicas_samples.values())[0])} samples"
        )
        trans_dgs[transformation] = dg_replicas_samples

    # for this protein, across the ligand transformation,
    # plot the information in improvement as a function of an increase in the number of replicas
    # ie for each case, calculate how much "range" is removed when increasing the replica
    sd_improvements = {}
    sd_improvement_range = []
    for system, data in trans_dgs.items():
        sds = []
        for rep_no, dgs in data.items():
            interval95 = st.t.interval(
                0.95, len(dgs) - 1, loc=np.mean(dgs), scale=st.sem(dgs)
            )
            sd = np.std(dgs)
            # print(f'{system} with {rep_no}, sd {sd}')
            # print(f'{system} with {rep_no}, interval range {np.abs(interval95[0]-interval95[1])}')
            sds.append(sd)
        sd_improvements[system] = sds

    return trans_dgs, sd_improvements


def plot_ddgs(
    exp,
    trans_dgs,
    sd_improvements,
    work_dir,
    filename,
    yrange=0.5,
    plots=[5, 5],
    experimental=None,
):
    plt.figure(figsize=(15, 10))  # , dpi=80) facecolor='w', edgecolor='k')
    plt.rcParams.update({"font.size": 9})

    counter = 1
    for tran, data in trans_dgs.items():
        # extract the name of the transformation
        _, left, right = tran.split("/")[-1].split("_")
        parsed_trans = left.upper() + " " + right.upper()

        # plot the boxplots
        plt.subplot(plots[0], plots[1], counter)
        plt.xlim([0.5, 20 + 0.5])
        plt.title(parsed_trans)

        for number_of_reps, ddGs in data.items():
            plt.boxplot(
                ddGs,
                positions=[
                    number_of_reps,
                ],
                showfliers=False,
            )
        plt.ylabel("$\\rm Bootstrapped~\Delta \Delta G $")
        plt.xlabel("$\\rm Replicas~\# ~/~ 20  $")
        plt.xticks(range(0, 20 + 1, 2), range(0, 20 + 1, 2))

        # plot the exp value
        exp_value = exp[parsed_trans]
        plt.axhline(y=exp_value)

        # ensure that the experimental line can be seen
        # ylim_min, ylim_max = plt.axes().get_ylim()
        # plt.ylim([min(exp_value, ylim_min), max(exp_value, ylim_max)])

        counter += 1

    plt.tight_layout()
    plt.savefig(work_dir / f"ddg_{filename}_y{yrange * 2:0.2f}.png", dpi=300)
    # plt.show()

    # ----------------------------------
    # for each system plot show how much the system decreases the uncertainty
    plt.figure(figsize=(15, 10))
    plt.rcParams.update({"font.size": 9})

    counter = 1
    for system, sds in sd_improvements.items():
        plt.subplot(plots[0], plots[1], counter)
        # plt.xlim([0.5, 20 + 0.5])
        plt.title(system.rsplit("/")[-1])
        plt.plot([5, 10, 15, 20], sds)
        plt.tight_layout()
        plt.xlim([-0.5, 20.5])
        plt.xlabel("Number of replicas")
        plt.ylabel("$\\rm  \Delta\Delta G ~ \sigma $")
        plt.savefig(work_dir / f"ddg_sds_decrease_{filename}.png", dpi=300)
        counter += 1


def merge_protein_cases(exp, data):
    """
    For the system or protein, take all the transformations t and
    F = merge each | ddGs(t) - exp | for each transformation t

    So subtract the experiment from each transformation, then take the absolute values of that distributions,
    and merge the distributions together
    """
    merged_abs_dists = {5: [], 10: [], 15: [], 20: []}
    for tran, data in trans_ddGs.items():
        # extract the name of the transformation
        _, left, right = tran.split("/")[-1].split("_")
        parsed_trans = left.upper() + " " + right.upper()

        exp_value = exp[parsed_trans]

        for number_of_reps, ddGs in data.items():
            d = np.abs(np.array(ddGs) - exp_value)
            merged_abs_dists[number_of_reps].extend(d)

    return merged_abs_dists


def get_exp_data():
    # load the data
    # columns 0 are names, and 13 are experiments
    data = genfromtxt(
        "/home/dresio/pubs/ties20/figures/energies.csv",
        delimiter=",",
        dtype=None,
        encoding="UTF-8",
        usecols=[0, 13],
    )

    # separate the proteins
    exp_data = {}
    tyk = data[4:15].T
    exp_data["tyk2"] = {t: float(ddg) for t, ddg in zip(tyk[0], tyk[1])}
    mcl = data[16:32].T
    exp_data["mcl1"] = {t: float(ddg) for t, ddg in zip(mcl[0], mcl[1])}
    thrombin = data[34:44].T
    exp_data["thrombin"] = {t: float(ddg) for t, ddg in zip(thrombin[0], thrombin[1])}
    ptp1b = data[46:55].T
    exp_data["ptp1b"] = {t: float(ddg) for t, ddg in zip(ptp1b[0], ptp1b[1])}
    cdk = data[57:63].T
    exp_data["cdk2"] = {t: float(ddg) for t, ddg in zip(cdk[0], cdk[1])}

    return exp_data


exp_data = get_exp_data()

systems = []
complexes = []
ddg_minus_exp = {}
ddgs = {}

root_work = Path("/home/dresio/ucl/ties/ties20/replica20")

# tyk2
# ligand
transformations = list(root_work.glob("tyk2/analysis/ddg_135_l*_l*_*.json"))
trans_ddGs, sds = combine(transformations)
# plot_ddgs(exp_data['tyk2'], trans_ddGs, sds, root_work, filename='tyk2', plots=[2, 3])
# ddg_minus_exp['tyk2'] = merge_protein_cases(exp_data['tyk2'], trans_ddGs)
ddgs["tyk2"] = trans_ddGs

# mcl1
# ligand
transformations = list(root_work.glob("mcl1/analysis/ddg_135_l*_l*_*.json"))
trans_ddGs, sds = combine(transformations)
# plot_ddgs(exp_data['mcl1'], trans_ddGs, sds, root_work, filename='mcl1', yrange=1, plots=[4, 2])
# ddg_minus_exp['mcl1'] = merge_protein_cases(exp_data['mcl1'], trans_ddGs)
ddgs["mcl1"] = trans_ddGs

# thrombin
# lig
transformations = list(root_work.glob("thrombin/analysis/ddg_135_l*_l*_*.json"))
trans_ddGs, sds = combine(transformations)
# plot_ddgs(exp_data['thrombin'], trans_ddGs, sds, root_work, filename='thrombin', yrange=0.5, plots=[2, 3])
# ddg_minus_exp['thrombin'] = merge_protein_cases(exp_data['thrombin'], trans_ddGs)
ddgs["thrombin"] = trans_ddGs

# ptp1b
# lig
transformations = list(root_work.glob("ptp1b/analysis/ddg_135_l*_l*_*.json"))
trans_ddGs, sds = combine(transformations)
# plot_ddgs(exp_data['ptp1b'], trans_ddGs, sds, root_work, filename='ptp1b', yrange=3, plots=[2, 3])
# ddg_minus_exp['ptp1b'] = merge_protein_cases(exp_data['ptp1b'], trans_ddGs)
ddgs["ptp1b"] = trans_ddGs


# cdk2
# lig
transformations = list(root_work.glob("cdk2/analysis/ddg_135_l*_l*_*.json"))
trans_ddGs, sds = combine(transformations)
# plot_ddgs(exp_data['cdk2'], trans_ddGs, sds, root_work, filename='cdk2', yrange=0.5, plots=[1, 3])
# ddg_minus_exp['cdk2'] = merge_protein_cases(exp_data['cdk2'], trans_ddGs)
ddgs["cdk2"] = trans_ddGs


symbol_mapping = {
    "tyk2": "x",
    "mcl1": ".",
    "thrombin": "*",
    "ptp1b": "+",
    "cdk2": "o",
}
colour_mapping = {
    "tyk2": "#E13D77",
    "mcl1": "#C9E13D",
    "thrombin": "#3DE1A7",
    "ptp1b": "#553DE1",
    "cdk2": "black",
}


# ------------------------
def ddgs_relExp_by_prot():
    # plot for each protein separately the progress
    plt.figure(figsize=(15, 7))
    plt.rcParams.update({"font.size": 15})

    protein_counter = 1
    for protein, dists in ddg_minus_exp.items():
        plt.subplot(2, 3, protein_counter)
        plt.title(protein.upper())
        if protein_counter in (1, 4):
            plt.ylabel("$\\rm |\Delta\Delta G - exp| $")
        if protein_counter in (3, 4, 5):
            plt.xlabel("# Replicas used in bootstrapping / 20")

        # add the experimental range?

        plt.ylim([0, 3])

        case_counter = 1
        for rep_no, ddgs in dists.items():
            plt.boxplot(
                ddgs,
                positions=[
                    case_counter,
                ],
                widths=0.13,
                showfliers=False,
            )
            case_counter += 1

        plt.xticks([1, 2, 3, 4], [5, 10, 15, 20])

        protein_counter += 1

    plt.tight_layout()
    plt.savefig(root_work / "ddgs_byprot.png")
    # plt.show()


# ddgs_relExp_by_prot()


def ties_ddgs_dists(ddgs):
    # plot for each protein separately the progress
    # show the progress what happens if you use 5 replicas sampling or 1 replica,
    # and how that affects the distributions
    fig = plt.figure(figsize=(7, 5))
    # fig, ax = plt.subplots(figsize=(7, 5))
    # https://www.delftstack.com/howto/matplotlib/how-to-set-the-figure-title-and-axes-labels-font-size-in-matplotlib/
    plt.rcParams.update({"xtick.labelsize": 5, "axes.labelsize": 7})

    def prep(bsrep):
        """
        Subtract the mean
        """
        bsrep = np.array(bsrep)
        # i have 5k bs samples
        assert len(bsrep) == 5000
        mu0 = bsrep - np.mean(bsrep)
        return mu0, np.std(mu0)

    # plot the first 3 example of sigma in data
    # meaning that the first 3 transformation from each protein
    # are used to show the distributions,
    limit_cases = 3
    protein_counter = 1
    for protein, dists in ddgs.items():
        case_counter = 0
        # plt.ylim([0, 3])

        for trans_name, ddg_bs in dists.items():
            if "l6_l14" in trans_name and protein == "ptp1b":
                print("skipping the bad case ptp1b l6-l14")
                continue

            # e.g. '/home/dresio/ucl/validation/replica20/tyk2/analysis/ddg_l15_l10'
            _, _, left, right = trans_name.split("/")[-1].upper().split("_")

            plt.subplot(limit_cases + 1, 5, protein_counter + case_counter * 5)
            # plt.title(f'{protein.upper()} {left}-{right}')
            if protein_counter in (1,):
                plt.ylabel("P(x)")
            # if case_counter in (limit_cases - 1, ):
            plt.xlabel(r" $\rm \Delta \Delta G ^ {TIES} - \mu $", va="center")
            # plt.xlabel(r' $\rm \Delta \Delta G ^ {TIES} - \mu $',
            #           trans.ScaledTranslation(0, 0, fig.dpi_scale_trans))
            if case_counter == 0:
                plt.title(f"{protein.upper()}")

            # take only the bootstrapped 5 replicas from the 20 replicas we had
            mu0_1rep, std1 = prep(ddg_bs[1])
            mu0_5rep, std5 = prep(ddg_bs[5])
            # take the mean away from the all samples
            # determine the number of beans (every 0.25 ddg)
            bin_no = int((max(mu0_1rep) - min(mu0_1rep)) / 0.1)
            # make even
            if bin_no % 2 != 0:
                bin_no += 1
            print(f"Number of bins: {bin_no}")

            # stretch the bis
            bins = np.linspace(min(mu0_1rep), max(mu0_1rep), bin_no)
            bins = bins + bins / 2
            # print('bins', bins)

            # print(f'Interval {intFrom1: .2f} to {intTo1: .2f} and {intFrom5: .2f} to {intTo5: .2f}')

            plt.hist(mu0_1rep, bins=bins, alpha=0.7, density=True, color="#253bff")
            plt.hist(mu0_5rep, bins=bins, alpha=0.7, density=True, color="#ff8e25")

            plt.text(
                plt.xlim()[0],
                plt.ylim()[1],
                "$\sigma_1$=" + f"{std1:.2f}",
                {"color": "#253bff"},
                va="top",
                ha="left",
            )
            ylim_by5 = plt.ylim()[1] / 4
            plt.text(
                plt.xlim()[0],
                ylim_by5 * 3,
                "$\sigma_5$=" + f"{std5:.2f}",
                {"color": "#ff8e25"},
                va="top",
                ha="left",
            )
            plt.text(
                plt.xlim()[1], plt.ylim()[1], f"{left}-{right}", va="top", ha="right"
            )
            plt.tick_params(
                axis="y",  # changes apply to the x-axis
                which="both",  # both major and minor ticks are affected
                left=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelleft=False,
            )
            # plt.ticklabel_format(axis='x', useOffset=0.1)
            case_counter += 1
            if case_counter > (limit_cases - 1):
                break

        # plt.xticks([1, 2, 3, 4], [5, 10, 15, 20])

        protein_counter += 1

    # extract the overall improvement in the distribution of the ddGs
    # ie the std of ddGs with 5 replicas  - std of ddGs with 1 replica
    sd5subSd1 = OrderedDict()
    for protein, dists in ddgs.items():
        sd5subSd1[protein] = []
        for trans_name, ddgs_bs in dists.items():
            mu0_1rep, std1 = prep(ddgs_bs[1])
            mu0_5rep, std5 = prep(ddgs_bs[5])
            sd5subSd1[protein].append(std1 / std5)

    # plot as the 5th plot
    # the following plots the overall improvement in dispersion
    protein_counter = 1
    for protein, d_sd in sd5subSd1.items():
        plt.subplot(limit_cases + 1, 5, protein_counter + 15)

        if protein_counter != 1:
            plt.tick_params(
                axis="y",  # changes apply to the x-axis
                which="both",  # both major and minor ticks are affected
                left=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelleft=False,
            )
        else:
            plt.ylabel("Count", va="top")
        plt.xlabel("$\\rm \sigma_1 / \sigma_5 $")

        plt.hist(d_sd)
        plt.ylim([0, 6])
        protein_counter += 1

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.0, hspace=0.4)
    # plt.subplot_tool()
    plt.savefig(root_work / "noexp_135_ddgs_dists.png", dpi=300)


def ties_ddgs_overall_dist_shape(ddgs):
    # take the distributions from all the sampling done with 1 replica,
    # however, merge all the distributions together to find out what the overall trend is,
    # so how do we compare the 55 TIES distributions?
    # so for each distribution we take Skewness and Kurtosis and keep track of them

    skews = []
    kurts = []
    for protein, dists in ddgs.items():
        for trans_name, ddg_bs in dists.items():
            if "l6_l14" in trans_name and protein == "ptp1b":
                print("skipping the bad case ptp1b l6-l14")
                continue

            # take only the bootstrapped 5 replicas from the 20 replicas we had
            skew = st.skew(ddg_bs[1])
            skews.append(skew)
            kurtosis = st.kurtosis(ddg_bs[1])
            kurts.append(kurtosis)

    fig = plt.figure(figsize=(7, 5))
    # https://www.delftstack.com/howto/matplotlib/how-to-set-the-figure-title-and-axes-labels-font-size-in-matplotlib/
    # plt.rcParams.update({'xtick.labelsize': 5, 'axes.labelsize': 7})

    plt.subplot(211)
    plt.hist(skews, density=True)
    plt.xlabel("Skewness")
    plt.ylabel("Frequency")
    plt.tick_params(
        axis="y",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        left=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelleft=False,
    )

    plt.subplot(212)
    plt.hist(kurts, density=True)
    plt.xlabel("Kurtosis")
    plt.ylabel("Frequency")
    plt.tick_params(
        axis="y",  # changes apply to the x-axis
        which="both",  # both major and minor ticks are affected
        left=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelleft=False,
    )

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.0, hspace=0.4)
    # plt.subplot_tool()
    plt.savefig(root_work / "ties20_normality.png", dpi=300)
    plt.show()


# ties_ddgs_dists(ddgs)
ties_ddgs_overall_dist_shape(ddgs)

# -----------------------------------
# lot the ddg-exp joined distributions as a bar plot
# plt.figure(figsize=(15, 10))
# plt.rcParams.update({'font.size': 9})
# plt.ylabel('$\\rm Bootstrapped~ |\Delta\Delta G - exp| $')
# plt.xlabel('Replica number')
#
# # rearrange the data, such that each internal list represents different proteins
#
# protein_factor = -0.3
# for protein, dists in ddg_minus_exp.items():
#     case_counter = 1
#     for rep_no, ddgs in dists.items():
#         plt.boxplot(ddgs, positions=[case_counter + protein_factor, ], widths=0.13, showfliers=False)
#         case_counter += 1
#
#
#     protein_factor += 0.15
#
# plt.xticks([1, 2, 3, 4], [5, 10, 15, 20])
# plt.show()

# ---------------------------------------------------------------
# combine all sds, so show the points and the distribution
# ligands first
# plt.figure(figsize=(15, 10))
# plt.rcParams.update({'font.size': 9})
# plt.subplot(1,2,1)
# plt.title('Ligands')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for lig in ligands:
#         for system, sds in lig.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color)
#
# # complex
# plt.subplot(1, 2, 2)
# plt.title('Complexes')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for complex in complexes:
#         for system, sds in complex.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color, label=name)
#
# # remove duplicate labels
# handles, labels = plt.gca().get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys())
#
# plt.tight_layout()
# plt.savefig(root_work / 'combined_sds.png')

# ---------------------------------------------------------------
# now we do the same thing but use the experimental data instead. We take the distribution, and subtract the exp value
# as a function of 5, 10, 15, 20 replicas,
# so we'll see if the values come closer to the experimental value, rather than just decreasing error

# combine all sds, so show the points and the distribution
# ligands first
# plt.figure(figsize=(15, 10))
# plt.rcParams.update({'font.size': 9})
# plt.title('ddG')
# plt.ylabel('$\\rm  \Delta G ~ \sigma $')
# plt.xlabel('Replica number')
#
# for i in range(20):
#     # extract from each case just the one index
#     for lig in systems:
#         for system, sds in lig.items():
#             if 'l6_l14' in system:
#                 continue
#             name = re.match(r".*/replica20/([a-z0-9]+)/analysis/.*", system).group(1)
#             sym = symbol_mapping[name]
#             color = colour_mapping[name]
#             plt.scatter(range(1, len(sds) + 1), sds, marker=sym, color=color)
#
#
# plt.tight_layout()
# plt.savefig(root_work / 'combined_sds.png')
