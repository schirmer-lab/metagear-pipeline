#!/usr/bin/env python3

import os
import click
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np


# return {msp_id:{sample_id:msp_abd(RPKM)}}
def _load_msp_profile(msp_profile_fp, sep="\t", has_column_header=False):
    """
    Load MSP abundance profile from file.

    Parameters:
        msp_profile_fp (str): Path to MSP profile file.
        sep (str): Delimiter used in the file.
        has_column_header (bool): If True, the first column of the header contains sample IDs.

    Returns:
        dict: Nested dictionary {msp_id: {sample_id: RPKM_value}}.
    """
    ret_abd_prev = {}
    ret_sample_lst = []
    with open(msp_profile_fp, "r") as f:
        for ii, line in enumerate(f):
            if ii == 0:
                if has_column_header:
                    ret_sample_lst = line.strip().split(sep)[1:]
                else:
                    ret_sample_lst = line.strip().split(sep)
                continue
            lst = line.strip().split(sep)
            ret_abd_prev[lst[0]] = [float(it) for it in lst[1:]]
    # reformat the abundance list by sample_id
    ret_abd = {}
    for msp_id in ret_abd_prev:
        ret_abd[msp_id] = dict(zip(ret_sample_lst, ret_abd_prev[msp_id]))
    return ret_abd


# return {spp:metaphlan_abd_lst(relative abd)},sample_id_lst
def _load_metaphlan_profile(metaphlan_profile_fp, sep="\t"):
    """
    Load MetaPhlAn relative abundance profile from combined output format.

    Parameters:
        metaphlan_profile_fp (str): Path to MetaPhlAn combined profile file.
        sep (str): Delimiter used in the file (default: tab-separated).

    Returns:
        dict: Nested dictionary {species: {sample_id: relative_abundance}}.
    """
    ret_abd_prev = {}
    ret_sample_lst = []
    with open(metaphlan_profile_fp, "r") as f:
        for ii, line in enumerate(f):
            # Skip header lines starting with #
            if line.startswith("#"):
                continue

            # First data line contains sample IDs
            if ii == 0 or (ii == 1 and len(ret_sample_lst) == 0):
                lst = line.strip().split(sep)
                # Skip first column (clade_name), get sample IDs
                ret_sample_lst = [
                    sample.replace("_microbial", "").replace("_paired", "")
                    for sample in lst[1:]
                ]
                continue

            lst = line.strip().split(sep)
            # First column is the species/taxa name
            spp_name = lst[0]
            # Only process species-level entries (containing 's__')
            if "s__" in spp_name:
                abundances = [float(abd) for abd in lst[1:]]
                ret_abd_prev[spp_name] = abundances

    # reformat the abundance list by sample_id
    ret_abd = {}
    for spp in ret_abd_prev:
        ret_abd[spp] = dict(zip(ret_sample_lst, ret_abd_prev[spp]))

    return ret_abd


# return the Coefficients,Interceptions,Mean squared error, and r2 values
def _linear_fit_Apair(cur_msp_id, cur_spp, metaphlan_abd, msp_abd, sample_id_lst):
    """
    Perform linear regression between one MSP and one species across shared samples.

    Parameters:
        cur_msp_id (str): MSP ID.
        cur_spp (str): Species name from MetaPhlAn.
        metaphlan_abd (dict): Species abundance {spp: {sample_id: value}}.
        msp_abd (dict): MSP abundance {msp_id: {sample_id: value}}.
        sample_id_lst (list): Sample IDs to include in the regression.

    Returns:
        list: [coefficient, intercept, MSE, R2 score]
    """
    cur_metaphlan_abd = np.array(
        [
            float(metaphlan_abd[cur_spp][cur_sample_id])
            for cur_sample_id in sample_id_lst
        ]
    )
    cur_metaphlan_abd = cur_metaphlan_abd.reshape((len(cur_metaphlan_abd), 1))
    cur_msp_abd = np.array(
        [float(msp_abd[cur_msp_id][cur_sample_id]) for cur_sample_id in sample_id_lst]
    )
    cur_msp_abd = cur_msp_abd.reshape((len(cur_msp_abd), 1))
    # linear model
    # Create linear regression object
    regr = linear_model.LinearRegression()

    # Fit a linear model
    regr.fit(cur_metaphlan_abd, cur_msp_abd)
    cur_msp_abd_pred = regr.predict(cur_metaphlan_abd)
    mse = mean_squared_error(cur_msp_abd, cur_msp_abd_pred)
    r2 = r2_score(cur_msp_abd, cur_msp_abd_pred)

    return [regr.coef_, regr.intercept_, mse, r2]


# run through all paired msp_id and metaphlan spp
# save to output_dir as a file named: msp_metaphlan_LM.full.txt
# columns as: "metaphlan_spp","msp_id","coefficients","interceptions","MSE","r2"
# modified from examples in https://scikit-learn.org/stable/auto_examples/linear_model/plot_ols.html
def linear_fit_all(full_res_sfp, msp_id_lst, spp_lst, metaphlan_abd, msp_abd):
    res = {}
    metaphlan_sample_id_lst = list(metaphlan_abd[spp_lst[0]].keys())
    msp_sample_id_lst = list(msp_abd[msp_id_lst[0]].keys())
    sample_id_lst = [it for it in metaphlan_sample_id_lst if it in msp_sample_id_lst]
    print("N samples: ", len(sample_id_lst))
    # print("msp missing samples:", [it for it in metaphlan_sample_id_lst if it not in msp_sample_id_lst])
    # print("metaphlan missing samples:", [it for it in msp_sample_id_lst if it not in metaphlan_sample_id_lst])

    for cur_msp_id in msp_id_lst:
        for cur_spp in spp_lst:
            res[(cur_spp, cur_msp_id)] = _linear_fit_Apair(
                cur_msp_id, cur_spp, metaphlan_abd, msp_abd, sample_id_lst
            )

    with open(full_res_sfp, "w") as sf:
        sf.write(
            "\t".join(
                [
                    "metaphlan_spp",
                    "msp_id",
                    "coefficients",
                    "interceptions",
                    "MSE",
                    "r2",
                ]
            )
            + "\n"
        )
        for cur_spp, cur_msp_id in res:
            sf.write(
                "\t".join(
                    [cur_spp, cur_msp_id]
                    + [str(float(val)) for val in res[(cur_spp, cur_msp_id)]]
                )
                + "\n"
            )
    return


# file format:
# "metaphlan_spp","msp_id","coefficients","interceptions","MSE","r2"
# return {("metaphlan_spp","msp_id"):["coefficients","interceptions","MSE","r2"]}
def load_msp_metaphlan_LM(fp, sep="\t", sel_msp_id_lst=[], sel_metaphlan_spp_lst=[]):
    ret = {}
    with open(fp, "r") as f:
        for ii, line in enumerate(f):
            if ii == 0:
                continue  # skip headers
            lst = line.strip().split(sep)
            metaphlan_spp, msp_id, coefficients, interceptions, MSE, r2 = lst
            if (len(sel_msp_id_lst) > 0) and (msp_id not in sel_msp_id_lst):
                continue
            if (len(sel_metaphlan_spp_lst) > 0) and (
                metaphlan_spp not in sel_metaphlan_spp_lst
            ):
                continue
            ret[(metaphlan_spp, msp_id)] = [
                float(it) for it in [coefficients, interceptions, MSE, r2]
            ]
    return ret


def _msp_taxaANN_metaphlanLM(
    output_dir, msp_profile_fp, metaphlan_profile_fp, metaphlan_version="v3"
):
    # load abd tables from both msp and metaphlan
    msp_abd = _load_msp_profile(msp_profile_fp)
    metaphlan_abd = _load_metaphlan_profile(metaphlan_profile_fp, sep="\t")

    msp_id_lst = list(msp_abd.keys())
    # Filter to only species-level taxa (containing 's__')
    spp_lst = [spp for spp in list(metaphlan_abd.keys()) if "s__" in spp]

    print(f"Loaded {len(msp_id_lst)} MSPs and {len(spp_lst)} species from MetaPhlAn")

    if metaphlan_version == "v3":
        full_res_sfp = os.path.join(output_dir, "msp_metaphlan_LM.full.txt")
        bestR2_res_sfp = os.path.join(output_dir, "msp_metaphlan_LM.bestR2.txt")
    elif metaphlan_version == "v4":
        full_res_sfp = os.path.join(output_dir, "msp_metaphlan4_LM.full.txt")
        bestR2_res_sfp = os.path.join(output_dir, "msp_metaphlan4_LM.bestR2.txt")
    else:
        print("unknown metaphlan version [v3 or v4]:", metaphlan_version)
        return
    # linear fit all against all
    # save everything
    linear_fit_all(
        full_res_sfp, msp_id_lst, spp_lst, metaphlan_abd, msp_abd
    )  # only run once for each cohort!

    # reload the res from file
    res_reload = load_msp_metaphlan_LM(full_res_sfp)
    # get the best metaphlan_spp for each msp_id according to the r2 values
    res_best_R2 = {}
    for item in res_reload:
        metaphlan_spp, msp_id = item
        coefficients, interceptions, MSE, r2 = res_reload[item]
        if msp_id not in res_best_R2:
            res_best_R2[msp_id] = [None, None, None, None]  # as a pseudo-min
        if res_best_R2[msp_id][-1] == None:
            res_best_R2[msp_id] = [metaphlan_spp, coefficients, interceptions, MSE, r2]
            continue
        if res_best_R2[msp_id][-1] < r2:
            res_best_R2[msp_id] = [metaphlan_spp, coefficients, interceptions, MSE, r2]
            continue
    # save the best match of each MSP

    with open(bestR2_res_sfp, "w") as sf:
        sf.write(
            "\t".join(
                [
                    "msp_id",
                    "metaphlan_spp",
                    "coefficients",
                    "interceptions",
                    "MSE",
                    "r2",
                ]
            )
            + "\n"
        )
        cur_msp_id_lst = list(res_best_R2.keys())
        cur_msp_id_lst.sort()
        for cur_msp_id in cur_msp_id_lst:
            cur_spp = res_best_R2[cur_msp_id][0]
            sf.write(
                "\t".join(
                    [cur_msp_id, cur_spp]
                    + [str(float(val)) for val in res_best_R2[cur_msp_id][1:]]
                )
                + "\n"
            )
    return  # done


# Example usage:
# output_dir = "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/mspminer"
# msp_profile_fp = "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/mspminer/msp_abundance.median.RPKM.txt"
# metaphlan_profile_fp = "/path/to/metaphlan_combined_output.txt"  # Tab-separated MetaPhlAn combined output
# _msp_taxaANN_metaphlanLM(output_dir,msp_profile_fp,metaphlan_profile_fp,metaphlan_version = "v3")


@click.command()
@click.option(
    "--msp-profile",
    required=True,
    type=click.Path(exists=True),
    help="Path to MSP abundance profile file (tab-separated, RPKM values)",
)
@click.option(
    "--metaphlan-profile",
    required=True,
    type=click.Path(exists=True),
    help="Path to MetaPhlAn combined profile file (tab-separated)",
)
@click.option(
    "--output-dir",
    required=True,
    type=click.Path(),
    help="Output directory for results",
)
@click.option(
    "--metaphlan-version",
    type=click.Choice(["v3", "v4"]),
    default="v3",
    help="MetaPhlAn version (affects output filenames, default: v3)",
)
def main(msp_profile, metaphlan_profile, output_dir, metaphlan_version):
    """
    Perform taxonomic annotation of MSPs using MetaPhlAn profiles through linear regression analysis.

    This tool performs linear regression between MSP abundance profiles and MetaPhlAn species
    abundance profiles to find the best taxonomic matches for each MSP.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Run the analysis
    click.echo(f"Starting MSP-MetaPhlAn taxonomic annotation analysis...")
    click.echo(f"MSP profile: {msp_profile}")
    click.echo(f"MetaPhlAn profile: {metaphlan_profile}")
    click.echo(f"Output directory: {output_dir}")
    click.echo(f"MetaPhlAn version: {metaphlan_version}")

    _msp_taxaANN_metaphlanLM(
        output_dir, msp_profile, metaphlan_profile, metaphlan_version
    )

    click.echo("Analysis completed successfully!")


if __name__ == "__main__":
    main()
