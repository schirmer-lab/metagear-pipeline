#!/usr/bin/env python

import os, glob, sys


# return {protein_id:[length,start,end,Pfam_id,Pfam_description,eval]}
def _load_FG_interproscan_Pfam(
    interproscan_fp,
    ann_tool,
    start_pos_idx=6,
    pfam_id_idx=4,
    eval_idx=8,
    prot_len_idx=2,
):
    ann_dic = {}
    cur_ann_tool = ann_tool
    with open(interproscan_fp, "r") as f:
        for line in f:
            if cur_ann_tool not in line:
                continue
            lst = line.strip().split("\t")
            protein_id = lst[0]
            start_pos = int(lst[start_pos_idx])
            end_pos = int(lst[start_pos_idx + 1])
            pfam_id = lst[pfam_id_idx]
            pfam_ann = lst[pfam_id_idx + 1]
            cur_eval = float(lst[eval_idx])
            prot_len = int(lst[prot_len_idx])
            if protein_id not in ann_dic:
                ann_dic[protein_id] = []
            ann_dic[protein_id].append(
                [prot_len, start_pos, end_pos, pfam_id, pfam_ann, cur_eval]
            )
    return ann_dic


def _calculate_functional_group(
    interproscan_fp,
    FG_sfp,
    ann_tool="Pfam",
    start_pos_idx=6,
    pfam_id_idx=4,
    eval_idx=8,
    prot_len_idx=2,
):
    ann_dic = _load_FG_interproscan_Pfam(
        interproscan_fp, ann_tool, start_pos_idx, pfam_id_idx, eval_idx, prot_len_idx
    )
    # return: {protein_id:[length,start,end,Pfam_id,Pfam_description,eval]}

    # add create functino group for each protein
    # (1) sort by start position for each protein
    # print([it for it in ann_dic if len(ann_dic[it])>1])
    # cur_protein_id =  'cohort_merged__MSMA26EN_k105_33921::5::4320::6767::+'
    FG_annotation_dic = {}
    for cur_protein_id in ann_dic:
        cur_ann_rec = ann_dic[cur_protein_id]
        cur_start_lst = [it[1] for it in cur_ann_rec]
        cur_start_lst__sorted_idx = sorted(
            range(len(cur_start_lst)), key=lambda i: cur_start_lst[i]
        )
        cur_ann_rec__sorted = [cur_ann_rec[idx] for idx in cur_start_lst__sorted_idx]
        cur_FG = ":::".join([it[3] for it in cur_ann_rec__sorted])
        FG_annotation_dic[cur_protein_id] = cur_FG
    # save to file
    with open(FG_sfp, "w") as sf:
        sf.write("\t".join(["protein_id", "FG"]) + "\n")
        for protein_id in FG_annotation_dic:
            sf.write("\t".join([protein_id, FG_annotation_dic[protein_id]]) + "\n")
    return


if __name__ == "__main__":
    interproscan_fp = sys.argv[1]  # merged interproscan file
    FG_sfp = sys.argv[2]  # output file

    _calculate_functional_group(interproscan_fp, FG_sfp)
