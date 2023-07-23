import os
import sys
import subprocess
import datetime
from dateutil import parser
from os import listdir
from os.path import isfile, join
from collections import OrderedDict
import pandas as pd
import numpy as np
from glob import glob
import random

###################################################
# Parse configfile 

# make sure required params are set
if not config:
    raise Exception("Please provide a configuration file with --configfile")

required_config_params = [
    "out_dir", 
    "reps", 
    "clades",
    "clade_size",
    "num_loci",
    "loc_len",
    "snp_per_locus",
    "tree_height",
    "relative_clade_height",
    "alpha", 
    "missing_proportion",
    "missing_strategy"
]

missing_params = [param for param in required_config_params if param not in config]
if missing_params:
    raise Exception(f"Missing required configuration parameters: {', '.join(missing_params)}")

# Optional parameters
config.setdefault("model", "gtrgamma")
config.setdefault("write_gene_alignments", False)

# Check if certain parameters are lists, convert if not
for param in ["tree_height", "relative_clade_height", "alpha", "missing_proportion", "missing_strategy"]:
    if not isinstance(config[param], list):
        config[param] = [config[param]]

# format out_dir
config["out_dir"] = os.path.abspath(config["out_dir"])
config["out_dir"] = os.path.normpath(config["out_dir"])

##################################################
# Define the output pattern
output_pattern = "rep{rep}_t{t}_c{c}_a{a}_{model}"

# Expand the output pattern for each parameter combination
final_outputs = expand(
    [os.path.join("{out_dir}", "datasets", output_pattern + ".phylip"),
    os.path.join("{out_dir}", "datasets", output_pattern + "_guidetree.tre"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.treefile"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.mlrate"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.iqtree"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.rooted.tre"),
    os.path.join("{out_dir}", "masks", output_pattern + "_p{p}_{strategy}_mask.npy"),
    os.path.join("{out_dir}", "masks", output_pattern + "_p{p}_{strategy}_original_missing_mask.npy"),
    os.path.join("{out_dir}", "metadata", "popmap.txt")],
    out_dir=config["out_dir"],
    rep=range(config["reps"]),
    t=config["tree_height"],
    c=config["relative_clade_height"],
    a=config["alpha"],
    model=[config["model"]],
    p=config["missing_proportion"],
    strategy=config["missing_strategy"]
)

##################################################
# SnakeMake workflow 

rule all:
    input:
        final_outputs

rule generate_popmap:
    output:
        os.path.join("{out_dir}", "metadata", "popmap.txt")
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], 'w') as f:
            for clade in range(config["clades"]):
                for ind in range(config["clade_size"]):
                    f.write(f"pop{clade}_{ind}\tpop{clade}\n")

rule simulate_sequences:
    output:
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}_guidetree.tre"),
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip")
    log:
        os.path.join(config["out_dir"], "logs", "tree_sim", "rep{rep}_t{t}_c{c}_a{a}_{model}.log")
    conda:
        os.path.join("envs", "tree_sim.yml")
    params:
        script_path = os.path.join(workflow.basedir, "scripts", "sim_treeparams.py"),
        output_prefix = os.path.join(config["out_dir"], "datasets", "rep"),
        num_loci = config["num_loci"],
        loc_len = config["loc_len"],
        num_clades = config["clades"],
        samples_per_clade = config["clade_size"],
        snp_per_locus = config["snp_per_locus"],
        seed = random.randint(1, 1e6)
    shell:
        """
        mkdir -p $(dirname {output[0]}) && \
        python3 {params.script_path} \
            --prefix {params.output_prefix}{wildcards.rep}_ \
            --alpha {wildcards.a} \
            --num_loci {params.num_loci} \
            --model {wildcards.model} \
            --tree_height {wildcards.t} \
            --relative_clade_height {wildcards.c} \
            --num_clades {params.num_clades} \
            --loc_length {params.loc_len} \
            --snps_per_locus {params.snp_per_locus} \
            --samples_per_clade {params.samples_per_clade} \
            --seed {params.seed} > {log} 2>&1
        """

rule run_iqtree:
    input:
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip")    
    output:
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip.treefile"),
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip.mlrate"),
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip.iqtree")
    conda:
        os.path.join("envs", "iqtree2.yml")
    log:
        os.path.join(config["out_dir"], "logs", "iqtree", "rep{rep}_t{t}_c{c}_a{a}_{model}.log")
    params:
        seed = random.randint(1, 1e6),
        threads = 1
    shell:
        """
        mkdir -p $(dirname {output[0]}) && \
        iqtree \
            -s {input[0]} \
            -m "GTR+I*G4" \
            -redo \
            --seed {params.seed} \
            -T {params.threads} \
            --mlrate > {log} 2>&1 && \
        rm -f {input[0]}*.mldist \
            {input[0]}*.ckp.gz \
            {input[0]}*.log \
            {input[0]}*.bionj \
            {input[0]}*.uniqueseq.phy \
            {input[0]}*.mldist
        """

rule reroot_tree:
    input:
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip.treefile")
    output:
        os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip.rooted.tre")
    conda:
        os.path.join("envs", "tree_sim.yml")
    log:
        os.path.join(config["out_dir"], "logs", "reroot_tree", "rep{rep}_t{t}_c{c}_a{a}_{model}.log")
    params:
        outgroup = "pop" + str(config["clades"] - 1) + "_",
        script_path = os.path.join(workflow.basedir, "scripts", "reroot_tree.py")
    shell:
        """
        mkdir -p $(dirname {output[0]}) && \
        python3 {params.script_path} \
            -i {input} \
            -o {output} \
            -O {params.outgroup} > {log} 2>&1
        """

rule generate_mask:
    input:
        phylip = os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip"),
        popmap = os.path.join(config["out_dir"], "metadata", "popmap.txt"),
        tree = os.path.join(config["out_dir"], "datasets", "rep{rep}_t{t}_c{c}_a{a}_{model}.phylip.rooted.tre"),
    output:
        missing = os.path.join(config["out_dir"], "masks", "rep{rep}_t{t}_c{c}_a{a}_{model}_p{p}_{strategy}_original_missing_mask.npy"),
        mask = os.path.join(config["out_dir"], "masks", "rep{rep}_t{t}_c{c}_a{a}_{model}_p{p}_{strategy}_mask.npy")
    params:
        script_path = os.path.join(workflow.basedir, "scripts", "generate_mask.py"),
        output_prefix = os.path.join(config["out_dir"], "masks", "rep{rep}_t{t}_c{c}_a{a}_{model}_p{p}_{strategy}_original_missing")
    log:
        os.path.join(config["out_dir"], "logs", "generate_mask", "rep{rep}_t{t}_c{c}_a{a}_{model}_p{p}_{strategy}.log")
    conda:
        os.path.join("envs", "pgsui.yml")
    shell:
        """
        mkdir -p $(dirname {output.mask}) && \
        python3 {params.script_path} \
            --input {input.phylip} \
            --popmap {input.popmap} \
            --tree {input.tree} \
            --missing {wildcards.p} \
            --strategy {wildcards.strategy} \
            --output {params.output_prefix} \
            > {log} 2>&1
        """
