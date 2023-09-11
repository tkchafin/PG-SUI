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
output_pattern = "rep{rep}-t{t}-c{c}-a{a}-{model}"

# Expand the output pattern for each parameter combination
final_outputs = expand(
    [os.path.join("{out_dir}", "datasets", output_pattern + ".phylip"),
    os.path.join("{out_dir}", "datasets", output_pattern + "_guidetree.tre"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.treefile"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.mlrate"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.iqtree"),
    os.path.join("{out_dir}", "datasets", output_pattern + ".phylip.rooted.tre"),
    os.path.join("{out_dir}", "masks", output_pattern + "-p{p}-{strategy}_mask.npy"),
    os.path.join("{out_dir}", "masks", output_pattern + "-p{p}-{strategy}_original_missing_mask.npy"),
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
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}_guidetree.tre"),
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip")
    log:
        os.path.join(config["out_dir"], "logs", "tree_sim", "rep{rep}-t{t}-c{c}-a{a}-{model}.log")
    conda:
        os.path.join("envs", "tree_sim.yml")
    params:
        script_path = os.path.join(workflow.basedir, "scripts", "sim_treeparams.py"),
        output_prefix = os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}"),
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
            --prefix {params.output_prefix} \
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
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip")    
    output:
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip.treefile"),
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip.mlrate"),
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip.iqtree")
    conda:
        "envs/iqtree2.yml" if "iqtree_bin" not in config else None
    log:
        os.path.join(config["out_dir"], "logs", "iqtree", "rep{rep}-t{t}-c{c}-a{a}-{model}.log")
    params:
        seed = random.randint(1, 1e6),
        threads = 1,
        iqtree_bin = config["iqtree_bin"] if "iqtree_bin" in config else "iqtree"
    shell:
        """
        mkdir -p $(dirname {output[0]}) && \
        {params.iqtree_bin} \
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
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip.treefile")
    output:
        os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip.rooted.tre")
    conda:
        os.path.join("envs", "tree_sim.yml")
    log:
        os.path.join(config["out_dir"], "logs", "reroot_tree", "rep{rep}-t{t}-c{c}-a{a}-{model}.log")
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
        phylip = os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip"),
        popmap = os.path.join(config["out_dir"], "metadata", "popmap.txt"),
        tree = os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip.rooted.tre"),
    output:
        missing = os.path.join(config["out_dir"], "masks", "rep{rep}-t{t}-c{c}-a{a}-{model}-p{p}-{strategy}_original_missing_mask.npy"),
        mask = os.path.join(config["out_dir"], "masks", "rep{rep}-t{t}-c{c}-a{a}-{model}-p{p}-{strategy}_mask.npy")
    params:
        script_path = os.path.join(workflow.basedir, "scripts", "generate_mask.py"),
        output_prefix = os.path.join(config["out_dir"], "masks", "rep{rep}-t{t}-c{c}-a{a}-{model}-p{p}-{strategy}")
    log:
        os.path.join(config["out_dir"], "logs", "generate_mask", "rep{rep}-t{t}-c{c}-a{a}-{model}-p{p}-{strategy}.log")
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


# rule mask_data

# rule plot_missingness

# rule imputation

# rule get_accuracy 
#   summarize accuracy 
#   also accuracy by site X some other metrics? MAF/ missingness ? 

# rule collect_accuracy


# rule imputation:
#     input:
#         phylip = os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip"),
#         popmap = os.path.join(config["out_dir"], "metadata", "popmap.txt"),
#         tree = os.path.join(config["out_dir"], "datasets", "rep{rep}-t{t}-c{c}-a{a}-{model}.phylip.rooted.tre"),
#         mask = os.path.join(config["out_dir"], "masks", "rep{rep}-t{t}-c{c}-a{a}-{model}-p{p}-{strategy}_mask.npy"),
#     params:
#         strategy = config["missing_strategy"],
#         prop_missing = config["missing_proportion"],
#         method = config["imputation_method"],
#     output:
#         os.path.join(config["out_dir"], "imputed", "rep{rep}-t{t}-c{c}-a{a}-{model}-p{p}-{strategy}-{method}_imputed.phy"),
#     log:
#         os.path.join(config["out_dir"], "logs", "imputation", "rep{rep}-t{t}-c{c}-a{a}-{model}-p{p}-{strategy}-{method}.log"),
#     run:
#         # Define classes map
#         classes_map = {
#             "ImputeKNN": ImputeKNN,
#             "ImputeRandomForest": ImputeRandomForest,
#             # Define the rest of the classes
#         }

#         # Read in the inputs
#         genotype_data = GenotypeData(
#             filename=input.phylip,
#             popmapfile=input.popmap,
#             guidetree=input.tree,
#         )

#         # Create a SimGenotypeDataTransformer instance and use it
#         # to simulate missing data
#         transformer = SimGenotypeDataTransformer(
#             genotype_data=genotype_data,
#             prop_missing=params.prop_missing,
#             strategy=params.strategy,
#         )
#         transformer.fit(genotype_data.genotypes_012(fmt="numpy"))
#         simulated_data = copy.deepcopy(genotype_data)

#         simulated_data.genotypes_012 = transformer.transform(
#             genotype_data.genotypes_012(fmt="numpy")
#         )

#         # Determine if we're using grid search or not
#         do_gridsearch = "_grid" in params.method
#         # Get the class name
#         class_name = params.method.replace("_grid", "")

#         # Get the class
#         class_instance = classes_map[class_name]

#         if do_gridsearch:
#             # Define your parameter grid
#             param_grid = {
#                 "n_estimators": [50, 100],
#                 # Add more parameters if needed
#             }
#         else:
#             param_grid = None

#         instance = class_instance(
#             simulated_data,
#             gridparams=param_grid,
#         )
#         imputed_data = instance.imputed.genotypes_012(fmt="numpy")

#         # Save imputed data
#         np.savetxt(output[0], imputed_data, fmt="%d")

#         # Log the results
#         with open(log[0], 'w') as f:
#             f.write("Imputation complete.\n")

