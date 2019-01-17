import argparse
import os
import sys
import pandas as pd

def main():
    description = """Apply the provided set of sample and locus  filters to the input de novo STR mutations file.
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--all-mutations-file",           type=str,   required=True,   dest="MUT_FILE",                         help="Input de novo STR mutations file. (- for stdin)")
    parser.add_argument("--filtered-mutations-file",      type=str,   required=False,  dest="FILTER_FILE",                      help="Input de novo STR mutations file. (- for stdin)")
    parser.add_argument("--log-file",                     type=str,   required=False,  dest="LOG_FILE",                         help="Input de novo STR mutations file. (- for stdin)")
    parser.add_argument("--filter-denovos-child",         type=int,   required=False,  dest="CHILD_THRESH",         default=0,  help=""])
    parser.add_argument("--filter-loc-denovos",           type=int,   required=False,  dest="LOCI_THRESH",          default=0,  help="")
    parser.add_argument("--filter-posterior",             type=float, required=False,  dest="POSTERIOR_THRESH",     default=0,  help="")
    parser.add_argument("--filter-step-size",             action="store_true",  required=False,      dest="STEP_SIZE",        default=False,          help="")
    # parser.add_argument("--min-denovos-child",       type=int,   required=False, dest="",                 default=0,        help="")
    # parser.add_argument("--max-denovos-child",       type=int,   required=False, dest="",             default=10000,        help="")
    # parser.add_argument("--min-loc-denovos",         type=int,   required=False, dest="",         default=0,        help=help_dict["--min-loc-denovos"])
    # parser.add_argument("--max-loc-denovos",         type=int,   required=False, dest="",         default=0,        help=help_dict["--max-loc-denovos"])
    args          = parser.parse_args()

    # Load all de novo mutation
    mutations = pd.read_table(os.path.join(MUT_FILE), sep="\t")
    mutations["chrom"] = mutations["chrom"].map(lambda x: x.lstrip("chr").asint())
    mutations["str_id"] = mutations.map(lambda x: x["chrom"]+"_"x["pos"])

    # Log summary statistics before filtering
    # Number de novo per person
    child_num_total = mutations.groupby(["family", "child", "phenotype"]).pos.count().reset_index(name='child_num_total')
    # Number de novo per loci
    str_num_total = mutations.groupby(["chrom", "pos", "period"]).child.count().reset_index(name='str_num_total')

    # Exclude dn STRs where mutation size is not multiple of period
    if STEP_SIZE:
        mutations["unit"] = mutations.apply(lambda x: x["mutsize"]%x["period"]==0, 1)
        mutations = mutations[mutations.unit]
        n_step_removed = len(mutations[~mutations.unit])

    # Exclue outlier childs based on CHILD_THRESH standard deviations from mean
    outlier_childs = []
    if CHILD_THRESH > 0:
        mean_muts = child_num_total.child_num_total.mean()
        std_muts = child_num_total.child_num_total.std()
        threshold = mean_muts + (CHILD_THRESH*std_muts)
        outlier_childs = child_num_total[child_num_total.child_num_total > threshold].child.tolist()
        mutations = mutations[~mutations.child.isin(outlier_childs)]

    # Exclude outlier STR
    outlier_loci = []
    if LOCI_THRESH > 0:
        mean_muts = str_num_total.str_num_total.mean()
        std_muts = str_num_total.str_num_total.std()
        threshold = mean_muts + (LOCI_THRESH*std_muts)
        outlier_loci = child_num_total[str_num_total.str_num_total > threshold].str_id.tolist()
        mutations = mutations[~mutations.str_id.isin(outlier_loci)]

    # Log summary statistics after filtering
    child_num_total_filt = mutations.groupby(["family", "child", "phenotype"]).pos.count().reset_index(name='child_num_total')
    str_num_total_filt = mutations.groupby(["chrom", "pos", "period"]).child.count().reset_index(name='str_num_total')

    # Exclude dn STRs based on posterior threshold
    if POSTERIOR_THRESH > 0:
        mutations = mutations[mutations.posterior >= POSTERIOR_THRESH]

    # Write out filtered mutations
    mutations.to_csv(os.path.join(FILTER_FILE), sep="\t", header=True, index=False)

    # Log summary statistics n based on posterior threshold
    child_num_total_pst = mutations.groupby(["family", "child", "phenotype"]).pos.count().reset_index(name='child_num_total')
    str_num_total_pst = mutations.groupby(["chrom", "pos", "period"]).child.count().reset_index(name='str_num_total')

    print("Before filtering:")
    print("Mean, STD # STR genotyped per child {}, {}".format(child_num_total.child_num_total.mean(), child_num_total.child_num_total.std()))
    print ("Mean, STD # child per STR loci {}, {}".format(str_num_total.str_num_total.mean(), str_num_total.str_num_total.std()))
    print("Removed {} outlier children:\n{}".format(len(outlier_childs),  outlier_childs))
    print("Removed {} step size STR loci".format(n_step_removed))
    print("Removed {} outlier STR loci".format(len(outlier_loci)))
    print("After filtering:")
    print("Mean, STD # STR genotyped per child {}, {}".format(child_num_total_filt.child_num_total.mean(), child_num_total_filt.child_num_total.std()))
    print ("Mean, STD # child per STR loci {}, {}".format(str_num_total_filt.str_num_total.mean(), str_num_total_filt.str_num_total.std()))
    print("Mean, STD # STR mutations per child {}, {}".format(child_num_total_pst.child_num_total.mean(), child_num_total_pst.child_num_total.std()))
    print ("Mean, STD # child per STR loci {}, {}".format(str_num_total_pst.str_num_total.mean(), str_num_total_pst.str_num_total.std()))


if __name__ == "__main__":
    main()
