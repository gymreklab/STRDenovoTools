import argparse
import os
import sys
import pandas as pd

def main():
    description = """Apply the provided set of sample and locus  filters to the input de novo STR mutations file.
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--all-mutations-file",           type=str,   required=True,   dest="MUT_FILE",                         help="Input de novo STR mutations file.")
    parser.add_argument("--filtered-mutations-file",      type=str,   required=False,  dest="FILTER_FILE",                      help="Output filtered de novo STR mutations file.")
    parser.add_argument("--log-file",                     type=str,   required=False,  dest="LOG_FILE",                         help="Output log file.")
    parser.add_argument("--filter-denovos-child",         type=int,   required=False,  dest="CHILD_THRESH",         default=0,  help="Use to exclude children with mutations above standard deviations from mean.(Recommend 5)")
    parser.add_argument("--filter-loc-denovos",           type=int,   required=False,  dest="LOCI_THRESH",          default=0,  help="Use to exclude STR loci with called mutateions above standard deviations from mean.(Recommend 5)")
    parser.add_argument("--filter-posterior",             type=float, required=False,  dest="POSTERIOR_THRESH",     default=0,  help="Poserior probability threshould. (Recommend 0.8)")
    parser.add_argument("--filter-step-size",             action="store_true",  required=False,      dest="STEP_SIZE",        default=False,          help="Use to exclude mutations not a unit of the STR motif size.")
    # parser.add_argument("--min-denovos-child",       type=int,   required=False, dest="",                 default=0,        help="")
    # parser.add_argument("--max-denovos-child",       type=int,   required=False, dest="",             default=10000,        help="")
    # parser.add_argument("--min-loc-denovos",         type=int,   required=False, dest="",         default=0,        help=help_dict["--min-loc-denovos"])
    # parser.add_argument("--max-loc-denovos",         type=int,   required=False, dest="",         default=0,        help=help_dict["--max-loc-denovos"])
    args          = parser.parse_args()

    # Load all de novo mutation
    mutations = pd.read_table(os.path.join(args.MUT_FILE), sep="\t")
    if isinstance(mutations["chrom"][0], str):
        mutations["chrom"] = mutations["chrom"].map(lambda x: x.lstrip("chr"))
    mutations["str_id"] = mutations.apply(lambda x: str(x["chrom"])+"_"+str(x["pos"]), axis=1) #apply fn to each row

    # Log summary statistics before filtering
    # Number de novo per person
    child_num_total = mutations.groupby(["family", "child", "phenotype"]).pos.count().reset_index(name='child_num_total')
    # Number de novo per loci
    str_num_total = mutations.groupby(["chrom", "pos", "period", "str_id"]).child.count().reset_index(name='str_num_total')

    # Exclude dn STRs where mutation size is not multiple of period
    n_step_removed=0
    if args.STEP_SIZE:
        mutations["unit"] = mutations["mutsize"]%mutations["period"] == 0
        n_step_removed = len(mutations[~mutations.unit])
        mutations = mutations[mutations.unit]


    # Exclude dn STRs based on posterior threshold
    if args.POSTERIOR_THRESH > 0:
        mutations = mutations[mutations.posterior >= args.POSTERIOR_THRESH]

    # Log summary statistics n based on posterior threshold
    child_num_total_pst = mutations.groupby(["family", "child", "phenotype"]).pos.count().reset_index(name='child_num_total')
    str_num_total_pst = mutations.groupby(["chrom", "pos", "period"]).child.count().reset_index(name='str_num_total')

    # Exclue outlier childs based on CHILD_THRESH standard deviations from mean
    outlier_childs = []
    if args.CHILD_THRESH > 0:
        mean_muts = child_num_total.child_num_total.mean()
        std_muts = child_num_total.child_num_total.std()
        threshold = mean_muts + (args.CHILD_THRESH*std_muts)
        outlier_childs = child_num_total[child_num_total.child_num_total > threshold].child.tolist()
        mutations = mutations[~mutations.child.isin(outlier_childs)]

    # Exclude outlier STR based on standard deviations from mean
    outlier_loci = []
    if args.LOCI_THRESH > 0:
        mean_muts = str_num_total.str_num_total.mean()
        std_muts = str_num_total.str_num_total.std()
        threshold = mean_muts + (args.LOCI_THRESH*std_muts)
        outlier_loci = str_num_total[str_num_total.str_num_total > threshold].str_id.tolist()
        mutations = mutations[~mutations.str_id.isin(outlier_loci)]

    # Log summary statistics after filtering
    child_num_total_filt = mutations.groupby(["family", "child", "phenotype"]).pos.count().reset_index(name='child_num_total')
    str_num_total_filt = mutations.groupby(["chrom", "pos", "period", ]).child.count().reset_index(name='str_num_total')

    # Exclude dn STRs based on posterior threshold
    if args.POSTERIOR_THRESH > 0:
        mutations = mutations[mutations.posterior >= args.POSTERIOR_THRESH]

    # Write out filtered mutations
    mutations.to_csv(os.path.join(args.FILTER_FILE), sep="\t", header=True, index=False)

    with open(os.path.join(args.LOG_FILE), 'w') as log:
        log.write("Before filtering:\n")
        log.write("Min={}, Max={}, Mean={:.2f}, STD={:.2f} STRs genotyped per child.\n".format(child_num_total.child_num_total.min(),
                                                                                     child_num_total.child_num_total.max(),
                                                                                     child_num_total.child_num_total.mean(),
                                                                                     child_num_total.child_num_total.std()))
        log.write ("Mean={:.2f}, STD={:.2f} children per STR loci.\n".format(str_num_total.str_num_total.mean(), str_num_total.str_num_total.std()))
        log.write("\nAfter excluding STRs with <{} posterior threshold:\n".format(args.POSTERIOR_THRESH))
        log.write("Min={}, Max={}, Mean={:.2f}, STD={:.2f} # STR mutations per child\n".format(child_num_total_pst.child_num_total.min(),
                                                                                     child_num_total_pst.child_num_total.max(),
                                                                                     child_num_total_pst.child_num_total.mean(),
                                                                                     child_num_total_pst.child_num_total.std()))
        log.write ("Mean={:.2f}, STD={:.2f} children per STR loci\n".format(str_num_total_pst.str_num_total.mean(), str_num_total_pst.str_num_total.std()))
        log.write("By CHR Mean: {}\n".format([round(str_num_total_pst[str_num_total_pst.chrom == c].str_num_total.mean(),2) for c in range(1,23)]))

        log.write("After excluding  >{} std children threshold and >{} std STRs threhold:\n".format(args.CHILD_THRESH, args.LOCI_THRESH))
        log.write("Min={}, Max={}, Mean={:.2f}, STD={:.2f} STRs genotyped per child.\n".format(child_num_total_filt.child_num_total.min(),
                                                                                     child_num_total_filt.child_num_total.max(),
                                                                                     child_num_total_filt.child_num_total.mean(),
                                                                                     child_num_total_filt.child_num_total.std()))
        log.write ("Mean={:.2f}, STD={:.2f} children per STR loci\n".format(str_num_total_filt.str_num_total.mean(), str_num_total_filt.str_num_total.std()))
        log.write("Removed {} outlier children:\n{}\n".format(len(outlier_childs),  outlier_childs))
        log.write("Removed {} step size STR loci\n".format(n_step_removed))
        log.write("Removed {} outlier STR loci\n".format(len(outlier_loci)))


if __name__ == "__main__":
    main()
