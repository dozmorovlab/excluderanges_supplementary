'''
BedCombine
by Brydon Wall
'''
import argparse
import os
from pybedtools import BedTool
import shutil

def main(args: argparse.Namespace):

    # Args
    in_hs_path = args.in_HS_path
    in_lm_path = args.in_LM_path
    in_cent_path = args.in_cent_path
    out_bed_path = args.out_bed_path
    d = int(args.d)
    #A = args.A

    print(f'in_hs_path: {in_hs_path}')
    print(f'in_lm_path: {in_lm_path}')
    print(f'in_cent_path: {in_cent_path}')
    print(f'out_bed_path: {out_bed_path}')
    print(f'd: {d}')

    os.mkdir('./temp')

    # Combine and merge all
    if not in_cent_path is None:
        lm_comb: BedTool = BedTool(in_cent_path).cat(BedTool(in_lm_path), postmerge = False).sort()
        combined: BedTool = BedTool(in_hs_path).cat(lm_comb, postmerge = False).sort()
    else:    
        combined: BedTool = BedTool(in_hs_path).cat(BedTool(in_lm_path), postmerge = False).sort()
    
    all: BedTool = combined.merge(d = d)
    all.sort()
    all.saveas('./temp/all.bed')

    hs: BedTool = BedTool('./temp/all.bed').intersect(BedTool(in_hs_path), wa = True, u = True)
    lm: BedTool = BedTool('./temp/all.bed').intersect(BedTool(in_lm_path), wa = True, u = True)
    if not in_cent_path is None:
        cm: BedTool = BedTool('./temp/all.bed').intersect(BedTool(in_cent_path), wa = True, u = True)

    # Both HS and LM
    both_hs_lm: BedTool = hs.intersect(lm, wa = True, u = True)

    # HS only
    hs_only: BedTool = hs.subtract(lm, A = True)

    # LM only
    lm_only: BedTool = lm.subtract(hs, A = True)

    if not in_cent_path is None:

        # CM only
        cm_only: BedTool = cm.subtract(hs, A = True).subtract(lm, A = True)
        cm_only = cm_only.each(
            lambda row: (row[0], row[1], row[2], 'Centromere')
        )
        cm_only.saveas('./temp/cm.bed')

        # HS and CM
        both_hs_cm: BedTool = hs.intersect(cm, wa = True, u = True).subtract(lm, A = True)
        both_hs_cm = both_hs_cm.each(
            lambda row: (row[0], row[1], row[2], 'Centromere, High Signal')
        )
        both_hs_cm.saveas('./temp/both_hs_cm.bed')

        # LM and CM
        both_lm_cm: BedTool = lm.intersect(cm, wa = True, u = True).subtract(hs, A = True)
        both_lm_cm = both_lm_cm.each(
            lambda row: (row[0], row[1], row[2], 'Centromere, Low Mappability')
        )
        both_lm_cm.saveas('./temp/both_lm_cm.bed')

        # CM, HS, and LM
        all_cm_hs_lm: BedTool = both_hs_lm.intersect(cm, wa = True, u = True)
        all_cm_hs_lm = all_cm_hs_lm.each(
            lambda row: (row[0], row[1], row[2], 'Centromere, High Signal, Low Mappability')
        )
        all_cm_hs_lm.saveas('./temp/all_cm_hs_lm.bed')

        # HS and LM
        both_hs_lm: BedTool = both_hs_lm.subtract(cm, A = True)

        # HS
        hs_only: BedTool = hs_only.subtract(cm, A = True)

        # LM
        lm_only: BedTool = lm_only.subtract(cm, A = True)

    both_hs_lm = both_hs_lm.each(
        lambda row: (row[0], row[1], row[2], 'High Signal, Low Mappability')
    )
    both_hs_lm.saveas('./temp/both_hs_lm.bed')
    
    hs_only = hs_only.each(
        lambda row: (row[0], row[1], row[2], 'High Signal')
    )
    hs_only.saveas('./temp/hs.bed')
    lm_only = lm_only.each(
        lambda row: (row[0], row[1], row[2], 'Low Mappability')
    )
    lm_only.saveas('./temp/lm.bed')

    # Combine final file
    final: BedTool = BedTool('./temp/both_hs_lm.bed').cat(
        BedTool('./temp/hs.bed'), postmerge = False
    ).cat(
        BedTool('./temp/lm.bed'), postmerge = False
    )

    if not in_cent_path is None:
        final = final.cat(
            BedTool('./temp/cm.bed'), postmerge = False
        ).cat(
            BedTool('./temp/both_hs_cm.bed'), postmerge = False
        ).cat(
            BedTool('./temp/both_lm_cm.bed'), postmerge = False
        ).cat(
            BedTool('./temp/all_cm_hs_lm.bed'), postmerge = False
        )

    final.sort().saveas(out_bed_path)

    shutil.rmtree('./temp')

if __name__ == '__main__':

    # Default Parameters
    in_hs_path = './test_HS.bed'
    in_lm_path = './test_LM.bed'
    out_bed_path = './test_HS_LM.bed'
    d = 0

    parser = argparse.ArgumentParser(description="Combine two beds and annotate (HS or LM)")

    # Add arguments
    parser.add_argument(
        '-s', "--in_HS_path", type=str, default=in_hs_path, help="Input high signal bed path"
    )
    parser.add_argument(
        '-m', "--in_LM_path", type=str, default=in_lm_path, help="Input low mappability bed path"
    )
    parser.add_argument(
        '-c', "--in_cent_path", type=str, default=None, help="Input centromere bed path"
    )
    parser.add_argument(
        '-o', "--out_bed_path", type=str, default=out_bed_path, help="Output bed path"
    )
    #parser.add_argument(
    #    '-A', "--intersect_whole_regions", type=bool, default=True, help="Same as -A in bedtools - Remove entire feature if any overlap.  That is, by default, only subtract the portion of A that overlaps B. Here, if any overlap is found, the entire feature is removed."
    #)
    parser.add_argument(
        "-d", type=int, default=d, help="""Maximum distance between features allowed for features
		to be merged.
		- Def. 0. That is, overlapping & book-ended features are merged.
		- (INTEGER)
		- Note: negative values enforce the number of b.p. required for overlap."""
    )
    # Parse arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args)
