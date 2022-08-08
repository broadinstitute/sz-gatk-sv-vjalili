#!/usr/bin/env python
# -*- coding: utf-8 -*-
# freq

"""
"""

import argparse
from collections import deque, defaultdict
import numpy as np
import pandas as pd
import pysam
import pybedtools as pbt
import svtk.utils as svu


def process_rdtest(rdtest):
    """Standardize rdtest column names"""

    # Drop metadata columns (available from VCF) and rename CNVID
    skip_cols = 'chr Start End SampleIDs Type'.split()
    rdtest = rdtest.drop(skip_cols, axis=1).rename(columns={'CNVID': 'name'})

    numeric_cols = 'Median_Power P 2ndMaxP Median_Rank Median_Separation'
    numeric_cols = numeric_cols.split()

    # Replace strings with NA
    for col in numeric_cols:
        repl = ['All_samples_called_CNV_no_analysis',
                'No_samples_for_analysis',
                'coverage_failure']
        rdtest[col] = rdtest[col].replace(repl, np.nan).astype(np.float)

    rdtest['log_pval'] = -np.log10(rdtest.P)
    rdtest['log_2ndMaxP'] = -np.log10(rdtest['2ndMaxP'])

    maxp = rdtest.loc[rdtest.log_pval != np.inf, 'log_pval'].max()
    max2p = rdtest.loc[rdtest.log_2ndMaxP != np.inf, 'log_2ndMaxP'].max()

    rdtest.loc[rdtest.log_pval == np.inf, 'log_pval'] = maxp + 5
    rdtest.loc[rdtest.log_2ndMaxP == np.inf, 'log_2ndMaxP'] = max2p + 5

    rdtest.log_pval = rdtest.log_pval.abs()
    rdtest.log_2ndMaxP = rdtest.log_2ndMaxP.abs()

    return rdtest


def process_srtest(srtest):
    metrics = 'SRQ SRCS pos'.split()
    srtest = srtest.pivot(index='name', values=metrics, columns='coord')
    srtest.columns = ['_'.join(col[::-1]).strip()
                      for col in srtest.columns.values]
    srtest = srtest.reset_index()
    return srtest


def process_petest(petest):
    return petest


def process_baftest(baftest):
    baftest['BAFDEL'] = -baftest.BAFDEL
    return baftest


def preprocess(df, dtype):
    if dtype == 'RD':
        return process_rdtest(df)
    elif dtype == 'SR':
        return process_srtest(df)
    elif dtype == 'PE':
        return process_petest(df)
    elif dtype == 'BAF':
        return process_baftest(df)
    else:
        return df


def _is_parent(s):
    return s.endswith('fa') or s.endswith('mo')


def _is_child(s):
    return s.endswith('p1') or s.endswith('s1') or s.endswith('pb')


def fam_info_readin(fam_file):
    fin = open(fam_file)
    # samp_pedi_hash = {}
    [fam, samp, fa, mo] = [[], [], [], []]
    for line in fin:
        pin = line.strip().split()
        fam.append(pin[0])
        samp.append(pin[1])
        fa.append(pin[2])
        mo.append(pin[3])
    fin.close()
    return [fam, samp, fa, mo]


def process_metadata(variants, rd_df):
    n_samples = len(variants.header.samples)
    called_counts = dict()
    called_samples = dict()

    # Calculate segdup coverage
    segdups = pbt.BedTool(args.segdups)

    # Check if endpoints are in repeat-masked sequence
    starts = metadata['chrom start end name'.split()].copy()
    starts['end'] = starts['start'] + 1
    ends = metadata['chrom start end name'.split()].copy()
    ends['start'] = ends['end'] - 1
    endpoints = pd.concat([starts, ends])
    rmsk = pbt.BedTool(args.rmsk)
    sect = bt.intersect(rmsk, u=True)
    rmasked_names = [i.fields[3] for i in sect.intervals]
    metadata['rmsk'] = metadata.name.isin(rmasked_names)

    for svtype in 'DEL DUP INV BND INS'.split():
        # Counts of variants per sample
        called_counts[svtype] = defaultdict(int)

        # List of variants specific to each sample
        called_samples[svtype] = defaultdict(list)

    metadata = deque()
    for variant in variants:
        chrom = variant.chrom
        start = variant.pos
        end = variant.stop
        called = svu.get_called_samples(variant)
        name = variant.id
        svtype = variant.info['SVTYPE']
        if svtype == 'BND':
            svlen = -1
        elif svtype == 'INS':
            svlen = variant.info.get('SVLEN', -1)
        else:
            svlen = end - start

        # Only use start/end for seg dup coverage. if it's a tloc,
        # we don't care so we can just set its "END" to pos + 1
        if end <= start:
            end = start + 1

        # Calculate VF
        vf = len(called) / n_samples

        # Increment counts of variants per sample
        for s in called:
            called_counts[svtype][s] += 1

        # Track called samples for outlier filtering
        called_samples[svtype][name] = set(called)

        cov = bt.coverage(segdups).to_dataframe()
        metadata['poor_region_cov'] = cov.thickStart

    metadata = np.array(metadata)
    cols = 'chrom start end name svtype svsize vf'.split()
    metadata = pd.DataFrame(metadata, columns=cols)

    # Flag variants specific to outlier samples
    metadata['is_outlier_specific'] = False
    for svtype in 'DEL DUP INV BND INS'.split():
        counts = pd.DataFrame.from_dict(called_counts[svtype], orient='index')\
                             .reset_index()\
                             .rename(columns={'index': 'sample', 0: 'var_count'})
        if counts.shape[0] == 0:
            continue

        q1 = counts.var_count.quantile(0.25)
        q3 = counts.var_count.quantile(0.75)
        thresh = q3 + 1.5 * (q3 - q1)
        outliers = counts.loc[counts.var_count >= thresh, 'sample'].values

        flagged = []
        for var_name, samples in called_samples[svtype].items():
            if samples.issubset(outliers):
                flagged.append(var_name)
        metadata.loc[metadata.name.isin(flagged), 'is_outlier_specific'] = True

    for col in 'start end svsize'.split():
        metadata[col] = metadata[col].astype(int)

    return metadata


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--variants', required=True, help='Default VCF')
    parser.add_argument('-r', '--RDtest')
    parser.add_argument('--segdups', required=True)
    parser.add_argument('--rmsk', required=True)
    parser.add_argument('fout')
    args = parser.parse_args()

    variants = pysam.VariantFile(args.variants)

    # Parse RDTest table
    rd_path = getattr(args, 'RDtest')
    if rd_path is not None:
        rd_df = pd.read_table(rd_path)
        rd_df = process_rdtest(rd_df)
        rd_df = rd_df.rename(columns=lambda c: 'RD_' + c if c != 'name' else c)
        rd_df = rd_df.set_index('name')
        #evidence = metadata.join(rd_df, how='outer', sort=True)
    else:
        rd_df = None

    metadata = process_metadata(variants, rd_df)

    evidence = deque()

    evidence = evidence.reset_index().rename(columns={'index': 'name'})

    # Replace infinite log-pvals
    LOG_CEIL = 300
    evidence = evidence.replace(np.inf, LOG_CEIL)

    evidence.to_csv(args.fout, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
