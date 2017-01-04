# -*- coding: utf-8 -*-
"""

Written by Shyam Saladi (saladi@caltech.edu), January 2016

Rewrite of Nanjiang Shu's calc_dG.pl, created 2010-08-27, updated 2010-09-13
Modified by Shyam Saladi, May 2012

"""

import sys
import warnings
import argparse

import numpy as np
import Bio.SeqIO

biological = {
    'A': (0.1267255, 0.0215152),
    'B': (1.5556351, 0.0133523),
    'C': (-0.0765051, 0.0994228),
    'D': (1.7939795, 0.0172643),
    'E': (1.4193720, 0.0089351),
    'F': (-0.2766953, 0.0010297),
    'G': (0.4813492, 0.0047210),
    'H': (1.1998590, 0.0080127),
    'I': (-0.4597384, 0.0181495),
    'K': (1.8485768, 0.0218446),
    'L': (-0.4282992, 0.0023804),
    'M': (-0.0774786, 0.0984413),
    'N': (1.3266132, 0.0092375),
    'P': (1.0860888, 0.0100568),
    'Q': (1.3336109, 0.0111996),
    'R': (1.6492534, 0.0512044),
    'S': (0.7023921, 0.0077661),
    'T': (0.5266550, 0.0311973),
    'U': (-0.0774786, 0.0984413),
    'V': (-0.2447218, 0.0979201),
    'W': (0.2909390, 0.0189282, -0.5479140, 0.0930222, 6.4736619),
    'X': (0.6348884, 0.0180273), # Undetermined/unknown characters
    'Y': (0.6275249, 0.0103896, -0.5744404, 0.0947821, 6.9164963),
    'Z': (1.3761092, 0.0099898)
}
"""dict of tuples: Biological deltaG of insertion scales by residue
"""

length_correction_coeff = (0.27045, 9.29274167549645,
                           -0.64513139783394, 0.00822196628688)

def scan_for_best_TM_dGraw(helix, profile=biological,
                           allow_sub=False,
                           with_length=False,
                           len_corr_coeff=length_correction_coeff):
    """Calculates dG values

    Parameters
    ----------
    helix : str

    profile : Optional[dict]
        Residue profile for deltaG calculation

    allow_sub : Optional[bool]
        Allow subsequences when searching for helix with lowest deltaG of
        insertion

    with_length : Optional[bool]
        Apply length correction to deltaG calculations

    length_correction_coeff : Optional[list(int)]
        Values for length correction

    Returns
    -------
    tuple
        float: deltaG of insertion score,
        int: Position of the start of the helix,
        int: Length of identified helix

    Raises
    ------
    None
    """
    # accommodate different profiles
    if not isinstance(profile, dict):
        raise ValueError("Profile type is not recognized")

    # set length correction coefficients
    if with_length:
        cur_corr_coeff = length_correction_coeff
    else:
        cur_corr_coeff = (length_correction_coeff[0], 0, 0, 0)

    if allow_sub and len(helix) > 9:
        len_min = min(9, len(helix))
    else:
        len_min = len(helix)

    lowest_start = -1
    lowest_length = -1
    lowest_dG = 1000000

    # lengths to scan over
    for thislen in range(len_min, len(helix) + 1):
        # windows to scan over
        for startidx in range(0, len(helix) - thislen + 1):
            # current helix:
            segdG = segment_dG(helix[startidx:startidx+thislen],
                               profile=profile,
                               len_corr=cur_corr_coeff)
            if segdG < lowest_dG:
                lowest_dG = segdG
                lowest_start = startidx
                lowest_length = thislen

    return (lowest_dG, lowest_start, lowest_length)


def segment_dG(helix, profile, len_corr):
    """Calculate deltaG of insertion for a single defined helix

    Parameters
    ----------
    helix : str
        Residues that make up the helix

    profile : dict of tuples
        Energetics profile to use for calculation. Passed to `pos_spec_dG`.

    len_corr : list of floats
        Length correction coefficients.

    Returns
    -------
    float
        delta G of insertion

    Raises
    ------
    None
    """
    dg_sum = 0
    dg_sin_sum = 0
    dg_cos_sum = 0

    for i, res in enumerate(helix):
        dg = pos_spec_dG(res, i+1, len(helix), profile)
        dg_sum += dg
        dg_sin_sum += dg * np.sin(100 * i * np.pi / 180)
        dg_cos_sum += dg * np.cos(100 * i * np.pi / 180)

    return dg_sum + len_corr[0] * np.sqrt(dg_sin_sum**2 + dg_cos_sum**2) + \
        len_corr[1] + len_corr[2]*len(helix) + len_corr[3]*len(helix)**2


def pos_spec_dG(aa, i, L, profile):
    """Calculate the delta G contribution for a single residue within a helix

    Parameters
    ----------
    aa : str
        Residue identity

    i : int
        position in helix

    L : int
        length of helix

    profile : dict
        Energetics profile to use for calculation

    Returns
    -------
    float
        delta G contribution

    Raises
    ------
    None
    """
    if aa not in profile:
        warnings.warn('%s not in profile: `X` substituted')
        aa = 'X'
    try:
        pos = 9 * (2 * (i-1)/(L-1) - 1)
    except ZeroDivisionError:
        warnings.warn("Calculation invalid: len(helix) = 1. Returning NaN")
        return np.nan
    if len(profile[aa]) == 2:
        return profile[aa][0] * np.exp(-1*profile[aa][1]*pos**2)
    elif len(profile[aa]) == 5:
        return profile[aa][0] * np.exp(-1*profile[aa][1]*pos**2) + \
            profile[aa][2] * \
            (np.exp(-1*profile[aa][3]*(pos-profile[aa][4])**2) +
             np.exp(-1*profile[aa][3]*(pos+profile[aa][4])**2))
    else:
        raise ValueError("Profile type for this residue is unrecognized")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate deltaG of TM insertion.')

    parser.add_argument('fna_filename',
        metavar='fna_file',
        type=str,
        default='-',
        help='FASTA-formatted nucleotide file of the protein sequences'
             '(stdin by default)')

    parser.add_argument('--output',
        metavar='output_file',
        type=str,
        default=sys.stdout,
        help='Specifies name of file to write to (stdout by default)')

    args = parser.parse_args()

    if args.fna_filename == '-':
        print("Reading sequences from stdin", file=sys.stderr)
        seq_records = Bio.SeqIO.parse(args.fna_file, "fasta")
    else:
        with open(args.fna_filename, 'r+') as fh:
            seq_records = list(Bio.SeqIO.parse(fh, "fasta"))

    for record in seq_records:
        print(scan_for_best_TM_dGraw(str(record.seq)))


if __name__ == '__main__':
    main()
