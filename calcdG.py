"""

Written by Shyam Saladi (saladi@caltech.edu), January 2016

Rewrite of Nanjiang Shu's calc_dG.pl, created 2010-08-27, updated 2010-09-13
Modified by Shyam Saladi, May 2012

"""
import math
import fileinput

import numpy as np


def main():
    for line in fileinput.input():
        line = line.strip()
        print(line, test_scan_for_best_TM_dGraw(helix=line))
    return


def test_scan_for_best_TM_dGraw():
    # run tests
    return


def scan_for_best_TM_dGraw(helix, profile='biological', allow_sub=False,
                           with_length=True):
    # accomidate different profiles
    if type(profile) is dict:
        profile = dict
    elif type(profile) is str and profile == 'biological':
        profile = biological
    else:
        raise ValueError

    # set length correction coefficients
    if with_length:
        c = (length_correction_coeff[0], 0, 0, 0)
    else:
        c = length_correction_coeff

    if allow_sub and len(helix) > 9:
        len_min = min(9, len(helix))
    else:
        len_min = len(helix)

    lowest = {'startidx': -1, 'length': -1, 'dG': 1000000}

    # lengths to scan over
    for thislen in range(len_min, len(helix) + 1):
        # windows to scan over
        for startidx in range(0, len(helix) - thislen):
            # intermediate calcs
            dg = 0
            dg_sum = 0
            dg_sin_sum = 0
            dg_cos_sum = 0
            # add for residue at each position
            for i in range(thislen):
                dg = pos_spec_dG(helix[startidx+i], i, thislen, profile)
                dg_sum += dg
                dg_sin_sum += dg * np.sin(100 * (i) * piover180)
                dg_cos_sum += dg * np.cos(100 * (i) * piover180)

            segment_dG = dg_sum + c[0] * np.sqrt(dg_sin_sum**2 + dg_cos_sum**2)

            # Correct for length
            segment_dG += c[1] + c[2]*L + c[3]*L*L

            if segment_dG < lowest['dG']:
                lowest['dG'] = segment_dG
                lowest['start'] = startidx
                lowest['length'] = thislen

    return lowest


def pos_spec_dG(aa, i, L, profile):
    pos = 9 * (2 * (i)/(L-1) - 1)  # check if consistent with paper (Shyam)
    if aa == "W" or aa == "Y":
        return profile[aa][0] * np.exp(-1*profile[aa][1]*pos**2) + \
            profile[aa][2] * \
                (np.exp(-1*profile[aa][3]*(pos-profile[aa][4])**2) +
                 np.exp(-1*profile[aa][3]*(pos+profile[aa][4])**2))
    else:
        return profile[aa][0] * np.exp(-1*profile[aa][1]*pos**2)


global length_correction_coeff, biological, piover180

piover180 = math.atan2(1, 1)/45
length_correction_coeff = (0.27045, 9.29274167549645,
                           -0.64513139783394, 0.00822196628688)

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
    'X': (0.6348884, 0.0180273),
    'Y': (0.6275249, 0.0103896, -0.5744404, 0.0947821, 6.9164963),
    'Z': (1.3761092, 0.0099898)
}


if __name__ == "__main__":
    main()

"""
# Unused here, but part of original script

opm = { 'A'  : (   -0.35564862900004,   0.01332130021828),
        'C'  : (   -0.03501895680136,   0.99999999999436),
        'D'  : (    1.59956132209071,   0.00947489754213),
        'E'  : (    1.44456437954037,   0.00991793633818),
        'F'  : (   -0.58372719807522,   0.01010636364230),
        'G'  : (   -0.21135828304947,   0.00345594635635),
        'H'  : (    0.98713634884570,   0.21148435891199),
        'I'  : (   -0.53139640795947,   0.02460014333285),
        'K'  : (    1.65710232354055,   0.01319508610467),
        'L'  : (   -0.40964877807885,   0.01877474391562),
        'M'  : (   -0.31236254433492,   0.01063405577742),
        'N'  : (    0.95392420347656,   0.01345897376904),
        'P'  : (    0.47894080618646,   0.02083850626139),
        'Q'  : (    1.13718143293910,   0.01049378516167),
        'R'  : (    1.45575556617996,   0.01812501297927),
        'S'  : (    0.16601192335102,   0.00036243232009),
        'T'  : (   -0.04699968531037,   0.05979821658957),
        'V'  : (   -0.45345841058687,   0.02575637636981),
        'W'  : (   -0.97553693183770,   0.00549871324932,   0.39792640886763,   0.05557263615780,    0.00000000350551),
        'X'  : (      0.34521774973501,   0.01890642845475),
        'Y'  : (   -0.47290915514074,   0.01054870058203,    0.51698866596097,   0.06638567043267,   -2.21478619615780)
        }

kdo = { 'A' : ( -1.8, 0 ),
    'B' : (  3.5, 0 ),
    'C' : ( -2.5, 0 ),
    'D' : (  3.5, 0 ),
    'E' : (  3.5, 0 ),
    'F' : ( -2.8, 0 ),
    'G' : (  0.4, 0 ),
    'H' : (  3.2, 0 ),
    'I' : ( -4.5, 0 ),
    'K' : (  3.9, 0 ),
    'L' : ( -3.8, 0 ),
    'M' : ( -1.9, 0 ),
    'N' : (  3.5, 0 ),
    'P' : (  1.6, 0 ),
    'Q' : (  3.5, 0 ),
    'R' : (  4.5, 0 ),
    'S' : (  0.8, 0 ),
    'T' : (  0.7, 0 ),
    'U' : ( -1.9, 0 ),
    'V' : ( -4.2, 0 ),
    'W' : (  0.9, 0, 0, 0, 0 ),
    'X' : (  0.5, 0 ),
    'Y' : (  1.3, 0, 0, 0, 0 ),
    'Z' : (  3.5, 0 )
    }

zhl = {'A' : (0.38 , 0 ),
    'B' : (2.45 , 0 ),
    'C' : (0.30 , 0 ),
    'D' : (3.27 , 0 ),
    'E' : (2.90 , 0 ),
    'F' : (1.98 , 0 ),
    'G' : (0.19 , 0 ),
    'H' : (1.44 , 0 ),
    'I' : (1.97 , 0 ),
    'K' : (3.46 , 0 ),
    'L' : (1.82 , 0 ),
    'M' : (1.40 , 0 ),
    'N' : (1.62 , 0 ),
    'P' : (1.44 , 0 ),
    'Q' : (1.84 , 0 ),
    'R' : (2.57 , 0 ),
    'S' : (0.53 , 0 ),
    'T' : (0.32 , 0 ),
    'U' : (1.40 , 0 ),
    'V' : (1.46 , 0 ),
    'W' : (1.53 , 0, 0, 0, 0 ),
    'X' : (1.55 , 0 ),
    'Y' : (0.49 , 0, 0, 0, 0 ),
    'Z' : (2,37 , 0 )
    }

aa2nr  =  { "A" : 1,"C" : 2,"D" : 3,"E" : 4,"F" : 5,"G" : 6,"H" : 7,"I" : 8,"K" : 9,"L" : 10,
          "M" : 11,"N" : 12,"P" : 13,"Q" : 14,"R" : 15,"S" : 16,"T" : 17,"V" : 18,"W" : 19, "Y" : 20}
"""
