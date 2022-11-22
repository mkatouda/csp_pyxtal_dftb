import sys
import os
import glob
import argparse
import subprocess
from time import time

import yaml
import numpy as np


critic2_cmd = 'critic2'
if os.getenv('CRITIC2_COMMAND'):
    critic2_cmd = os.getenv('CRITIC2_COMMAND')


def get_parser():
    class customHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawTextHelpFormatter):
        pass

    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
        description='Crystal structure comarison tool using critic2'
    )
    parser = argparse.ArgumentParser(
        formatter_class=customHelpFormatter,
    )
    parser.add_argument(
        '-i', '--inp', type=str,
        help = 'yaml style input file, overwriting argument values',
    )
    parser.add_argument(
        '-ie', '--iext', type=str, default='cif',
        help = 'File extention of crystal structure input files: [cif|res|gen|vasp|xsf|cube]'
    )
    parser.add_argument(
        '-r', '--refstruc', type=str, default=None,
        help = 'Reference crystal structure input filename'
    )
    parser.add_argument(
        '-d', '--diffmat', type=str, default='diffmat.csv',
        help = 'Crystal structure differnce output filename'
    )
    parser.add_argument(
        '-s', '--struclist', type=str, default='struclist.csv',
        help = 'Crystal structure list output filename'
    )
    parser.add_argument(
        '--comparison', type=str, default='POWDER',
        help = 'Reference crystal structure input filename: [POWDER|RDF|AMD]'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help = 'Verbose output'
    )

    return parser.parse_args()

def set_config(args):
    # Read config yaml file
    if args.inp is not None and os.path.isfile(args.inp):
        with open(args.inp, 'r') as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}

    # Set up default config values from program arguments
    conf_def = vars(args).copy()
    [conf.setdefault(k, v) for k, v in conf_def.items()]

    return conf

def run_cmd(cmd, verbose=False):
    print(' '.join(cmd))
    results = subprocess.run(cmd, capture_output=True, check=True, text=True)
    if verbose:
        print('stdout:\n', results.stdout, '\nstderr:\n', results.stderr)
    return results

def parse_diffmat(results, ncrys, verbose=False):

    if ncrys <= 1: return None

    outlist = results.stdout.split('\n')

    diff_find = False
    diff_mat = np.empty((ncrys, ncrys))

    if ncrys >= 3:
        #print('Parser:')
        jbgn = 0
        for line in outlist:
            if 'DIFF' in line and 'DIFF' in line.split()[0]:
                diff_find = True
                i = 0
                jend = jbgn + 5
                if jend > ncrys:
                    jend = ncrys
                    #print('jbgn:', jbgn, 'jend:', jend)
            elif diff_find:
                #print(i, line.split())
                a = line.split()[1:]
                for j in range(jbgn, jend):
                    diff_mat[i, j] = float(a[j-jbgn])
                    #print('i:', i, 'j:', j, a[j-jbgn])
                i += 1
                if i == ncrys:
                    diff_find = False
                    jbgn += 5

    else:
        diff_mat[:,:] = 0.0
        for line in outlist:
            if 'DIFF' in line and 'DIFF' in line.split()[1]:
                diff_mat[0, 1] = diff_mat[1, 0] = float(line.split()[3])
                break

    if verbose: print('diff_mat:', diff_mat.shape, '\n', diff_mat)

    return diff_mat

def compare_crys(infmt='cif', refstrucfile=None, comparison='POWDER', 
                 diffmatfile='diffmat.csv', struclistfile='struclist.csv',
                 verbose=False):
    critic2in = 'cryscompare.cri'
    critic2out = 'cryscompare.cro'
    critic2err = 'cryscompare.cre'

    strucfiles = glob.glob("*." + infmt)
    ncrys = len(strucfiles)
    if ncrys <= 1: return

    s1 = "compare \ \n"
    s2 = ''
    for strucfile in strucfiles:
        s1 += "{} \ \n".format(strucfile)
        s2 += "{},\n".format(strucfile)
    if refstrucfile is not None:
        ncrys += 1
        s1 += "{} \ \n".format(refstrucfile)
        s2 += "{},\n".format(refstrucfile)
    s1 += '{}\n'.format(comparison.upper())
    s1 += 'END\n'
    with open(critic2in, mode='w') as f:
        f.write(s1)
    with open(struclistfile, mode='w') as f:
        f.write(s2)

    print('Numbers of structures in critic2 crystal comparison:', ncrys)

    print('Start critic2 crystal comparison:')
    t1 = time()
    cmd = [critic2_cmd, critic2in]
    results = run_cmd(cmd, verbose=verbose)
    tc = time() - t1
    print('Finish critic2 crystal comparison:\n Elapsed time (s):', tc)

    with open(critic2out, mode='w') as f:
        f.write(results.stdout)
    with open(critic2err, mode='w') as f:
        f.write(results.stderr)

    diff_mat = parse_diffmat(results, ncrys, verbose=verbose)
    np.savetxt(diffmatfile, diff_mat, delimiter=',', fmt='%.7f')

    return diff_mat

def main():
    args = get_parser()
    print(args)

    conf = set_config(args)

    print('======= Input configulations =======')
    for k, v in conf.items():
        print('{}: {}'.format(k, v))
    print('====================================')

    infmt = conf['iext']
    refstrucfile = None
    if conf['refstruc'] is not None:
        refstrucfile = conf['refstruc']
        print('Compare also with reference crystal struture:', refstrucfile)
        if os.path.isfile(refstrucfile):
            refstrucfile = os.path.abspath(refstrucfile)
    comparison = conf['comparison']
    diffmatfile = conf['diffmat']
    struclistfile = conf['struclist']
    verbose = conf['verbose']

    compare_crys(infmt=infmt, refstrucfile=refstrucfile, comparison=comparison,
                 diffmatfile=diffmatfile, struclistfile=struclistfile,
                 verbose=verbose)

if __name__ == '__main__':
    main()
