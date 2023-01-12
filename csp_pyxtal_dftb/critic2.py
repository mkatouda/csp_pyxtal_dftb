import sys
import os
import re
import glob
import argparse
import subprocess
import csv
from time import time

import yaml
import numpy as np
import matplotlib.pyplot as plt
from rdkit.ML.Cluster import Butina


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
        '-cc', '--clustercsv', type=str, default='cluster_profile.csv',
        help = 'Clustring profile of crystal structure output csv filename'
    )
    parser.add_argument(
        '-cp', '--clusterpng', type=str, default='cluster_profile.png',
        help = 'Clustring profile of crystal structure output png filename'
    )
    parser.add_argument(
        '--comparison', type=str, default='POWDER',
        help = 'Reference crystal structure input filename: [POWDER|RDF|AMD]'
    )
    parser.add_argument(
        '--cluster-cutoff', type=float, default=0.5,
        help = 'cutoff of distant matrix of crystals for clustering'
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

def plot_cluster_profile(clusters, pngfile):
    fig = plt.figure(1, figsize=(10, 8))
    plt1 = plt.subplot(111)
    plt.axis([0, len(clusters), 0, len(clusters[0])+1])
    plt.xlabel('Cluster index', fontsize=20)
    plt.ylabel('Number of molecules', fontsize=20)
    plt.tick_params(labelsize=16)
    plt1.bar(range(1, len(clusters)), [len(c) for c in clusters[:len(clusters)-1]], lw=0)
    plt.savefig(pngfile)
    plt.clf()
    plt.close()

def compare_crys(infmt='cif', refstrucfile=None, comparison='POWDER',
                 cluster_cutoff=0.5, diffmatfile='diffmat.csv',
                 struclistfile='struclist.csv', clustercsvfile='cluster_profile.csv',
                 clusterpngfile='cluster_profile.png', verbose=False):
                 
    critic2in = 'cryscompare.cri'
    critic2out = 'cryscompare.cro'
    critic2err = 'cryscompare.cre'

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [atoi(c) for c in re.split(r'(\d+)', text)]

    strucfiles = sorted(glob.glob("*." + infmt), key=natural_keys)
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

    t1 = time()
    diff_mat = parse_diffmat(results, ncrys, verbose=verbose)
    np.savetxt(diffmatfile, diff_mat, delimiter=',', fmt='%.7f')

    dist_mat_tri = []
    for i in range(ncrys):
        for j in range(i):
            dist_mat_tri.append(diff_mat[i, j])

    if verbose:
        print('Distance matrix:')
        ij = 0
        for i in range(ncrys):
            for j in range(i):
                print(i, j, dist_mat_tri[ij])
                ij += 1

    t1 = time()
    print('Start Butina crystal clustering using critic2 results:')
    clusters_tmp = Butina.ClusterData(dist_mat_tri, ncrys, cluster_cutoff,
                                      isDistData=True, reordering=True)
    tc = time() - t1
    print('Finish Butina crystal clustering:\n Elapsed time (s):', tc)

    clusters = []
    for i in range(len(clusters_tmp)):
        clusters.append([c+1 for c in clusters_tmp[i]])

    with open(clustercsvfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(clusters)

    print('Number of clusters:', len(clusters))
    if verbose:
        for i in range(len(clusters)):
            print('Cluster_{}:'.format(i), len(clusters[i]), clusters[i])

    plot_cluster_profile(clusters, clusterpngfile)

    return clusters

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
    cluster_cutoff = conf['cluster_cutoff']
    diffmatfile = conf['diffmat']
    struclistfile = conf['struclist']
    clustercsvfile = conf['clustercsv']
    clusterpngfile = conf['clusterpng']
    verbose = conf['verbose']

    compare_crys(infmt=infmt, refstrucfile=refstrucfile, comparison=comparison,
                 cluster_cutoff=cluster_cutoff, diffmatfile=diffmatfile,
                 struclistfile=struclistfile, clustercsvfile=clustercsvfile,
                 clusterpngfile=clusterpngfile, verbose=verbose)

if __name__ == '__main__':
    main()
