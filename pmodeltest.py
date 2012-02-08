#!/usr/bin/python
"""
14 May 2010

Tests evolutionary models like jModelTest but also for proteins.
'bitfast' option allows to run topology optimization on most weighted models.
Return 'best' tree with SH branch support value.

All based on the behavior of jModelTest:
 jModelTest: phylogenetic model averaging.
 Posada D
 Mol Biol Evol25p1253-6(2008 Jul)

and ProtTest:
 ProtTest: selection of best-fit models of protein evolution.
 Abascal F, Zardoya R, Posada D
 Bioinformatics21p2104-5(2005 May 1)

WARNING: only working with phyml version: v3.0_360-500M
         available at: http://www.atgc-montpellier.fr/phyml/
"""

__author__  = "francois serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "1.04"
__title__   = "pmodeltest v%s" % __version__

from argparse import ArgumentParser
from subprocess import Popen, PIPE
from re import sub
from sys import stderr as STDERR
from sys import stdout as STDOUT
from cmath import exp
from time import sleep

def get_job_list(algt, wanted_models, speed=True, verbose=False, protein=False,
              support=False, sequential=False, rerun=False):
    '''
    runs a list of models and returns a dictionary of results.
    If rerun is True, will have into account +I+G information
    and only run the wanted ones.
    '''
    opt = 'l'  if speed   else 'tl'
    sup = '-4' if support else '0'
    typ = 'aa' if protein else 'nt'
    job_list = {}
    for model in models [typ]:
        for freq in freqs[typ].keys():
            if modelnames[typ][model+freq][0] not in wanted_models:
                continue
            model_name  = modelnames[typ][model + freq][0]
            model_param = modelnames[typ][model + freq][1]
            for inv in invts.keys():
                for gam in gamma.keys():
                    job = model_name + inv + gam + freq
                    job_list [job] = {
                        'cmd': ['phyml', '--sequential'*sequential,
                                '-i', algt,'-d',
                                'aa' if protein else 'nt',
                                '-n', '1',
                                '-b', sup,
                                '-o', opt,
                                '-m', model,
                                '--run_id', job] + freqs[typ][freq] + invts[inv] + gamma[gam],
                        'params': model_param + int (opt=='tl') + (inv != '') + (gam != ''),
                        'algt': algt
                    }
    if verbose:
        for job in job_list:
            print job + ': ' + ' '.join(job_list[job]['cmd'])
    return job_list

def run_jobs(job_list, nprocs=1, refresh=2):
    procs = {}
    done = 0
    todo = len (job_list)
    jobs = job_list.keys()[:]
    try:
        while True:
            if len(procs)<nprocs and jobs:
                job = jobs.pop()
                procs[job] = {'p': Popen(job_list[job]['cmd'], stderr=PIPE, stdout=PIPE),
                              'job': job}
                sleep(refresh)
                continue
            for p in procs:
                if procs[p]['p'].poll() is None:
                    continue
                if procs[p]['p'].returncode == -9:
                    print ' WAHOOO!!! this was killed:'
                    print procs[p]
                    return
                out, err = procs[p]['p'].communicate()
                if 'Err: ' in out:
                    exit ('ERROR: problem running phyml: '+out)
                del procs[p]
                done += 1
                job_list[job]['out'] = out
                job_list[job]['err'] = err
                break
            sleep(refresh)
            if done == todo:
                break
    except Exception as e:
        print 'ERROR at', job
        print e
    return job_list

def parse_jobs(job_list, algt):
    for job in job_list:
        (numspe, lnl, dic) = parse_stats(job_list[job]['algt'] + '_phyml_stats_%s.txt' % job)
        numparam = job_list[job]['params'] + numspe*2-3
        aic = 2*numparam-2*lnl
        job_list[job]['AIC' ] = aic
        job_list[job]['lnL' ] = lnl
        job_list[job]['K'   ] = numparam
        job_list[job]['dic' ] = dic
        job_list[job]['tree'] = get_tree (algt + '_phyml_tree_%s.txt' % job)
    return job_list

def aic_calc(results, speed):
    '''
    compute and displays AICs etc... 
    '''
    ord_aic = sorted ([ [results[x]['AIC'], x] for x in results.keys() ])
    ord_aic = [ x[1] for x in ord_aic ]
    min_aic = results[ord_aic[0]]['AIC']
    for model in ord_aic:
        results[model]['deltar'] =  results[model]['AIC'] - min_aic
        results[model]['weight'] = exp (-0.5 * results[model]['deltar']).real
    sumweight = sum ([ results[x]['weight'] for x in results.keys() ])
    cumweight = 0
    good_models = []
    for model in ord_aic:
        results[model]['weight'] = results[model]['weight']/sumweight
        cumweight += results[model]['weight']
        results[model]['cumweight'] = cumweight
        if results[model]['cumweight'] < 0.9999:
            good_models.append([model, int(1000*results[model]['weight']+0.5)])
    print >> STDOUT,  '\n\n*************************************'
    print >> STDOUT,      '         TABLE OF WEIGHTS'
    if speed:
        print >> STDOUT,  ' (WARNING: no topology optimization)'
    print >> STDOUT,      '*************************************\n'
    print >> STDOUT,  '   %-14s | %-4s | %-9s | %-8s | %-9s | %-6s | %-5s ' % \
          ('MODEL', 'K', 'lnL', 'AIC', 'delta AIC',
           'Weight', 'Cumulative weights')
    print >> STDOUT, '   ' + '-'*80
    print >> STDOUT, \
          '\n'.join ([ '   %-15s|'%(x) + \
                       ' %-4s |'   % str   (results[x]['K'])         +\
                       ' %-9.2f |' % float (results[x]['lnL'])       +\
                       ' %-8.2f |' % float (results[x]['AIC'])       +\
                       ' %-9.3f |' % float (results[x]['deltar'])    +\
                       ' %-6.3f |' % float (results[x]['weight'])    +\
                       ' %-5.3f'   % float (results[x]['cumweight']) \
                       for x in ord_aic ])
    print >> STDOUT, '\n'
    return results, ord_aic

def re_run(job_list, algt, cutoff=0.95, nprocs=1,refresh=2):
    for job in job_list.keys()[:]:
        if job_list[job]['cumweight'] > 0.95:
            del job_list[job]
    print >> STDOUT,  '\nREFINING...\n    doing the same but computing topologies' + \
                      ' only for models that sums a weight of 0.95\n\n    ' + \
                      '\n    '.join(job_list.keys()) + '\n'
    job_list = run_jobs(job_list, nprocs=nprocs,refresh=refresh)
    return parse_jobs(job_list, algt)
    
def main():
    '''
    main function when called by command line.
    infile must be in phylip format.
    '''
    opts = get_options()
    # remove gamma inv and frequencies if not wanted
    if opts.nogam:
        del gamma ['+G']
    if opts.noinv:
        del invts ['+I']
    if opts.nofrq:
        del freqs ['nt']['+F']
        del freqs ['aa']['+F']

    # first run of models
    results = get_job_list(opts.algt, opts.models, opts.speedy, opts.verb,
                           opts.protein, opts.support, sequential=opts.sequential)
    results = run_jobs(results,nprocs=opts.nprocs,refresh=opts.refresh)
    results = parse_jobs(results, opts.algt)
    results, ord_aic = aic_calc(results, opts.speedy)

    # if bit fast, second run with wanted models (that sums weight of 0.95)
    if opts.medium:
        results = re_run(results, opts.algt, cutoff=0.95, refresh=opts.refresh, nprocs=opts.nprocs)
        results, ord_aic = aic_calc(results, False)

    print >> STDOUT,  '\n\n*************************************************'
    results[ord_aic[0]]['cmd'][results[ord_aic[0]]['cmd'].index ('-o') + 1] += 'r'
    print >> STDOUT,\
          'Re-run of best model with computation of rates and support...'
    cmd = results[ord_aic[0]]['cmd']
    # add best tree search and support to phyml command line
    cmd [cmd.index ('-b')+1] = '-4'
    cmd += ['-s', 'BEST']
    print >> STDOUT, '  Command line = ' + ' '.join (cmd) + '\n'
    # run last phyml
    (out, err) = Popen('echo "end" | ' + ' '.join (cmd),
                       stdout=PIPE, shell=True).communicate()
    if err is not None or 'Err: ' in out:
        exit ('ERROR: problem at last run of phyml: '+out)
    tree = get_tree   (opts.algt + '_phyml_tree_%s.txt' % ord_aic[0])
    print >> STDOUT, '\n Corresponding estimations of rates/frequencies:\n'
    print_model_estimations (parse_stats (opts.algt + '_phyml_stats_%s.txt' % ord_aic[0])[2])
    print >> STDOUT, '\nTree corresponding to best model, '\
          + ord_aic[0] + ' (with SH-like branch supports alone)\n'
    print >> STDOUT,  tree
    if opts.outfile:
        open (opts.outfile, 'w').write (tree)
    if opts.outtrees:
        out_t = open (opts.outtrees, 'w')
        for run in results:
            out_t.write ('command: ' + \
                         ' '.join (results [run]['cmd']) + \
                         '\ntree (nw):    ' + results [run]['tree'] + '\n')
        out_t.close ()

    print >> STDOUT, "Done."

def print_model_estimations (dic):
    '''
    prints table with estimation of rates/frequencies done by phyML
    '''
    for key in filter (lambda x: not x.startswith ('freq') and \
                       not x.startswith ('rate'), dic):
        print >> STDOUT, '     %-20s = %s' % (key, dic [key])
    print >> STDOUT, ''
    for key in filter (lambda x: x.startswith ('freq'), dic):
        print >> STDOUT, '     %-20s = %s' % (key, dic [key])
    print >> STDOUT, ''
    for key in filter (lambda x: x.startswith ('rate'), dic):
        print >> STDOUT, '     %-20s = %s' % (key, dic [key])
    print >> STDOUT, ''

def parse_stats(path):
    '''
    parse stats file of phyml, to extract the likelyhood value
    '''
    dic = {}
    for line in open(path, "rU"):
        if line.startswith('. Log-likelihood:'):
            lnl          = float (line.strip().split()[-1])
        elif line.startswith('. Number of taxa:'):
            numspe       = int (line.strip().split()[-1])
        elif line.startswith('  - f(A)= '):
            dic['frequency (A)']    = float (line.strip().split()[-1])
        elif line.startswith('  - f(T)= '):
            dic['frequency (T)']    = float (line.strip().split()[-1])
        elif line.startswith('  - f(G)= '):
            dic['frequency (G)']    = float (line.strip().split()[-1])
        elif line.startswith('  - f(C)= '):
            dic['frequency (C)']    = float (line.strip().split()[-1])
        elif line.startswith('  A <-> C '):
            dic['rate A <-> C']    = float (line.strip().split()[-1])
        elif line.startswith('  A <-> G'):
            dic['rate A <-> G']    = float (line.strip().split()[-1])
        elif line.startswith('  A <-> T'):
            dic['rate A <-> T']    = float (line.strip().split()[-1])
        elif line.startswith('  C <-> G'):
            dic['rate C <-> G']    = float (line.strip().split()[-1])
        elif line.startswith('  C <-> T'):
            dic['rate C <-> T']    = float (line.strip().split()[-1])
        elif line.startswith('  C <-> G'):
            dic['rate C <-> G']    = float (line.strip().split()[-1])
        elif line.startswith('. Proportion of invariant:'):
            dic['prop. of invariant']   = float (sub ('.*([0-9]+\.[0-9]+).*',
                                       '\\1', line.strip()))
        elif line.startswith('  - Gamma shape parameter:'):
            dic['gamma shape'] = float (sub ('.*([0-9]+\.[0-9]+).*',
                                       '\\1', line.strip()))
    return (numspe, lnl, dic)

def get_tree(path):
    '''
    just returns the tree generated by phyml
    '''
    for line in open(path):
        if line.startswith('('):
            return line

def get_options():
    '''
    parse option from command line call
    '''
    model_list = {'dna': ['JC', 'K80', 'TrNef', 'TPM1', 'TPM2', 'TPM3',
                          'TIM1ef', 'TIM2ef', 'TIM3ef', 'TVMef', 'SYM',
                          'F81', 'HKY', 'TrN', 'TPM1uf', 'TPM2uf', 'TPM3uf',
                          'TIM1', 'TIM2', 'TIM3', 'TVM', 'GTR'],
                  'aa' :  ['LG', 'WAG', 'JTT', 'MtREV', 'Dayhoff', 'DCMut',
                           'RtREV', 'CpREV', 'VT', 'Blosum62', 'MtMam',
                           'MtArt', 'HIVw', 'HIVb'] }
    class Choose_model(list):
        def __init__(self,vals):
            self.vals = reduce(lambda x,y:x+y, vals)
        def __contains__(self,x):
            return all ([True if i in self.vals else False for i in x.split(',')])
        def __iter__(self):
            return self.vals.__iter__()
    choose_model = Choose_model(model_list.values())
    parser = ArgumentParser(
        version=__title__,
        description="""Reads sequences from file fasta format, and align according to translation.
        DNA models availabe are: %s.
        AA models available are: %s""" % (','.join(model_list['dna']), ','.join(model_list['aa'])))
    parser.add_argument('-i', dest='algt', type=str, required=True,
                        help='path to input file in phylip format')
    parser.add_argument('-o', dest='outfile', type=str, 
                        help='name of outfile tree (newick format)')
    parser.add_argument('-O', dest='outtrees', metavar="PATH",
                        help='name of outfile with all trees (newick format)')
    parser.add_argument('--support', action='store_true',
                        dest='support', default=False,
                        help='[%(default)s] compute SH-like branch support for each model (slower).')
    parser.add_argument('--fast', action='store_true',
                        dest='speedy', default=False,
                        help='[%(default)s] Do not do topology optimization.')
    parser.add_argument('--protein', action='store_true',
                        dest='protein', default=False,
                        help='[%(default)s] working with amino-acid sequences.')
    parser.add_argument('--sequential', action='store_true',
                        dest='sequential', default=False,
                        help='[%(default)s] Phylip sequential format.')
    parser.add_argument('--bitfast', action='store_true',
                        dest='medium', default=False,
                        help=\
                        '''[%(default)s] Same as fast, but reruns models with
                        topology optimization for best models (the ones with
                        cumulative weight => 0.95)''')
    parser.add_argument('--noinv', action='store_true',
                        dest='noinv', default=False,
                        help='[%(default)s] Do not check for invariant sites.')
    parser.add_argument('--nogam', action='store_true',
                        dest='nogam', default=False,
                        help='[%(default)s] Do not check for gamma distribution.')
    parser.add_argument('--nofrq', action='store_true',
                        dest='nofrq', default=False,
                        help='[%(default)s] Do not check for differences in rate frequencies.')
    parser.add_argument('--nprocs', metavar='INT',
                        dest='nprocs', default=2,
                        help='[%(default)s] Number of CPUs to use.')
    parser.add_argument('--sleep', metavar='FLOAT',
                        dest='refresh', default=0.1,
                        help='[%(default)s] Refresh rate in seconds.')
    parser.add_argument('--verbose', action='store_true',
                        dest='verb', default=False,
                        help='[%(default)s] Displays information about PhyML command line.')
    parser.add_argument('-m', metavar='LIST',
                        dest='models',
                        default='all',
                        choices=choose_model,
                        help= '[%(default)s] DNA/AA models. e.g.: -m "JC,TrN,GTR"')
    opts = parser.parse_args()
    typ = 'aa' if opts.protein else 'dna'
    if not opts.algt:
        exit(parser.print_help())
    if opts.medium:
        opts.speedy = True
    if opts.models == 'all':
        opts.models = model_list[typ]
    else:
        opts.models = opts.models.split(',')
    return opts


#########
# GLOBALS

global models
global freqs
global invts
global gamma
global modelnames

models = {'nt':
          ['000000',
           '010010',
           '010020',
           '012210',
           '010212',
           '012012',
           '012230',
           '010232',
           '012032',
           '012314',
           '012345'],
          'aa':
          ['LG'      ,
           'WAG'     ,
           'JTT'     ,
           'MtREV'   ,
           'Dayhoff' ,
           'DCMut'   ,
           'RtREV'   ,
           'CpREV'   ,
           'VT'      ,
           'BLOSUM2' ,
           'MtMam'   ,
           'MtArt'   ,
           'HIVw'    ,
           'HIVb'    ]
          }


freqs = {'nt': {'': ['-f', '0.25,0.25,0.25,0.25'], '+F': ['-f', 'm']},
         'aa': {'': ['-f', 'm'], '+F': ['-f', 'e']}}

invts = {'': ['-v', '0'], '+I': ['-v', 'e']}

gamma = {'': ['-c', '1'], '+G': ['-c', '4', '-a', 'e']}

# phyml model names, real names, and number of extra parameters
modelnames = { 'nt': { '000000' + ''    : ['JC'      , 0 ],
                       '010010' + ''    : ['K80'     , 1 ],
                       '010020' + ''    : ['TrNef'   , 2 ],
                       '012210' + ''    : ['TPM1'    , 2 ],
                       '010212' + ''    : ['TPM2'    , 2 ],
                       '012012' + ''    : ['TPM3'    , 2 ],
                       '012230' + ''    : ['TIM1ef'  , 3 ],
                       '010232' + ''    : ['TIM2ef'  , 3 ],
                       '012032' + ''    : ['TIM3ef'  , 3 ],
                       '012314' + ''    : ['TVMef'   , 4 ],
                       '012345' + ''    : ['SYM'     , 5 ],
                       '000000' + '+F'  : ['F81'     , 3 ],
                       '010010' + '+F'  : ['HKY'     , 4 ],
                       '010020' + '+F'  : ['TrN'     , 5 ],
                       '012210' + '+F'  : ['TPM1uf'  , 5 ],
                       '010212' + '+F'  : ['TPM2uf'  , 5 ],
                       '012012' + '+F'  : ['TPM3uf'  , 5 ],
                       '012230' + '+F'  : ['TIM1'    , 6 ],
                       '010232' + '+F'  : ['TIM2'    , 6 ],
                       '012032' + '+F'  : ['TIM3'    , 6 ],
                       '012314' + '+F'  : ['TVM'     , 7 ],
                       '012345' + '+F'  : ['GTR'     , 8 ]},
               'aa': { 'LG'       + ''  : ['LG'      , 0 ],
                       'WAG'      + ''  : ['WAG'     , 0 ],
                       'JTT'      + ''  : ['JTT'     , 0 ],
                       'MtREV'    + ''  : ['MtREV'   , 0 ],
                       'Dayhoff'  + ''  : ['Dayhoff' , 0 ],
                       'DCMut'    + ''  : ['DCMut'   , 0 ],
                       'RtREV'    + ''  : ['RtREV'   , 0 ],
                       'CpREV'    + ''  : ['CpREV'   , 0 ],
                       'VT'       + ''  : ['VT'      , 0 ],
                       'BLOSUM2'  + ''  : ['Blosum62', 0 ],
                       'MtMam'    + ''  : ['MtMam'   , 0 ],
                       'MtArt'    + ''  : ['MtArt'   , 0 ],
                       'HIVw'     + ''  : ['HIVw'    , 0 ],
                       'HIVb'     + ''  : ['HIVb'    , 0 ],
                       'LG'       + '+F': ['LG'      , 19],
                       'WAG'      + '+F': ['WAG'     , 19],
                       'JTT'      + '+F': ['JTT'     , 19],
                       'MtREV'    + '+F': ['MtREV'   , 19],
                       'Dayhoff'  + '+F': ['Dayhoff' , 19],
                       'DCMut'    + '+F': ['DCMut'   , 19],
                       'RtREV'    + '+F': ['RtREV'   , 19],
                       'CpREV'    + '+F': ['CpREV'   , 19],
                       'VT'       + '+F': ['VT'      , 19],
                       'BLOSUM2'  + '+F': ['Blosum62', 19],
                       'MtMam'    + '+F': ['MtMam'   , 19],
                       'MtArt'    + '+F': ['MtArt'   , 19],
                       'HIVw'     + '+F': ['HIVw'    , 19],
                       'HIVb'     + '+F': ['HIVb'    , 19]}
               }



if __name__ == "__main__":
    exit(main())
