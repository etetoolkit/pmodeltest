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

"""

__author__  = "francois serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "1.03"
__title__   = "pmodeltest v%s" % __version__

import os
from optparse import OptionParser
from subprocess import Popen, PIPE
from re import sub
from sys import stderr as STDERR
from sys import stdout as STDOUT
from cmath import exp

def extend_dict(source, target):
    for k,v in source.iteritems():
        target[k] = v

def run_parallel(cmds):
    for c in cmds:
        os.system(c)

def get_test_thread(wanted_models, speed, protein,
                    support, sequential=True):
    '''Generate a list of phyml executions needed to perform the whole
    test. The complete test thread is returned as a list of
    dictionaries. Each dictionary contains the needed phyml arguments.
    '''
    opt = 'l'  if speed   else 'tl'
    sup = '-4' if support else '0'
    typ = 'aa' if protein else 'nt'
    test_thread = []

    for model in models[typ]:
        for freq in freqs[typ].keys():
            mname  = modelnames[typ][model + freq][0]
            mparams = modelnames[typ][model + freq][1]
            if modelnames[typ][model+freq][0] not in wanted_models:
                continue

            for inv in invts.keys():
                for gam in gamma.keys():
                    phyml_args = {
                        '-d': 'aa' if protein else 'nt',
                        '-n': '1',
                        '-b': sup,
                        '-o': opt,
                        '-m': model,
                    }
                    extend_dict(freqs[typ][freq], phyml_args)
                    extend_dict(invts[inv], phyml_args)
                    extend_dict(gamma[gam], phyml_args)
                    test_thread.append([phyml_args, mname, mparams, freq, inv, gam, opt])
    return test_thread

def get_cmd(args):
    cmd = "phyml " + ' '.join(["%s %s" %(k,v)
                               for k,v in args.iteritems()])
    return cmd 

def iter_test_commands(alg_path, test_thread):
    for args, mname, mparams, freq, inv, gam, opt in test_thread:
        test_alg_path = alg_path + ".%s_%s_%s_%s" %(mname,freq,inv,gam)

        os.symlink(os.path.basename(alg_path), test_alg_path)
        
        args["-i"] = test_alg_path
        yield get_cmd(args)

def collect_results(alg_path, test_thread):
    results = {}
    for args, mname, mparams, freq, inv, gam, opt in test_thread:
        test_alg_path = alg_path + ".%s_%s_%s_%s" %(mname,freq,inv,gam)

        numspe, lnl, dic = parse_stats(test_alg_path + '_phyml_stats.txt')
        # num of param = X (nb of branches) + 1(topology) + Y(model)
        numparam = mparams + int (opt=='tl') + \
                   (inv != '') + (gam != '') + numspe*2-3
        aic = 2 * numparam-2 * lnl
        #log += '  K = '+str (numparam)+', lnL = '+str(lnl) + \
        #       ', AIC = ' + str (aic)
        #if err is not None or 'Err: ' in out:
        #    exit ('ERROR: problem running phyml: '+out)
        results[mname + inv + gam + freq] =  {
            'AIC' : aic,
            'lnL' : lnl,
            'K'   : numparam,
            'dic' : dic,
            'tree': get_tree(test_alg_path + '_phyml_tree.txt'),
            'cmnd': args}
        #if verb:
        #    print >> STDOUT, log
    return results

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


def main():
    '''
    main function when called by command line.
    infile must be in phylip format.
    '''
    opts = get_options()
    # remove gamma inv and frequencies if not wanted
    if opts.nogam:
        del gamma['+G']
    if opts.noinv:
        del invts['+I']
    if opts.nofrq:
        del freqs['nt']['+F']
        del freqs['aa']['+F']

    # first run of models
    #results = run_phyml(opts.algt, opts.models, opts.speedy, opts.verb,
    #                    opts.protein, opts.support, sequential=opts.sequential)

    alg_path = os.path.realpath(opts.algt)
    if not os.path.exists("pmodeltest_results"):
        os.mkdir("pmodeltest_results")
    alg_link_path = os.path.join("pmodeltest_results",
                                 os.path.basename(opts.algt))
    if not os.path.exists(alg_link_path):
        os.symlink(alg_path, alg_link_path)
        print alg_path, alg_link_path

    wanted_models = set(sub('\+.*', '', opts.models).split(','))        
    test_thread = get_test_thread(wanted_models, opts.speedy,
                                  opts.protein, opts.support, opts.sequential)
    run_parallel([cmd for cmd in iter_test_commands(alg_link_path, test_thread)])
    results = collect_results(alg_link_path, test_thread)
    results, ord_aic = aic_calc(results, opts.speedy)
    print ord_aic

    # if bit fast, second run with wanted models (that sums weight of 0.95)
    if opts.medium:
        best_results = {}
        for model in ord_aic:
            best_results[model] = results[model]
            if results[model]['cumweight'] > 0.95:
                break

        print >> STDOUT,  '\nREFINING...\n    doing the same but computing topologies' + \
              ' only for models that sums a weight of 0.95\n\n    '

        refine_cmds = []
        for r in best_results.itervalues():
            args = r["cmnd"]
            args["-o"] = "lr"
            refine_cmds.append(get_cmd(args))
        run_parallel(refine_cmds)
        results, ord_aic = aic_calc(best_results, False)

    print >> STDOUT,  '\n\n*************************************************'
    #results[ord_aic[0]]['cmnd'][results[ord_aic[0]]['cmnd'].index ('-o') + 1] += 'r'
    winner_args = results[ord_aic[0]]['cmnd']
    winner_args["-o"] = "tlr"
    winner_args["-b"] = "-4"
    winner_args["-s"] = "BEST"
    winner_args["-i"] = alg_path
    run_parallel([get_cmd(winner_args)])
    print >> STDOUT,\
          'Re-run of best model with computation of rates and support...'
    
    tree = get_tree   (alg_path + '_phyml_tree.txt')
    print >> STDOUT, '\n Corresponding estimations of rates/frequencies:\n'
    print_model_estimations (parse_stats (alg_path + '_phyml_stats.txt')[2])
    print >> STDOUT, '\nTree corresponding to best model, '\
          + ord_aic[0] + ' (with SH-like branch supports alone)\n'
    print >> STDOUT,  tree
    if opts.outfile:
        open (opts.outfile, 'w').write (tree)
    if opts.outtrees:
        out_t = open (opts.outtrees, 'w')
        for run in results:
            out_t.write ('command: ' + \
                         ' '.join (results [run]['cmnd']) + \
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
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]",
        description="""\
Reads sequences from file fasta format, and align according to translation.                                      
.                                                                           .
********************************************                                            
"""
        )

    model_list = {'dna': ['JC', 'K80', 'TrNef', 'TPM1', 'TPM2', 'TPM3',
                          'TIM1ef', 'TIM2ef', 'TIM3ef', 'TVMef', 'SYM',
                          'F81', 'HKY', 'TrN', 'TPM1uf', 'TPM2uf', 'TPM3uf',
                          'TIM1', 'TIM2', 'TIM3', 'TVM', 'GTR'],
                  'aa' :  ['LG', 'WAG', 'JTT', 'MtREV', 'Dayhoff', 'DCMut',
                           'RtREV', 'CpREV', 'VT', 'Blosum62', 'MtMam',
                           'MtArt', 'HIVw', 'HIVb'] }
    parser.add_option('-i', dest='algt', metavar="PATH", \
                      help='path to input file in phylip format')
    parser.add_option('-o', dest='outfile', metavar="PATH", \
                      help='name of outfile tree (newick format)')
    parser.add_option('-O', dest='outtrees', metavar="PATH", \
                      help='name of outfile with all trees (newick format)')
    parser.add_option('--support', action='store_true', \
                      dest='support', default=False, \
                      help='[%default] compute SH-like branch support for each model (slower).')
    parser.add_option('--fast', action='store_true', \
                      dest='speedy', default=False, \
                      help='[%default] Do not do topology optimization.')
    parser.add_option('--protein', action='store_true', \
                      dest='protein', default=False, \
                      help='[%default] working with amino-acid sequences.')
    parser.add_option('--sequential', action='store_true', \
                      dest='sequential', default=False, \
                      help='[%default] Phylip sequential format.')
    parser.add_option('--bitfast', action='store_true', \
                      dest='medium', default=False, \
                      help=\
                      '''[%default] Same as fast, but reruns models with
                      topology optimization for best models (the ones with
                      cumulative weight => 0.95)''')
    parser.add_option('--noinv', action='store_true', \
                      dest='noinv', default=False, \
                      help='[%default] Do not check for invariant sites.')
    parser.add_option('--nogam', action='store_true', \
                      dest='nogam', default=False, \
                      help='[%default] Do not check for gamma distribution.')
    parser.add_option('--nofrq', action='store_true', \
                      dest='nofrq', default=False, \
                      help='[%default] Do not check for differences in rate frequencies.')
    parser.add_option('--verbose', action='store_true', \
                      dest='verb', default=False, \
                      help=\
                      '''[%default] Displays information about PhyML command
                      line.''')
    parser.add_option('-m', metavar='LIST',
                      dest='models',
                      default=(' '*80).join([ m  +': ' + \
                                              ','.join (model_list[m]) \
                                              for m in model_list ]),
                      help=\
                      '''[%default] DNA/AA models.                            
                      e.g.: -m "JC,TrN,GTR"
                      ''')
    opts = parser.parse_args()[0]
    typ = 'aa' if opts.protein else 'dna'
    if not opts.algt:
        exit(parser.print_help())
    if opts.medium:
        opts.speedy = True
    if opts.models == (' '*80).join([ m  +': ' + ','.join(model_list[m]) \
                                      for m in model_list ]):
        opts.models = ','.join (model_list[typ])
    if len (set (opts.models.split(',')) - set (model_list[typ])) > 0:
        print >> STDERR, 'ERROR: those models are not in list of ' + \
              'allowed models: \n   '+ \
               ', '.join (list (set (opts.models.split(',')) - \
                                set (model_list[typ]))) + '\n\n'
        exit(parser.print_help())
    return opts


#########
# GLOBALS

global models
global freqs
global invts
global gamma
global modelnames

models = {'nt':
          [ '000000',
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
          [ 'LG',
            'WAG',
            'JTT',
            'MtREV',
            'Dayhoff',
            'DCMut',
            'RtREV',
            'CpREV',
            'VT',
            'Blosum62',
            'MtMam',
            'MtArt',
            'HIVw',
            'HIVb']
          }


freqs = {'nt': {'': {'-f': '0.25 0.25 0.25 0.25'}, '+F': {'-f':'m'}},
         'aa': {'': {'-f': 'm'}, '+F': {'-f': 'e'}}}

invts = {'': {}, '+I': {'-v': 'e'}}

gamma = {'': {'-c': '1', '-a': '1.0'}, '+G': {'-c': '4', '-a': 'e'}}

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
                       'Blosum62' + ''  : ['Blosum62', 0 ],
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
                       'Blosum62' + '+F': ['Blosum62', 19],
                       'MtMam'    + '+F': ['MtMam'   , 19],
                       'MtArt'    + '+F': ['MtArt'   , 19],
                       'HIVw'     + '+F': ['HIVw'    , 19],
                       'HIVb'     + '+F': ['HIVb'    , 19]}
               }



if __name__ == "__main__":
    exit(main())
