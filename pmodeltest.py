#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/14 08:58:48

from optparse import OptionParser
from subprocess import Popen, PIPE
from re import sub
from sys import stderr as STDERR
from numpy import exp
#from consensus import consensus

__version__ = "0.50"
__title__   = "pmodeltest v%s" % __version__


# global variables
models = [ ['-m', '000000'],
           ['-m', '010010'],
           ['-m', '010020'],
           ['-m', '012210'],
           ['-m', '010212'],
           ['-m', '012012'],
           ['-m', '012230'],
           ['-m', '010232'],
           ['-m', '012032'],
           ['-m', '012314'],
           ['-m', '012345'] ]

freqs = { 'ef': ['-f', '0.25 0.25 0.25 0.25'],
          'uf': ['-f', 'e'                  ] }

invts = { ''  : [],
          '+I': ['-v', 'e'  ] }

gamma = { ''  : ['-a', '1.0'],
          '+G': ['-c', '4', '-a', 'e'] }

modelnames = { '000000' + 'ef': ['JC'     , 0],
               '010010' + 'ef': ['K80'    , 1],
               '010020' + 'ef': ['TrNef'  , 2],
               '012210' + 'ef': ['TPM1'   , 2],
               '010212' + 'ef': ['TPM2'   , 2],
               '012012' + 'ef': ['TPM3'   , 2],
               '012230' + 'ef': ['TIM1ef' , 3],
               '010232' + 'ef': ['TIM2ef' , 3],
               '012032' + 'ef': ['TIM3ef' , 3],
               '012314' + 'ef': ['TVMef'  , 4],
               '012345' + 'ef': ['SYM'    , 5],
               '000000' + 'uf': ['F81'    , 3],
               '010010' + 'uf': ['HKY'    , 4],
               '010020' + 'uf': ['TrN'    , 5],
               '012210' + 'uf': ['TPM1uf' , 5],
               '010212' + 'uf': ['TPM2uf' , 5],
               '012012' + 'uf': ['TPM3uf' , 5],
               '012230' + 'uf': ['TIM1'   , 6],
               '010232' + 'uf': ['TIM2'   , 6],
               '012032' + 'uf': ['TIM3'   , 6],
               '012314' + 'uf': ['TVM'    , 7],
               '012345' + 'uf': ['GTR'    , 8]
               }

def run_phyml(algt, wanted_models, speed, verb, rerun=False):
    '''
    runs a list of models and returns a dictionnary of results.
    If rerun is True, will have into account +I+G information
    and only run the wanted ones.
    '''
    if speed:
        opt = 'lr'
    else:
        opt = 'tlr'
    results = {}
    for model in models:
        for freq in freqs.keys():
            if not rerun and modelnames[model[1]+freq][0] not in \
                   sub('\+.*', '', wanted_models).split(','):
                continue
            for inv in invts.keys():
                for gam in gamma.keys():
                    if rerun:
                        if modelnames[model[1]+freq][0] + inv + gam \
                               not in wanted_models.split(','):
                            continue
                    model_name  = modelnames[model[1] + freq][0]
                    model_param = modelnames[model[1] + freq][1]
                    log = '\nModel ' + \
                          model_name + inv + gam + '\n'
                    log += 'Command line = ' +algt+ ' ' +\
                           ' '.join(model + freqs[freq] + invts[inv] + \
                                    gamma[gam]) + ' -t ' + opt + '\n'
                    (out, err) = Popen(['bin/phyml', '--sequential', 
                                        '-i', algt,
                                        '-d', 'nt',
                                        '-n', '1',
                                        '-b', '0',
                                        '-o', opt] + \
                                       model + freqs[freq] + \
                                       invts[inv] + gamma[gam],
                                       stdout=PIPE).communicate()
                    (numspe, lnl, dic) = parse_stats(algt + '_phyml_stats.txt')
                    tree          = get_tree   (algt + '_phyml_tree.txt') 
                    # num of param = X (nb of branches) + 1(topology) + Y(model)
                    numparam = model_param + \
                               (inv != '') + (gam != '') + numspe*2-3 + 1
                    aic = 2*numparam-2*lnl
                    log += 'K = '+str (numparam)+', lnL = '+str(lnl) + \
                           '\nAIC = ' + str (aic) + '\n'
                    log += '-----------------------------------\n'
                    if err is not None:
                        exit ('problem running phyml: '+out)
                    results[model_name + inv + gam] =  { 'AIC' : aic,
                                                         'lnL' : lnl,
                                                         'K'   : numparam,
                                                         'tree': tree,
                                                         'dic' : dic}
                    if verb:
                        print log

    return results


def aic_calc(results, speed):
    '''
    compute and displays AICs etc... 
    '''
    ord_aic = sorted (map (lambda x: [results[x]['AIC'], x], results.keys()))
    ord_aic = map (lambda x: x[1], ord_aic)
    min_aic = results[ord_aic[0]]['AIC']
    
    for model in ord_aic:
        results[model]['deltar'] =  results[model]['AIC'] - min_aic
        results[model]['weight'] = exp (-0.5 * results[model]['deltar'])

    sumweight = sum (map (lambda x: results[x]['weight'], results.keys()))

    cumweight = 0
    good_models = []
    for model in ord_aic:
        results[model]['weight'] = results[model]['weight']/sumweight
        cumweight += results[model]['weight']
        results[model]['cumweight'] = cumweight
        if results[model]['cumweight'] < 0.9999:
            good_models.append([model, int(1000*results[model]['weight']+0.5)])

    print '\n\n*************************************'
    print     '         TABLE OF WEIGHTS'
    if speed:
        print ' (WARNING: no topology optimization)'
    print     '*************************************\n'
    print '   MODEL       | AIC         | delta AIC   | Weights           | Cumulative weights'
    print '   ---------------------------------------------------------------------------------'
    print '\n'.join (map (lambda x: '   %-12s|'%(x) + \
                          ' %-11s |' % (str (results[x]['AIC'])       )+\
                          ' %-11s |' % (str (results[x]['deltar'])    ) +\
                          ' %-17s |' % (str (results[x]['weight'])    ) +\
                          ' %-15s' % (str (results[x]['cumweight']) ) \
                          , ord_aic))
    print '\n'
    return results, ord_aic


def main():
    '''
    main function when called by command line.
    infile must be in phylip format.
    '''
    opts = get_options()
    results = run_phyml(opts.algt, opts.models, opts.speedy, opts.verb)
    results, ord_aic = aic_calc(results, opts.speedy)

    if opts.medium:
        wanted_models = []
        for model in ord_aic:
            if results[model]['cumweight'] < 0.95:
                wanted_models.append(model)
            else:
                wanted_models.append(model)
                break
        wanted_models = ','.join(wanted_models)
        print '\nREFINING...\n    doing the same but computing topologies' + \
              ' only fo models that sums a weight > 0.95\n' + \
              wanted_models + '\n'
        results = run_phyml(opts.algt, wanted_models, \
                            False, opts.verb, rerun=True)
        results, ord_aic = aic_calc(results, False)

    print '\n\n*************************************************'
    print '  Tree corresponding to best model, '+ ord_aic[0] + '\n'
    print results[ord_aic[0]]['tree']


    #    filter (results[x][])
    #    results, log = run_phyml(opts.algt, opts.models, opts.speedy)
    #    if opts.verb:
    #        print log
    #print '\n\n Models keeped to build consensus: \n' + \
    #      ', '.join (map(lambda x: x[0], good_models))
    #
    #trees = []
    #tree_file = open(opts.algt + '_intree', 'w')
    #for model, weight in good_models:
    #    for i in range (weight):
    #        tree_file.write(results[model]['tree'])
    #        trees.append(results[model]['tree'])
    #tree_file.close()
    #
    #
    ##Popen(['yes | bin/consense'], stdout=PIPE)
    #
    ##final_tree   = get_tree(path + 'outtree')
    #better_model = ord_aic[0]
    #
    #print consensus(trees)
    # FINI!!!! YUJUUUUUUUUUUUUUUUUU

def parse_stats(path):
    '''
    parse stats file of phyml, to extract the likelyhood value
    '''
    dic = {}
    for line in open(path):
        if line.startswith('. Log-likelihood:'):
            lnl          = float (line.strip().split()[-1])
        elif line.startswith('. Number of taxa:'):
            numspe       = int (line.strip().split()[-1])
        elif line.startswith('  - f(A)= '):
            dic['fA']    = float (line.strip().split()[-1])
        elif line.startswith('  - f(T)= '):
            dic['fT']    = float (line.strip().split()[-1])
        elif line.startswith('  - f(G)= '):
            dic['fG']    = float (line.strip().split()[-1])
        elif line.startswith('  - f(C)= '):
            dic['fC']    = float (line.strip().split()[-1])
        elif line.startswith('  A <-> C '):
            dic['rAC']    = float (line.strip().split()[-1])
        elif line.startswith('  A <-> G'):
            dic['rAG']    = float (line.strip().split()[-1])
        elif line.startswith('  A <-> T'):
            dic['rAT']    = float (line.strip().split()[-1])
        elif line.startswith('  C <-> G'):
            dic['rCG']    = float (line.strip().split()[-1])
        elif line.startswith('  C <-> T'):
            dic['rCT']    = float (line.strip().split()[-1])
        elif line.startswith('  C <-> G'):
            dic['rCG']    = float (line.strip().split()[-1])
        elif line.startswith('. Proportion of invariant:'):
            dic['inv']   = float (line.strip().split()[-1])
        elif line.startswith('  - Gamma shape parameter:'):
            dic['gamma'] = float (line.strip().split()[-1])
    return (numspe, lnl, dic)

def get_tree(path):
    '''
    juste return the tree
    '''
    for line in open(path):
        if line.startswith('('):
            return line

def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]",
        description="""\
Reads sequeneces from file fasta format, and align acording to translation.                                      
.                                                                           .
********************************************                                      
TODO:                                                                                     
      (- fix averaging topology)                                                   
      (- averaging specific parameters)                                          
      - support for proteins                                                       
********************************************                                            
"""
        )

    model_list = ['JC', 'K80', 'TrNef', 'TPM1', 'TPM2', 'TPM3', 'TIM1ef', \
                  'TIM2ef', 'TIM3ef', 'TVMef', 'SYM', 'F81', 'HKY', 'TrN', \
                  'TPM1uf', 'TPM2uf', 'TPM3uf', 'TIM1', 'TIM2', 'TIM3', 'TVM', \
                  'GTR']
    parser.add_option('-i', dest='algt', metavar="PATH", \
                      help='path to input file in fasta format')
    parser.add_option('-o', dest='outfile', metavar="PATH", \
                      help='path to output file in fasta format')
    parser.add_option('-t', '--trimseqs', action='store_true', \
                      dest='trimseq', default=False, \
                      help='[%default] remove bad sequences (uses trimAl).')
    parser.add_option('--fast', action='store_true', \
                      dest='speedy', default=False, \
                      help='[%default] Do not do topology optimization.')
    parser.add_option('--bitfast', action='store_true', \
                      dest='medium', default=False, \
                      help=\
                      ''' [%default] Same as fast, but reruns models with
                      topology optimization for best models (the ones with
                      cumulative weight => 0.95)''')
    parser.add_option('--verbose', action='store_true', \
                      dest='verb', default=False, \
                      help=\
                      '''[%default] Displays information about PhyML command
                      line.''')
    parser.add_option('-m', metavar='LIST', \
                      dest='models', default=','.join(model_list), \
                      help=\
                      '''[%default] DNA models.                            
                      e.g.: -m "JC,TrN,GTR"
                      ''')
    opts = parser.parse_args()[0]
    if not opts.algt:
        exit(parser.print_help())
    if opts.medium:
        opts.speedy = True
    if len (set (opts.models.split(',')) - set (model_list)) > 0:
        print >>STDERR, 'ERROR: those models are not in list of ' + \
              'allowed models: \n   '+ \
               ', '.join (list (set (opts.models.split(',')) - \
                                set (model_list))) + '\n\n'
        exit(parser.print_help())
    return opts

if __name__ == "__main__":
    exit(main())
