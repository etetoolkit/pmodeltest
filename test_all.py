__author__  = "francois serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "1.03"
__title__   = "pmodeltest v%s" % __version__

from subprocess import Popen, PIPE
from re import sub
from pmodeltest import *


def main():
    errors = []
    data  = '/home/francisco/Box/utils/pmodeltest/examples/dna.phy'
    models = ['JC', 'K80', 'TrNef', 'TPM1', 'TPM2', 'TPM3',
              'TIM1ef', 'TIM2ef', 'TIM3ef', 'TVMef', 'SYM',
              'F81', 'HKY', 'TrN', 'TPM1uf', 'TPM2uf', 'TPM3uf',
              'TIM1', 'TIM2', 'TIM3', 'TVM', 'GTR']
    job_list = get_job_list(data, models, speed=True, protein=False,
                            support=False)
    job_list = run_jobs(job_list, nprocs=2, refresh=0.01)
    job_list = parse_jobs(job_list, data)
    job_list, ord_aic = aic_calc(job_list, True)
    if not ord_aic == ["TPM1uf+G+F","TPM3uf+G+F","TIM1+G+F","TPM1uf+I+G+F","TIM3+G+F","HKY+G+F","TPM3uf+I+G+F","TVM+G+F","TIM1+I+G+F","TrN+G+F","TIM3+I+G+F","HKY+I+G+F","TPM2uf+G+F","GTR+G+F","TVM+I+G+F","TPM1uf+I+F","TPM3uf+I+F","TrN+I+G+F","TIM2+G+F","TPM2uf+I+G+F","GTR+I+G+F","TIM1+I+F","TIM3+I+F","HKY+I+F","TIM2+I+G+F","TVM+I+F","TrN+I+F","TPM2uf+I+F","GTR+I+F","TIM2+I+F","TPM3uf+F","TPM1uf+F","HKY+F","TIM3+F","TIM1+F","TVM+F","TrN+F","TPM2uf+F","GTR+F","TIM2+F","F81+G+F","F81+I+G+F","SYM+G","SYM+I+G","F81+I+F","TVMef+G","TVMef+I+G","SYM+I","TIM2ef+G","TIM2ef+I+G","TVMef+I","TPM2+G","TPM2+I+G","TIM3ef+G","TIM2ef+I","TIM1ef+G","TIM3ef+I+G","TIM1ef+I+G","TPM2+I","TPM3+G","TPM1+G","TIM3ef+I","TPM3+I+G","TPM1+I+G","TIM1ef+I","TrNef+G","TrNef+I+G","TPM3+I","TPM1+I","K80+G","TrNef+I","K80+I+G","K80+I","F81+F","SYM","TIM2ef","TVMef","TIM3ef","TIM1ef","TPM2","TrNef","TPM3","TPM1","K80","JC+G","JC+I+G","JC+I","JC"]:
        print 'Test ordering models: ERROR'
        errors.append('ordering models')
    else:
        print 'Test ordering models: OK'
    job_list = re_run(job_list, data, cutoff=0.95, refresh=0.01, nprocs=2)
    if sorted(job_list.keys()) == sorted(["TIM3+G+F","TrN+I+G+F","TPM3uf+I+F","TIM1+G+F","GTR+G+F","TPM1uf+I+F","TPM3uf+I+G+F","TVM+G+F","TrN+G+F","TPM1uf+G+F","TPM2uf+G+F","TVM+I+G+F","HKY+I+G+F","TIM1+I+G+F","TPM3uf+G+F","TIM3+I+G+F","TPM1uf+I+G+F","HKY+G+F"]):
        print 'Test better models: OK'
    else:
        print 'Test better models: ERROR'
        errors.append('better models')
    job_list, ord_aic = aic_calc(job_list, False)
    if ord_aic == ['TPM1uf+G+F','TPM3uf+G+F','TIM1+G+F','TPM1uf+I+G+F','TIM3+G+F','HKY+G+F','TPM3uf+I+G+F','TVM+G+F','TIM1+I+G+F','TrN+G+F','TIM3+I+G+F','HKY+I+G+F','TPM2uf+G+F','GTR+G+F','TVM+I+G+F','TPM1uf+I+F','TPM3uf+I+F','TrN+I+G+F']:
        print 'Test ordering better models: OK'
    else:
        print 'Test ordering better models: ERROR'
        errors.append('ordering better models')
    
    expected = {'GTR+G+F'     : 6051.0123,
                'HKY+G+F'     : 6048.42834,
                'HKY+I+G+F'   : 6050.24228,
                'TIM1+G+F'    : 6047.84272,
                'TIM1+I+G+F'  : 6049.61996,
                'TIM3+G+F'    : 6048.40994,
                'TIM3+I+G+F'  : 6050.20164,
                'TPM1uf+G+F'  : 6046.19426,
                'TPM1uf+I+F'  : 6051.35108,
                'TPM1uf+I+G+F': 6047.98534,
                'TPM2uf+G+F'  : 6050.42124,
                'TPM3uf+G+F'  : 6046.75778,
                'TPM3uf+I+F'  : 6051.40928,
                'TPM3uf+I+G+F': 6048.56366,
                'TVM+G+F'     : 6049.37956,
                'TVM+I+G+F'   : 6051.17726,
                'TrN+G+F'     : 6050.12206,
                'TrN+I+G+F'   : 6051.92446}
    
    for j in expected:
        if not round(expected[j],2) == round(job_list[j]['AIC'],2):
            print 'Test AIC values: ERROR'
            errors.append('AIC values')
            break
    else:
        print 'Test AIC values: OK'
    
    cmd = job_list[ord_aic[0]]['cmd']
    cmd [cmd.index ('-b')+1] = '-4'
    cmd += ['-s', 'BEST']
    (out, err) = Popen('echo "end" | ' + ' '.join (cmd),
                       stdout=PIPE, shell=True).communicate()
    tree = get_tree   (data + '_phyml_tree_%s.txt' % ord_aic[0])
    print_model_estimations (parse_stats (data + '_phyml_stats_%s.txt' % ord_aic[0])[2])
    
    expected_tree = '((((((Macaca_mulatta:0.0117313526,(Alouata_seniculus:0.0190656000,(Pongo_pygmaeus_old_3:0.0030481589,Gorilla_gorilla:0.0020294397)0.9770000000:0.0075362735)0.0000000000:0.0004306512)0.9960000000:0.0172550253,Pan_troglodytes:0.0223821735)0.8640000000:0.0093079991,Trachypithecus_johnii:0.0229373226)1.0000000000:0.0776907115,Procolobus_badius:0.1966057310)0.9970000000:0.0650031999,Trachypithecus_geii:0.0048258029)0.8290000000:0.0030733844,Semnopithecus_entellus:0.0052455542,Callithrix_jacchus:0.0081475480);\n'
    expected_tree = sub('(:[0-9]+\.[0-9]{3})[0-9]*','\\1', expected_tree)
    expected_tree = sub('(\))[0-9]+\.[0-9]*','\\1', expected_tree)
    tree = sub('(:[0-9]+\.[0-9]{3})[0-9]*','\\1', tree)
    tree = sub('(\))[0-9]+\.[0-9]*','\\1', tree)
    if expected_tree == tree:
        print 'Testing final topology: OK'
    else:
        print 'Testing final topology: ERROR'
        errors.append('Different final trees')
        print expected_tree
        print tree

    print '\n\n\n TEST FINISHED\n\n'
    if errors:
        print 'ERROR founds:'
        for e in errors:
            print e
            
if __name__ == "__main__":
    exit(main())
