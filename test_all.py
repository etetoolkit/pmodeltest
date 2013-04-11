__author__  = "francois serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "1.04"
__title__   = "pmodeltest v%s" % __version__

from re import sub
import pmodeltest as pmt

pmt.PHYML = '/usr/bin/phyml'

def test_nt():
    print '\n START TEST NUCLEOTIDS'
    print ' *********************\n\n'
    errors = []
    data  = '/home/francisco/Box/utils/pmodeltest/examples/dna.phy'
    models = ['JC', 'K80', 'TrNef', 'TPM1', 'TPM2', 'TPM3',
              'TIM1ef', 'TIM2ef', 'TIM3ef', 'TVMef', 'SYM',
              'F81', 'HKY', 'TrN', 'TPM1uf', 'TPM2uf', 'TPM3uf',
              'TIM1', 'TIM2', 'TIM3', 'TVM', 'GTR']
    job_list = pmt.get_job_list(data, models, speed=True, protein=False,
                            verbose=False)
    job_list = pmt.run_jobs(job_list, nprocs=2, refresh=0.01)
    job_list = pmt.parse_jobs(job_list, data)
    pmt.clean_all (job_list, data)
    job_list, ord_aic = pmt.aic_calc(job_list, True, verbose=True)
    if not ord_aic == ['TPM1uf+G+F','TPM3uf+G+F','TIM1+G+F','TPM1uf+I+G+F',
                       'TIM3+G+F','HKY+G+F','TPM3uf+I+G+F','TVM+G+F',
                       'TIM1+I+G+F','TrN+G+F','TIM3+I+G+F','HKY+I+G+F',
                       'TPM2uf+G+F','GTR+G+F','TVM+I+G+F','TPM1uf+I+F',
                       'TPM3uf+I+F','TrN+I+G+F','TIM2+G+F','TPM2uf+I+G+F',
                       'GTR+I+G+F','TIM1+I+F','TIM3+I+F','HKY+I+F','TIM2+I+G+F',
                       'TVM+I+F','TrN+I+F','TPM2uf+I+F','GTR+I+F','TIM2+I+F',
                       'TPM3uf+F','TPM1uf+F','HKY+F','TIM3+F','TIM1+F','TVM+F',
                       'TrN+F','TPM2uf+F','GTR+F','TIM2+F','F81+G+F',
                       'F81+I+G+F','SYM+G','SYM+I+G','TVMef+G','TVMef+I+G',
                       'SYM+I','TIM2ef+G','TIM2ef+I+G','TVMef+I','TPM2+G',
                       'TPM2+I+G','TIM3ef+G','TIM2ef+I','TIM1ef+G','TIM3ef+I+G',
                       'TIM1ef+I+G','TPM2+I','TPM3+G','TPM1+G','TIM3ef+I',
                       'TPM3+I+G','TPM1+I+G','TIM1ef+I','TrNef+G','TrNef+I+G',
                       'F81+I+F','TPM3+I','TPM1+I','K80+G','TrNef+I','K80+I+G',
                       'K80+I','F81+F','SYM','TIM2ef','TVMef','TIM3ef','TIM1ef',
                       'TPM2','TrNef','TPM3','TPM1','K80','JC+G','JC+I+G','JC+I','JC']:
        print 'Test ordering models: ERROR'
        errors.append('ordering models')
    else:
        print 'Test ordering models: OK'
    job_list = pmt.re_run(job_list, data, cutoff=0.95, refresh=0.01, nprocs=2,
                          verbose=True)
    job_list = pmt.parse_jobs(job_list, data)
    pmt.clean_all (job_list, data)
    if sorted(job_list.keys()) == sorted(['GTR+G+F','HKY+G+F','HKY+I+G+F','TIM1+G+F',
                                          'TIM1+I+G+F','TIM3+G+F','TIM3+I+G+F',
                                          'TPM1uf+G+F','TPM1uf+I+F','TPM1uf+I+G+F',
                                          'TPM2uf+G+F','TPM3uf+G+F','TPM3uf+I+F',
                                          'TPM3uf+I+G+F','TVM+G+F','TVM+I+G+F',
                                          'TrN+G+F','TrN+I+G+F']):
        print 'Test better models: OK'
    else:
        print 'Test better models: ERROR'
        errors.append('better models')
    job_list, ord_aic = pmt.aic_calc(job_list, False)
    if ord_aic == ['TPM1uf+G+F','TPM3uf+G+F','TIM1+G+F','TPM1uf+I+G+F','TIM3+G+F',
                   'HKY+G+F','TPM3uf+I+G+F','TVM+G+F','TIM1+I+G+F','TrN+G+F',
                   'TIM3+I+G+F','HKY+I+G+F','TPM2uf+G+F','GTR+G+F','TVM+I+G+F',
                   'TPM1uf+I+F','TPM3uf+I+F','TrN+I+G+F']:
        print 'Test ordering better models: OK'
    else:
        print 'Test ordering better models: ERROR'
        errors.append('ordering better models')
    
    expected = {'TrN+I+G+F'    : 6051.92212,
                'TPM1uf+I+G+F' : 6047.98258,
                'TVM+I+G+F'    : 6051.17612,
                'TIM1+I+G+F'   : 6049.6176 ,
                'TVM+G+F'      : 6049.37942,
                'HKY+G+F'      : 6048.42828,
                'TPM3uf+G+F'   : 6046.7577 ,
                'TIM3+G+F'     : 6048.40978,
                'TPM3uf+I+G+F' : 6048.5627 ,
                'TIM3+I+G+F'   : 6050.20106,
                'GTR+G+F'      : 6051.01224,
                'TrN+G+F'      : 6050.122  ,
                'TPM2uf+G+F'   : 6050.421  ,
                'TPM1uf+G+F'   : 6046.19414,
                'TPM3uf+I+F'   : 6051.40916,
                'TIM1+G+F'     : 6047.84258,
                'HKY+I+G+F'    : 6050.24122,
                'TPM1uf+I+F'   : 6051.35094}
    for j in expected:
        if not round(expected[j],2) == round(job_list[j]['AIC'],2):
            print 'Test AIC values: ERROR'
            errors.append('AIC values')
            break
    else:
        print 'Test AIC values: OK'


    tree = pmt.re_run_best(ord_aic[0], job_list[ord_aic[0]]['cmd'], data, verbose=True) 
    job_list = pmt.parse_jobs({ord_aic[0]: job_list[ord_aic[0]]}, data)
    pmt.clean_all ({ord_aic[0]: job_list[ord_aic[0]]}, data)
    
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

    print '\n\n\n TEST NUCLEOTIDES FINISHED\n\n'
    if errors:
        print 'ERROR founds:'
        for e in errors:
            print e
    else:
        print '  -> ALL TEST OK'

    print '\n\nas last check, have a look to jmodeltest result in test folder'
    #print ''.join([l for l in open('test/jmodeltest0.1_out.txt')])

def test_aa():
    print '\n START TEST AMINO-ACIDS'
    print ' **********************\n\n'
    errors = []
    data  = '/home/francisco/Box/utils/pmodeltest/examples/protein.phy'
    models = ['LG', 'WAG', 'JTT', 'MtREV', 'Dayhoff', 'DCMut',
              'RtREV', 'CpREV', 'VT', 'Blosum62', 'MtMam',
              'MtArt', 'HIVw', 'HIVb']
    job_list = pmt.get_job_list(data, models, speed=True, protein=True,
                            verbose=False) # here verbose prints all command lines
    job_list = pmt.run_jobs(job_list, nprocs=2, refresh=0.01)
    job_list = pmt.parse_jobs(job_list, data)
    pmt.clean_all (job_list, data)
    job_list, ord_aic = pmt.aic_calc(job_list, True, verbose=True)
    if not ord_aic == ['Dayhoff+G','DCMut+G','Dayhoff+I+G','DCMut+I+G',
                       'Dayhoff+I','DCMut+I','Dayhoff+G+F','DCMut+G+F',
                       'LG+G+F','Dayhoff+I+G+F','DCMut+I+G+F','LG+I+G+F',
                       'Dayhoff+I+F','DCMut+I+F','LG+I+F','WAG+G','WAG+I+G',
                       'WAG+I','Dayhoff','DCMut','WAG+G+F','RtREV+G+F','LG+G',
                       'WAG+I+G+F','WAG+I+F','RtREV+I+G+F','LG+I+G','RtREV+I+F',
                       'JTT+G+F','JTT+I+G+F','LG+I','JTT+I+F','JTT+G','JTT+I+G',
                       'Dayhoff+F','DCMut+F','CpREV+G+F','LG+F','JTT+I',
                       'CpREV+I+G+F','CpREV+I+F','WAG','WAG+F','MtArt+G+F',
                       'MtArt+I+G+F','RtREV+F','MtREV+G+F','JTT+F','MtREV+I+G+F',
                       'LG','HIVb+G+F','HIVb+I+G+F','MtREV+I+F','Blosum62+G+F',
                       'MtArt+I+F','JTT','Blosum62+I+G+F','Blosum62+I+F','CpREV+F',
                       'RtREV+G','HIVb+I+F','RtREV+I+G','CpREV+G','Blosum62+G',
                       'VT+G','RtREV+I','CpREV+I+G','Blosum62+I+G','VT+G+F',
                       'VT+I+G','CpREV+I','Blosum62+I','VT+I+G+F','VT+I','VT+I+F',
                       'MtMam+G+F','MtMam+I+G+F','Blosum62+F','MtREV+F','Blosum62',
                       'RtREV','MtMam+I+F','VT+F','VT','CpREV','MtArt+F','HIVb+F',
                       'HIVb+G','HIVb+I+G','HIVb+I','HIVw+G+F','HIVw+I+G+F',
                       'MtMam+F','HIVw+I+F','HIVb','MtREV+G','MtREV+I+G','MtArt+G','MtArt+I+G','MtREV+I','HIVw+F','MtArt+I','MtMam+G','MtMam+I+G','MtREV','HIVw+G','HIVw+I+G','MtMam+I','HIVw+I','MtArt','HIVw','MtMam']:
        print 'Test ordering models: ERROR'
        errors.append('ordering models')
    else:
        print 'Test ordering models: OK'
    job_list = pmt.re_run(job_list, data, cutoff=0.95, refresh=0.01, nprocs=2, verbose=True)
    job_list = pmt.parse_jobs(job_list, data)
    pmt.clean_all (job_list, data)
    if sorted(job_list.keys()) == sorted(['Dayhoff+G', 'DCMut+I+G', 'DCMut+G', 'Dayhoff+I+G']):
        print 'Test better models: OK'
    else:
        print 'Test better models: ERROR'
        errors.append('better models')
    job_list, ord_aic = pmt.aic_calc(job_list, False)
    if ord_aic == ['Dayhoff+G', 'DCMut+G', 'Dayhoff+I+G', 'DCMut+I+G']:
        print 'Test ordering better models: OK'
    else:
        print 'Test ordering better models: ERROR'
        errors.append('ordering better models')

    expected = {'DCMut+G'    : 3389.69506,
                'DCMut+I+G'  : 3391.46964,
                'Dayhoff+G'  : 3389.14838,
                'Dayhoff+I+G': 3390.9257}
    for j in expected:
        if not round(expected[j],2) == round(job_list[j]['AIC'],2):
            print 'Test AIC values: ERROR'
            print j, job_list[j]['AIC'], round(job_list[j]['AIC'],2)
            errors.append('AIC values')
            break
    else:
        print 'Test AIC values: OK'
    
    tree = pmt.re_run_best(ord_aic[0], job_list[ord_aic[0]]['cmd'], data, verbose=True)
    pmt.clean_all ({ord_aic[0]: job_list[ord_aic[0]]}, data)

    expected_tree = '((human:0.0563158672,(rat:0.2056794657,rabbit:0.0839307491)0.7670000000:0.0369524233)0.6450000000:0.0289065340,marsupial:0.2993530194,cow:0.1082479462);\n'
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

    print '\n\n\n TEST PROTEIN FINISHED\n\n'
    if errors:
        print 'ERROR founds:'
        for e in errors:
            print e
    else:
        print '  -> ALL TEST OK'
    
    print '\n\nas last check, have a look to jmodeltest result in test folder'
    #print ''.join([l for l in open('test/prottest3_out.txt')])
    
    
def main():
    test_nt()
    test_aa()
            
if __name__ == "__main__":
    exit(main())
