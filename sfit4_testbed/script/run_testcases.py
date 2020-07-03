#!/usr/bin/python2.7

# goes through all testcases, runs hbin, sfit4 and prints out some
# numbers.  Intended to test a setup on a new computer, new compiler,
# new compiler flags or after code changes.
#
# test.ctl is required and contains some information for running the test cases.
#
# Mathias Palm <mathias.palm@uni-bremen.de> January 2019

import os,sys,string

from libs import sfit4_ctl,summary,statevec, read_from_file
from shutil import copy
from pathlib2 import Path

import numpy as np
import subprocess

class test_sfit4:

    def __init__(self, cfgfile):

        fid = open(cfgfile)

        self.results = {}
        self.results_orig = {}
        
        for l in fid:
            l = l.split('#')[0]
            key = l.split('=')[0].strip()
            subkeys = key.rsplit('.')
            if len(key) == 0:
                continue
            if len(l) < 2:
                print ('No value for key '+key+'?')
            if key.lower() == 'sfit4_dir':
                self.sfit4_dir = l.rsplit('=')[1].strip()
                if self.sfit4_dir == '<SFIT4-DIR>' or not Path(self.sfit4_dir).is_dir():
                    self.sfit4_dir = raw_input('Please specify the sfit-core-code directory\n')
                while not Path(self.sfit4_dir).is_dir():
                    self.sfit4_dir = raw_input('{} does not exist or is no directory. Please try again.\n>'.format(self.sfit4_dir))
                self.sfit4_dir += '/'
                continue
            if key.lower() =='linelist_dir':
                self.linelist_dir = l.rsplit('=')[1].strip()
                if self.linelist_dir == '<LINELIST-DIR>' or not Path(self.linelist_dir).is_dir():
                    self.linelist_dir = str(raw_input('Please specify the linelistdir\n'))
                while not Path(self.linelist_dir).is_dir():
                    self.linelist_dir = raw_input('{} does not exist or is no directory. Please try again.\n>'.format(self.linelist_dir))
                self.linelist_dir += '/'
                continue
            if key.lower() == 'testdir':
                self.testcase_dir = l.rsplit('=')[1].strip()
                if self.testcase_dir == '<TESTBED-DIR>' or not Path(self.testcase_dir).is_dir():
                    self.testcase_dir = raw_input('Please specify the testcase directory\n')
                while not Path(self.testcase_dir).is_dir():
                    self.testcase_dir = raw_input('{} does not exist or is no directory. Please try again.\n>'.format(self.testcase_dir))
                self.testcase_dir += '/'
                continue
            if key.lower() =='origtestcases_dir':
                ll = l.rsplit('=')[1].strip()+'/'
                self.origtestcases_dir = os.path.join(self.testcase_dir,ll)
                continue
            if key.lower() == 'resultfile':
                ll = l.rsplit('=')[1].strip()
                self.resultfile = os.path.join(self.testcase_dir,ll)
                continue
            if key.lower() == 'hbinfile':
                ll = l.rsplit('=')[1].strip()
                self.hbinfile = os.path.join(self.testcase_dir,ll)
                continue
            if key.lower() == 'tips':
                self.tips = l.rsplit('=')[1].strip()
                continue
            if key.lower() == 'gas':
                self.gases = l.rsplit('=')[1:][0].split()
                print (l)
                for tc in self.gases:
                    tc = tc.strip()
                    self.results.update({tc:{}})
                    self.results_orig.update({tc:{}})



            if len(subkeys) > 1:

                if subkeys[0].strip() == 'gas' and len(self.results.keys()) == 0:
                    print('Gases must be defined before details')
                    return(None)
                gas = subkeys[1].strip()
                if self.results.keys().count(gas) == 0:
                    continue
                for sk in subkeys:
                    lval = False
                    if l.rsplit('=')[1:][0].strip() == 'T':
                        lval = True
                    if sk.lower() == 'dir':
                        self.results[gas].update({'dir':l.rsplit('=')[1:][0].strip()})
                    if sk.lower() == 'run_hbin':
                        self.results[gas].update({'hbin':lval})
                    if sk.lower() == 'run_sfit4':
                        self.results[gas].update({'sfit':lval})


        extra_gas = filter(lambda x: self.gases.count(x) == 0, self.results.keys())

        if len(extra_gas) > 0:
            print ('No details found for gas(es) '+ string.join(extra_gas))
        
    def run_sfit4_in_testcase(self, sfit4=True, hbin=True, tips=False, error = True):

        ctl = sfit4_ctl()
        hbin_ctl = sfit4_ctl()
        
        for tc in self.results.keys():
            print('Entering %s'%(self.results[tc]['dir']))            
            os.chdir(os.path.join(self.testcase_dir, self.results[tc]['dir']))
            
            flag = False
            if hbin and self.results[tc]['hbin']:
                print('Calling hbin')
                copy(self.hbinfile,'.')
                hbin_ctl.replace_in_file('hbin.ctl','file.in.linelist', self.linelist_dir)
                chbin = os.path.join(self.sfit4_dir,'src','hbin')
                if not Path(chbin).is_file():
                    print('{} does not exist'.format(chbin))
                    exit()
                    
                rhbin = subprocess.Popen(chbin,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
                if len(rhbin[1]):
                    print ('hbin failed in %s'%(tc))

                for l in rhbin[0].split('\n'):
                    if l.find('Saving binary')>-1:
                        flag = True
                        fhbin = l.split(':')[1].strip()
                        print('hbin file in %s: %s')%(tc,fhbin)
                if not flag:
                    print('Something wrong in %s'%tc)
                    print('hbin terminated with:')
                    print(rhbin[0])
                    self.results[tc].update({hbin:False})
                    
                ctl.replace_in_file('sfit4.ctl','file.in.linelist',fhbin)

            if sfit4 and self.results[tc]['sfit']:
                if error:
                    ctl.replace_in_file('sfit4.ctl','kb','T')
                else:
                    ctl.replace_in_file('sfit4.ctl','kb','F')
                if tips:
                    ctl.replace_in_file('sfit4.ctl','fw.tips','T')
                else:
                    ctl.replace_in_file('sfit4.ctl','fw.tips','F')
                print('Calling sfit4')
                csfit = os.path.join(self.sfit4_dir,'src','sfit4')
                rsfit = subprocess.Popen(csfit,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
                flag = False
                for l in rsfit[0].split('\n'):
                    if l.find('CVRG')>-1:
                        flag = True
                        if l[l.index('CVRG')+5]=='T':
                            print('sfit converged in %s'%tc)
                        else:
                            print('sfit did not converge in %s'%tc)
                            self.results[tc].update({sfit4:False})

                if not flag:
                    print('Something wrong in %s'%tc)
                    print('sfit4 terminated with:')
                    print(rsfit[1])
                    self.results[tc].update({sfit4:False})
                            #
#            sf4l0args = ['-i',tcpath,
#                         '-b',os.path.join(self.sfit4_dir,'src'),
#                         '-f',flag]
#            sfit4Layer0(sf4l0args)


    def read_summaries(self):
        
        # Store values from summary
        for tc in self.results:

            sum_file = os.path.join(self.origtestcases_dir,
                                    'summary.%s'%(tc.lower()))
            if os.path.exists(sum_file):
                sum_orig = summary(sum_file)
                result = {'summary': True,
                          'apriori':sum_orig.apriori[0],
                          'retriev':sum_orig.retriev[0],
                          'chi_y_2':sum_orig.chi_y_2
                }
                self.results_orig[tc].update(result)
            else:
                self.results_orig[tc].update({'summary':False})

            sum_file = os.path.join(self.testcase_dir,
                                    self.results[tc]['dir'],'summary')
            if os.path.exists(sum_file):
                sum_new = summary(sum_file)
                result = {'summary': True,
                          'apriori':sum_new.apriori[0],
                          'retriev':sum_new.retriev[0],
                          'chi_y_2':sum_new.chi_y_2,
                          'converged':sum_new.converged
                }
                self.results[tc].update(result)
            else:
                self.results[tc].update({'summary':False})

    def read_statevectors(self):
        #Store values from statevec
        state_orig = statevec(os.path.join(self.origtestcases_dir,'statevec.%s'%(tc[1])))
        result = {'ret_profile': state_orig.rt_vmr[0]}
        results_orig[tc[0]].update(result)
        
        state = statevec(os.path.join(tcpath,'statevec'))
        result = {'ret_profile': state.rt_vmr[0]}
        results[tc[0]].update(result)

    def print_summary(self):
        str = ''
        for rs in self.results.keys():
            if self.results[rs]['converged']:
                str += 'Testcase {0}: RUN OK '.format(rs)
            else:
                str += 'Testcase {0}: RUN NOT OK '.format(rs)

            diverge = self.results[rs]['chi_y_2'] - self.results_orig[rs]['chi_y_2']
            diverge *= 2
            diverge /= self.results[rs]['chi_y_2'] + self.results_orig[rs]['chi_y_2']
            if diverge < 0.01:
                str += 'RESULTS OK \n'.format(rs)
            else:
                str += 'CHI_2_Y DIVERGES BY {1:1%} %\n'.format(rs, diverge)

        str += '\n'
        print(str)
        
    def print_results(self):

        with open(self.resultfile, 'w') as fid:
            fid.write('Quantity\tThis run\t Original run\t Difference in Percent\n')

            for rs in self.results.keys():
                if self.results[rs]['summary'] and self.results_orig[rs]['summary']:
                    fid.write('Gas {}\n'.format(rs))
                    fid.write('Apriori\t\t')
                    fid.write('{}\t'.format(self.results[rs]['apriori']))
                    fid.write('{}\t'.format(self.results_orig[rs]['apriori']))
                    fid.write('{} %\n'.format(200*(self.results[rs]['apriori']-self.results_orig[rs]['apriori'])/(self.results[rs]['apriori']+self.results_orig[rs]['apriori'])))
                    
                    fid.write('Retrieved:\t')
                    fid.write('{}\t'.format(self.results[rs]['retriev']))
                    fid.write('{}\t'.format(self.results_orig[rs]['retriev']))
                    fid.write('{} %\n'.format(200*(self.results[rs]['retriev']-self.results_orig[rs]['retriev'])/(self.results[rs]['retriev']+self.results_orig[rs]['retriev'])))

                    
                    fid.write('CHI_Y_2:\t')
                    fid.write('{}\t'.format(self.results[rs]['chi_y_2']))
                    fid.write('{}\t'.format(self.results_orig[rs]['chi_y_2']))
                    fid.write('{} %\n'.format(200*(self.results[rs]['chi_y_2']-self.results_orig[rs]['chi_y_2'])/(self.results[rs]['chi_y_2']+self.results_orig[rs]['chi_y_2'])))
                    #            print 'MEAN SQARE Diff. RETRIEVED VMR:', np.sqrt(np.mean((self.results[rs]['chi_y_2']-self.results_orig[rs]['chi_y_2'])**2))
                else:
                    fid.write('No new summary or original summary file found')
    
                
if __name__ == '__main__':

    tc = test_sfit4('test.cfg')
    runsfit = True
    runhbin = True
    error = True
    tips = True
    if sys.argv.count('--nosfit4') > 0:
        runsfit = False
    if sys.argv.count('--nohbin') > 0:
        runhbin = False
    if sys.argv.count('--noerror') > 0:
        error = False
    if sys.argv.count('--notips') > 0:
        tips = False
    script_path = os.path.dirname(os.path.realpath(__file__))
    tc.run_sfit4_in_testcase(sfit4=runsfit,hbin=runhbin,tips=tips,error=error)
    tc.read_summaries()
    tc.print_summary()
 #   tc.read_statevectors()
    tc.print_results()
