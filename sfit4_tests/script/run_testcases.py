#!/usr/bin/python 

# goes through all testcases, runs hbin, sfit4 and prints out some
# numbers.  Intended to test a setup on a new computer, new compiler,
# new compiler flags or after code changes.
#
# test.ctl is required and contains some information for running the test cases.
#
# Mathias Palm <mathias.palm@uni-bremen.de> January 2019

import os,sys,string

from libs import sfit4_ctl,summary,statevec, read_from_file

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
                print 'No value for key '+key+'?'
            if key.lower() == 'sfit4_dir':
                self.sfit4_dir = l.rsplit('=')[1].strip()
                continue
            if key.lower() == 'testdir':
                self.testcase_dir = l.rsplit('=')[1].strip()
                continue
            if key.lower() == 'resultfile':
                self.resultfile = l.rsplit('=')[1].strip()
                continue
            if key.lower() == 'gas':
                self.gases = l.rsplit('=')[1:][0].split()
                print l
                for tc in self.gases:
                    tc = tc.strip()
                    self.results.update({tc:{}})
                    self.results_orig.update({tc:{}})



            if len(subkeys) > 1:

                if subkeys[0].strip() == 'gas' and len(self.results.keys()) == 0:
                    print 'Gases must be defined before details'
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
            print 'No details found for gas(es) '+ string.join(extra_gas)
        
    def run_sfit4_in_testcase(self, sfit4=True, hbin=True, error = True):

        ctl = sfit4_ctl()
        for tc in self.results.keys():
            print('Entering %s'%(self.results[tc]['dir']))            
            os.chdir(os.path.join(self.testcase_dir, self.results[tc]['dir']))
            flag = False
            if hbin and self.results[tc]['hbin']:
                print('Calling hbin')
                chbin = os.path.join(self.sfit4_dir,'src','hbin')
                rhbin = subprocess.Popen(chbin,stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
                if len(rhbin[1]):
                    print 'hbin failed in %s'%(tc)

                for l in rhbin[0].split('\n'):
                    if l.find('Saving binary')>-1:
                        flag = True
                        fhbin = l.split(':')[1].strip()
                        print('hbin file in %s: %s')%(tc,fhbin)
                if not flag:
                    print('Something wrong in %s'%tc)
                    print('hbin terminated with:')
                    print rhbin[0]
                    self.results[tc].update({hbin:False})
                    
                ctl.replace_in_file('sfit4.ctl','file.in.linelist',fhbin)

            if sfit4 and self.results[tc]['sfit']:
                if error:
                    ctl.replace_in_file('sfit4.ctl','kb','T')
                else:
                    ctl.replace_in_file('sfit4.ctl','kb','F')
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
                    print rsfit[1]
                    self.results[tc].update({sfit4:False})
                            #
#            sf4l0args = ['-i',tcpath,
#                         '-b',os.path.join(self.sfit4_dir,'src'),
#                         '-f',flag]
#            sfit4Layer0(sf4l0args)


    def read_summaries(self):
        
        # Store values from summary
        for tc in self.results:

            sum_file = os.path.join(self.testcase_dir,
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
                          'chi_y_2':sum_new.chi_y_2
                }
                self.results[tc].update(result)
            else:
                self.results[tc].update({'summary':False})

    def read_statevectors(self):
        #Store values from statevec
        state_orig = statevec(os.path.join(orig_testcases,'statevec.%s'%(tc[1])))
        result = {'ret_profile': state_orig.rt_vmr[0]}
        results_orig[tc[0]].update(result)
        
        state = statevec(os.path.join(tcpath,'statevec'))
        result = {'ret_profile': state.rt_vmr[0]}
        results[tc[0]].update(result)


    def print_results(self):
        
        print '\t\t\tThis run\t', 'Orinignal run\t', 'Difference in Percent'

        for rs in self.results.keys():
            print 'Testcase for:', rs
            if self.results[rs]['summary'] and self.results_orig[rs]['summary']:
                print 'Target Apriori:\t\t', self.results[rs]['apriori'], '\t', self.results_orig[rs]['apriori'], '\t', 200*(self.results[rs]['apriori']-self.results_orig[rs]['apriori'])/(self.results[rs]['apriori']+self.results_orig[rs]['apriori']), '%'
                print 'Target Retrieved:\t', self.results[rs]['retriev'], '\t', self.results_orig[rs]['retriev'], '\t', 200*(self.results[rs]['retriev']-self.results_orig[rs]['retriev'])/(self.results[rs]['retriev']+self.results_orig[rs]['retriev']), '%'
                print 'CHI_Y_2:\t\t', self.results[rs]['chi_y_2'], '\t', self.results_orig[rs]['chi_y_2'], '\t', 200*(self.results[rs]['chi_y_2']-self.results_orig[rs]['chi_y_2'])/(self.results[rs]['chi_y_2']+self.results_orig[rs]['chi_y_2']), '%'
                #            print 'MEAN SQARE Diff. RETRIEVED VMR:', np.sqrt(np.mean((self.results[rs]['chi_y_2']-self.results_orig[rs]['chi_y_2'])**2))
            else:
                print 'No new summary or original summary file found'

if __name__ == '__main__':

    tc = test_sfit4('test.cfg')
    runsfit = True
    runhbin = True
    error = True
    if sys.argv.count('--nosfit4') > 0:
        runsfit = False
    if sys.argv.count('--nohbin') > 0:
        runhbin = False
    if sys.argv.count('--noerror') > 0:
        error = False
    tc.run_sfit4_in_testcase(sfit4=runsfit,hbin=runhbin,error=error)
    tc.read_summaries()
#    tc.read_statevectors()
    tc.print_results()
