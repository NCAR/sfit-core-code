#!/usr/bin/python3

import os, sys
from sfit4_ctl import sfit4_ctl
import tarfile

def create_testcase_rep(dir='.'):

    if not os.path.exists(os.path.join(dir,'sfit4.ctl')):
        print('The directory {} is not a sfit4 directory.'.format(dir))
        return

    ctl = sfit4_ctl()
    ctl.read_ctl_file(os.path.join(dir,'sfit4.ctl'))
    
    tarf = tarfile.open('testcase.tgz', 'w:gz')

    
    if os.path.exists(os.path.join(dir,'sfit4.ctl')):
        tarf.add(os.path.join(dir,'sfit4.ctl'), arcname='sfit4.ctl')
    else:
        print ('dir {} does not contain an sfit4.ctl'.format(os.path.abspath(dir)))
        tarf.close()
        return()
    if os.path.exists(os.path.join(dir,'hbin.ctl')):
        tarf.add(os.path.join(dir,'hbin.ctl'), arcname='hbin.ctl')

    # necessary files
    nec_files = {'file.in.spectrum':'always',
                 'file.in.refprofile':'always',
                 'file.in.stalayers':'always',
                 'file.in.solarlines':'fw.solar_spectrum'}
    
    for necf in nec_files:
        if nec_files[necf] == 'always' or ctl.get_value(nec_files[necf]) == 'T':
            necfile = os.path.join(dir, ctl.get_value(necf))
            if os.path.exists(necfile):
                tarf.add(necfile, arcname=os.path.basename(necfile))
            else:
                print ('file  {} does not exists'.format(necfile))
                tarf.close()
                return()

    # is an SA matrix needed?

    for gas in ctl.get_value('gas.profile.list').split():
        if ctl.get_value('gas.profile.{}.correlation'.format(gas)) == 'T':
            type = ctl.get_value('gas.profile.{}.correlation.type'.format(gas))
            if type == '4' or type == '5':
                necfile = os.path.join(dir, ctl.get_value('file.in.sa_matrix'))
                if os.path.exists(necfile):
                    tarf.add(necfile, arcname=os.path.basename(necfile))
                else:
                    print ('file  {} does not exists'.format(necfile))
                    tarf.close()
                    return()

    # Isotope input needed?

    if ctl.get_value('fw.isotope_separation') == 'T':
        iso_file = os.path.join(dir, ctl.get_value('file.in.isotope'))
        if os.path.exists(iso_file):
            tarf.add(iso_file, arcname=os.path.basename(iso_file))
        else:
            print ('file  {} does not exists'.format(iso_file))
            tarf.close()
            return()
        
    # Apodisation and/or phase function file needed?

    if ctl.get_value('fw.apod_fcn') == 'T':
        apod_file = os.path.join(dir, ctl.get_value('file.in.modulation_fcn'))
        if os.path.exists(apod_file):
            tarf.add(apod_file, arcname=os.path.basename(apod_file))
        else:
            print ('file  {} does not exists'.format(apod_file))
            tarf.close()
            return()
        
    if ctl.get_value('fw.phase_fcn') == 'T':
        phase_file = os.path.join(dir, ctl.get_value('file.in.phase_fcn'))
        if os.path.exists(apod_file):
            tarf.add(phase_file, arcname=os.path.basename(phase_file))
        else:
            print ('file  {} does not exists'.format(phase_file))
            tarf.close()
            return()

# hbin.ctl existing?
        hbin_file = os.path.join(dir, 'hbin.ctl')
        if os.path.exists(hbin_file):
            tarf.add(hbin_file, arcname=os.path.basename(hbin_file))
        else:
            print ('file  {} does not exists'.format(apod_file))
            tarf.close()
            return()
    

    tarf.close()
               
if __name__ == '__main__':

    if len(sys.argv) == 2:
        direc = sys.argv[1]
    else:
        direc = '.'
    direc = os.path.abspath(direc)
    create_testcase_rep(direc)
               
