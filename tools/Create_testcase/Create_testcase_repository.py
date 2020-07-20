import os, sys
from sfit4_ctl import sfit4_ctl
import tarfile

def create_testcase_rep(dir='.'):

    ctl = sfit4_ctl()
    ctl.read_ctl_file(os.path.join(dir,'sfit4.ctl'))
    
    tarf = tarfile.open('testcase.tgz', 'w:gz')

    
    if os.path.exists(os.path.join(dir,'sfit4.ctl')):
        tarf.add(os.path.join(dir,'sfit4.ctl'))
    else:
        print ('dir {} does not contain an sfit4.ctl'.format(os.path.abspath(dir)))
        tarf.close()
        return()
    if os.path.exists(os.path.join(dir,'hbin.ctl')):
        tarf.add(os.path.join(dir,'hbin.ctl'))

    # necessary files
    nec_files = {'file.in.spectrum':'always',
                 'file.in.refprofile':'always',
                 'file.in.stalayers':'always',
                 'file.in.solarlines':'fw.solar_spectrum',
                 'file.in.isotope':'fw.isotope_separation'}
    
    for necf in nec_files:
        if nec_files[necf] == 'always' or ctl.get_value(nec_files[necf]) == 'T':
            necfile = os.path.join(dir, ctl.get_value(necf))
            print(necfile)
            if os.path.exists(necfile):
                tarf.add(necfile)
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
                    tarf.add(necfile)
                else:
                    print ('file  {} does not exists'.format(necfile))
                    tarf.close()
                    return()
                
            
    tarf.close()
               
if __name__ == '__main__':
    create_testcase_rep(sys.argv[1])
               
