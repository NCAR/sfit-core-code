# scripts copied from sfit-processing-environment/Lib_MP in Dev_Mathias

import string,re
import numpy as np

class sfit4_ctl:

    def __init(self):
        pass


    def read_ctl_file(self, ctlfile):

        fid = open(ctlfile, 'r')
        
        self.value = {}
        old_tag = ''
        for line in fid:
            line = line.strip()
#            pdb.set_trace()
            if line.find('#')>-1:
                line = line[0:line.find('#')].strip()
            if len(line)==0:
                if len(old_tag) > 0:
                    self.value[old_tag] += ' ' + line.strip().strip()
                continue;
            if line.find('=') == -1:
                self.value[old_tag] += ' ' + line.strip().strip()
                continue
            tags = line.split('=')
            old_tag = tags[0].strip()
            if len(tags) == 2:
                self.value[tags[0].strip()] = tags[1].strip()


        fid.close()

    def get_value(self, tag):

        tag = tag.strip()
        if self.value.has_key(tag):
            return self.value[tag].strip()
        else:
            print 'key ' + tag.strip() + ' not found'
            return(-1)

    def get_keys(self, level=''):
        # gets all keys under key level
        # e.g. for key file.out.summary 
        # get_keys('file.out') it returns ['summary']
        # if level is not given, it returns all keys
        res = []
        if level.endswith('.'):
            level = level[0:-1]
        else:
            level = level+'.'            
        for v in self.value.keys():
            if level == '.':
                key = v
            else:
                key = v.partition(level)[2]
            if key != '':
                res.append(key)
        return(res)
        
    def replace(self, tag, value):

        tag = tag.strip()
        if self.value.has_key(tag):
            self.value[tag] = value.strip()
        else:
            print 'key ' + tag.strip() + ' not found'

    def replace_in_file(self, ctlfile, newtag, newvalue):
        fido = open(ctlfile, 'r')
        value = ''
        tag = ''
        newlines = []
        flag = False
        for line in fido:
            line = line.strip()
            comment = ''
            if line.find('#')>-1:
                #comment = line[line.find('#'):]
                line = line[0:line.find('#')].strip()
            if len(line)==0:
                if len(value) > 0:
                    newlines.append(tag + ' = ' + value + 
                                    ' ' + comment + '\n')
                    value = ''
                    comment = ''
                elif len(comment)>0:
                    newlines.append(comment + '\n')
                    comment = ''
                continue
            old_line = line
            if len(value)> 0 and line[0].isalpha():
                if (tag == newtag):
 #                   print tag
                    newlines.append(tag + ' = ' + newvalue +
                                ' ' + comment + '\n')
                else:
#                    print tag
                    newlines.append(tag + ' = ' + value + 
                                ' ' + comment + '\n')
#                print tag + ' = '  + value + '\n'
                value = ''
                comment = ''
            if line.find('=') == -1: 
                if tag == newtag:
                    dum = line.strip().strip()
                else:
                    value  += ' ' + line.strip().strip()
                    value = value.strip()
#                    print value 
                continue
            tags = line.split('=')
            tag = tags[0].strip()
            value = tags[1].strip()
            if tag.lower() == newtag.strip().lower():
                value = newvalue
                flag = True

        # last line has not yet been written!
        if len(value) > 0:
            newlines.append(tag + ' = ' + value + 
                            ' ' + comment + '\n')
            

        if not flag:
            newlines.append(newtag + ' = ' + newvalue + '\n')
        fido.close()
        fidn = open(ctlfile, 'w')
        fidn.writelines(newlines)
        fidn.close()


    def write(self, ctlfile):
        fid = open(ctlfile, 'w')
        for k,v in self.value.iteritems():
            fid.write(k + ' = ' + v + '\n')

        fid.close()

class summary:
    def __init__(self, filename):
        sumf = read_from_file(filename)
        self.header = sumf.get_line()
        nr_spc = int(sumf.next(1).pop(0))
        sumf.skipline(nr_spc)

        nr_retgas = int(sumf.next(1).pop(0))
        head1 = sumf.get_line().split()
        self.gas = []
        self.prf = []
        self.cell = []
        self.apriori = np.zeros(nr_retgas)
        self.retriev = np.zeros(nr_retgas)
        for nr in range(0,nr_retgas):
            int(sumf.next(1).pop(0))
            self.gas.extend(sumf.next(1))
            self.prf.extend(sumf.next(1))
            if head1.count('IFCELL') > 0:
                 self.cell.extend(sumf.next(1))
            self.apriori[nr] = sumf.next(1).pop(0)
            self.retriev[nr] = sumf.next(1).pop(0)

        nr_band = int(sumf.next(1).pop(0))
        sumf.skipline(1)
        self.mw_start = np.zeros(nr_band)
        self.mw_stop = np.zeros(nr_band)
        self.spac = np.zeros(nr_band)
        self.points = np.zeros(nr_band)
        self.max_opd = np.zeros(nr_band)
        self.fov = np.zeros(nr_band)
        self.snr_ret = np.zeros(nr_band)
        self.snr_apr = np.zeros(nr_band)

        try:
            for nr in range(0,nr_band):
                sumf.next(1)
                self.mw_start[nr] = sumf.next(1).pop(0)
                self.mw_stop[nr] = sumf.next(1).pop(0)
                self.spac[nr] = sumf.next(1).pop(0)
                self.points[nr] = sumf.next(1).pop(0)
                self.max_opd[nr] = sumf.next(1).pop(0)
                self.fov[nr] = sumf.next(1).pop(0)
                sumf.next(3)
                self.snr_apr[nr] = sumf.next(1).pop(0)
                self.snr_ret[nr] = sumf.next(1).pop(0)
                sumf.skip_reminder_of_line()
        except:
            pass

        sumf.skipline(1)
        try:
            self.rms = sumf.next(1).pop(0)
            self.chi_y_2 = sumf.next(1).pop(0)
            sumf.next(1)
            self.dofs = sumf.next(1).pop(0)
            sumf.next(1)
            self.iter = sumf.next(1).pop(0)
            self.iter_max = sumf.next(1).pop(0)
            self.converged = sumf.next(1)
            self.div_warn = sumf.next(1)
            del sumf
        except:
            pass


class statevec :
    def __init__(self, filename):

        stvf = read_from_file(filename)
        
        stvf.skipline()
        self.nr_layer = stvf.next(1).pop(0)
        stvf.next(1).pop(0)
        stvf.next(1).pop(0)
        self.tretflag = stvf.next(1).pop(0)
        stvf.next(1).pop(0)
        stvf.next(1).pop(0)
        stvf.skipline()
        self.Z = stvf.next(self.nr_layer)
        stvf.skipline()
        self.P = stvf.next(self.nr_layer)
        stvf.skipline()
        self.T = stvf.next(self.nr_layer)
        if self.tretflag == 'T':
            stvf.skipline()
            self.Tret = stvf.next(self.nr_layer)
        
        nr_gas = stvf.next(1).pop(0)
        self.ap_col = []
        self.ap_vmr = []
        self.rt_col = []
        self.rt_vmr = []
        self.gas = []
        for gas in range(0,nr_gas):
            stvf.skipline()
            self.gas.append(stvf.next(1).pop(0))
            self.ap_col.extend(stvf.next(1))
            self.ap_vmr.append(stvf.next(self.nr_layer))
            stvf.skipline()
            stvf.skipline()
            self.rt_col.extend(stvf.next(1))
            self.rt_vmr.append(stvf.next(self.nr_layer))

        nr_aux = stvf.next(1).pop(0)
        self.aux = []
        self.ap_aux = []
        self.rt_aux = []
        for aux in range(0,nr_aux):
            self.aux.append(stvf.next(1).pop(0))
        for aux in range(0,nr_aux):
            self.ap_aux.append(stvf.next(1).pop(0))
        for aux in range(0,nr_aux):
            self.rt_aux.append(stvf.next(1).pop(0))
        del stvf


class read_from_file:
    def __init__(self, filename):
        self.file = open(filename,'r')
        self.line = self.file.readline().split()
        self.count = len(self.line)
        self.flag = 1

    def __del__(self):
        self.file.close()

    def next(self, nr):
        parse = lambda tmp: tmp.isdigit() and int(tmp) \
            or \
            len(set('.,').intersection(tmp)) == 1 and \
            (tmp.count('.') + tmp.count(',')) == 1 and \
            len(set(string.ascii_letters).intersection(tmp)) == 0 and\
            float(tmp) or\
            len(set('.,').intersection(tmp)) <= 1 and \
            (tmp.count('.') + tmp.count(',')) <= 1 and \
            len(set(string.ascii_letters).intersection(tmp)) == 1 and\
            len(set('eEdD').intersection(tmp)) == 1 and \
            (tmp.count('e') + tmp.count('E') + tmp.count('d') + tmp.count('D')) <= 1 and \
            float(tmp) \
            or tmp            
        
        nrs = []
        for n in range(0,nr):
            if self.count < 1 :
                self.nextline()
            if self.flag < 0:
                return()
            tmp = parse(self.line.pop(0))
            # parse function does not work if number is zero.
            if tmp == '0':
                tmp = int(tmp)
            # and checking for floats
            if type(tmp) == str:
                try:
                    tmp = float(tmp)
                except:
                    pass
            nrs.extend([tmp])
            self.count = self.count - 1
        return(nrs)

    def skip_reminder_of_line(self):
        self.count = 0
        self.nextline()

    def get_line(self):
        if self.count < 1:
            self.nextline()
        line = string.join(self.line, ' ')
        while len(line.strip()) == 0:
            self.nextline()
            if self.stat:
                return()
            line = string.join(self.line, ' ')
        self.nextline()
        return(line)

    def nextline(self):
        self.count = 0
        while self.count == 0:

            sep = re.compile('[ ,]')
            self.line = sep.split(self.file.readline().strip())
            self.line = filter(bool,self.line)
            self.count = len(self.line)
            if self.count < 1:
                # Check for EOF if empty line. May be more elegantly solved using read_lines(), 
                # but dont know yet how to do this
                if len(self.file.read(1)) == 0:
                    self.flag = -1
                    break
                self.file.seek(-1,1)
                
    def skipline(self, nr=1):
        for n in range(0,nr):
            self.get_line()
#        self.nextline()

    def stat(self):
        return(self.flag)
        
