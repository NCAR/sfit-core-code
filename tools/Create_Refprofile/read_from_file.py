import string, re
import pdb
import matplotlib.pyplot as plt

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
        line = ' '.join(self.line)
        while len(line.strip()) == 0:
            self.nextline()
            if self.stat:
                return()
            line = ' '.join(self.line)
        self.nextline()
        return(line)

    def nextline(self):
        self.count = 0
        while self.count == 0:

            sep = re.compile('[ ,]')
            self.line = sep.split(self.file.readline().strip())
            self.line = list(filter(bool,self.line))
            self.count = len(self.line)
            if self.count < 1:
                # Check for EOF if empty line. May be more elegantly solved using read_lines(), 
                # but dont know yet how to do this
                if len(self.file.read(1)) == 0:
                    self.flag = -1
                    break
                self.file.seek(self.file.tell()-1,0)
                
    def skipline(self, nr=1):
        for n in range(0,nr):
            self.get_line()
#        self.nextline()

    def stat(self):
        return(self.flag)
