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
        if tag in self.value:
            return self.value[tag].strip()
        else:
            print ('key ' + tag.strip() + ' not found')
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
        if tag in self.value:
            self.value[tag] = value.strip()
        else:
            print ('key ' + tag.strip() + ' not found')

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
