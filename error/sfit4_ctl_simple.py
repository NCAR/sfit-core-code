# Written by Matthias Palm, distrtibuted at Error Analysis Workshop in Boulder 2013
# replace_in_file function has been removed in this 'simple' version.
# Updated by Stephanie Conway for the sfit4 workshop in Tsukuba, June 2013.

# This program opens the sfit4.ctl file and sorts them into an object consisting of 
# tag, value pairs (read_ctrl_file).  Also includes functions for changing values 
# in sfit4.ctl object (replace), writting the sfit.ctl object to file (write).  It
# should be noted that after reading the sfit4.ctl file in, writing it removes all 
# comments. 

class sfit4_ctl:

    def __init(self):
        pass

    def read_ctl_file(self, ctlfile):
    # This function opens the sfit4.ctl file, and sorts the entries into tag, value 
    # pairs into an object called value 

    # Open sfit4.ctl file
      with open(ctlfile) as fopen:
        self.value = {}
        old_tag = ''
    # Iterate through lines in sfit4.ctl file
        for line in fopen:
          line = line.strip()
    # Remove all comments (anything after a '#' in the file)
          if line.find('#')>-1:
            line = line[0:line.find('#')].strip()
    # Handle lines where the value may be on multiple lines (i.e. apriori values)
          if len(line)==0:
            if len(old_tag) > 0:
              self.value[old_tag] += ' ' + line.strip().strip()
            continue;
          if line.find('=') == -1:
            self.value[old_tag] += ' ' + line.strip().strip()
            continue
     # Split the entry at the '=' sign
          tags = line.split('=')
          old_tag = tags[0].strip()
          if len(tags) == 2:
            self.value[tags[0].strip()] = tags[1].strip()

    def replace(self, tag, value):
    # This function replaces the value associated with the tag  
      tag = tag.strip()
      if self.value.has_key(tag):
        self.value[tag] = value.strip()
      else:
        print 'key ' + tag.strip() + ' not found'

    def write(self, ctlfile):
    # This function writes the tag, value pairs seperated by an '=' sign, 
    # effectively overwriting the old sfit4.ctl file
      fid = open(ctlfile, 'w')
      for k,v in self.value.iteritems():
        fid.write(k + ' = ' + v + '\n')

      fid.close()


