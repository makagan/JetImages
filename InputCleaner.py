#!/usr/bin/env python
'''Library to ingest and clean lines from file/stdin;
they are evaled into python objects and returned
'''

import time
__author__ = 'Josh Cogan'

def smart_open(fname, flag='r'):
    '''Handle normal and gzipped files seamlessly'''
    if fname[-6:] == 'tar.gz' or fname[-4:] == '.tgz':
        import gzip
        return gzip.open(fname, flag)
    return open(fname, flag)


class InputCleanerGen:
    def __init__(self, filenames=[], lastResort=None, wait=False, search=None, repl=None):
        self.listOfLines = []
        self.filenames   = filenames
        self.lastResort  = lastResort
        self.wait        = wait
        self.search      = search

        if repl is None:
            self.repl = r'\1'
        else:
            self.repl      = repl
        self._open_file_ = None

        old = time.time()
        self.nLines = sum(1 for ff in filenames for line in smart_open(ff))
        #print 'It took %d seconds to count the %d lines' % (time.time() - old, self.nLines)

    def size(self):
        '''Max possible, but if search is not none, actual lines will be
        less'''
        return self.nLines

    def getGen(self, maxReturned=-1):
        count = 0

        if isinstance(maxReturned, float): maxReturned = int(maxReturned)
        if not isinstance(maxReturned, int): maxReturned = -1


        sandbox = {}
        for filename in self.filenames:
            self._open_file_ = smart_open(filename, 'r')
            for line in self._open_file_:
                safe = InputCleanerGen.madeSafe(line, self.search, self.repl)
                if safe is None:
                    continue

                if maxReturned >= 0 and count >= maxReturned:
                    break

                count += 1
                yield InputCleanerGen.typify(eval(safe, sandbox))


            self._open_file_.close()
            self._open_file_ = None

            if maxReturned >= 0 and count >= maxReturned:
                return

    def __del__(self):
        if self._open_file_ is not None:
            self._open_file_.close()


    @staticmethod
    def madeSafe(dirty, search, repl):
        '''Output from this function is "safe" for eval'''
        import sys
        import re

        if search is not None and None is re.search(search, dirty):
            return None

        if search is not None:
            dirty = re.sub(search, repl, dirty)
        words = set(re.findall('(\w+)', dirty))

        if len( words & set(('import', 'os', 'sys'))) > 0:
            print 'WOW DANGEROUS INPUT\n%s\nGoodbye' % dirty
            sys.exit(1)

        dirty = dirty.strip()
        return dirty if len(dirty) > 0 else None

    @staticmethod
    def typify(di):
        typed = {}
        for k,v in di.iteritems():
            try:
                val = float(v)
                if val == 0:
                    val = 0
                elif (val - int(val))/val < 1e-6:
                    val = int(val)
            except ValueError:
                val = v
            except TypeError:
                val = v
    
            typed[k] = val
    
        return typed

class InputCleaner:
    cls_search = r'^(.*)$'
    cls_repl   = r'\1'

    def __init__(self, filenames=[], lastResort=None, wait=False, search=cls_search, repl=cls_repl):
        self.listOfLines = []
        self.filenames   = filenames
        self.lastResort  = lastResort
        self.wait        = wait
        self.search      = search
        self.repl        = repl

        if not self.wait:
            self.go()

    def go(self):
        for filename in self.filenames:
            fh = smart_open(filename, 'r')
            self._ingest(fh)
            fh.close()

        if len(self.filenames) == 0:
            print 'Falling back on stdin!'
            self._ingest(self.lastResort) # assumes its sys.stdin

    def add(self, fHandle):
        self.inputs.append(fHandle)

    def _ingest(self, fHandle):#, mapAndFilter=None):
        '''Here we extract all the science from a unix file'''
        sandbox = {}
        safeLines = [InputCleaner.madeSafe(line, self.search, self.repl)
                        for line in fHandle]
        listOfDicts = [InputCleaner.typify(eval(line, sandbox))
                        for line in safeLines if line is not None]
        self.listOfLines.extend([di for di in listOfDicts])

    def clear(self):
        self.inputs      = None
        self.lastResort  = None
        self.flush()

    def flush(self):
        self.listOfLines = []
    
    def excreteList(self):
        ret = list(self.listOfLines)
        return ret

    def excreteList(self):
        ret = list(self.listOfLines)
        return ret


    @staticmethod
    def madeSafe(dirty, search='#.*', repl=''):
        '''Output from this function is "safe" for eval'''
        import sys
        import re

        if search is not None and None is re.search(search, dirty):
            return None

        if search is not None:
            dirty = re.sub(search, repl, dirty)
        words = set(re.findall('(\w+)', dirty))

        if len( words & set(('import', 'os', 'sys'))) > 0:
            print 'WOW DANGEROUS INPUT\n%s\nGoodbye' % dirty
            sys.exit(1)
    
        dirty = dirty.strip()
        return dirty if len(dirty) > 0 else None
            
    @staticmethod
    def typify(di):
        typed = {}
        for k,v in di.iteritems():
            try:
                val = float(v)
                if val == 0:
                    val = 0
                elif (val - int(val))/val < 1e-6:
                    val = int(val)
            except ValueError:
                val = v
            except TypeError:
                val = v
    
            typed[k] = val
    
        return typed

if __name__ == '__main__':

    import sys
    ic = InputCleaner(sys.argv[1:], sys.stdin)

    print ic.excreteList()
