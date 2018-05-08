import pprint
from StringIO import StringIO
import csv
def pivotlike(p):
    return all(
        all(isinstance(x, str) for x in v) for v in p.values()
        )
class pivot(object):
    def __init__(self):
        self.piv = {}
    def read(self, infile):
        with open(infile) as f:
            s = f.read()
            s = s.split('\n')
            s = filter(
                lambda x: x if not x.startswith('#') else None,
                s
            )
            s = '\n'.join(s)
            rows = list(csv.reader(StringIO(s), delimiter=','))
            self.piv={}
            for r in rows:
                if len(r) < 2:
                    continue
                self.piv[r[0]] = r[1:]
            
    def flattened(self):
        ret = [] 
        for k, proteins in self.piv.items():
            for p in proteins:
                ret.append([k, p])
        return ret
    def uniq(self):
        for k in self.piv.keys():
            self.piv[k] = list(set(self.piv[k]))
    def write(self, outfile):
        with open(outfile,'w') as f:
            for k,v in self.piv.items():
                f.write(','.join([k]+v))
                f.write('\n')
        return 
    def update(self, dictlike):
        assert(pivotlike(dictlike))
        self.piv.update(dictlike)
    def push(self, key, value):
        self.update({key: self.piv.get(key, []) + [value]})
    def __str__(self):
        pp = pprint.PrettyPrinter(indent=4)
        return pp.pformat(self.piv) 
    def sort(self):
        for k in self.piv.keys():
            self.piv[k] = sorted(self.piv[k])

if __name__ == '__main__':
    p = pivot()
    p.read('interactions.csv')
    print p
    p.update({'limin':'buzy'})
    print p
