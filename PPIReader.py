from StringIO import StringIO
import csv
import pivotHandler as pH
class PPIReader(pH.pivot):
    def participant(self):
        major = self.piv.keys()
        partner = []
        for v in self.piv.values():
            partner += v
            
        return list(set(major + partner))
    def majors(self):
        return self.piv.keys()
    def __getitem__(self, key):
        return self.piv[key]
if __name__ == '__main__':
    import sys
    pr = PPIReader()
    pr.read(sys.argv[1])
    print '\n'.join(map(lambda x: ','.join(x), pr.flattened()))

