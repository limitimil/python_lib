from StringIO import StringIO
import csv
import pivotHandler as pH
class PPIReader(pH):
    def participant(self):
        major = self.ppi.keys()
        partner = []
        for v in self.ppi.values():
            partner += v
            
        return list(set(major + partner))
if __name__ == '__main__':
    import sys
    pr = PPIReader(open(sys.argv[1]))
    print '\n'.join(map(lambda x: ','.join(x), pr.flattened()))

