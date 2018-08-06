import BLASTReader as br
import PPIReader 

import csv, os, re, copy
import pandas as pd
class ppir(PPIReader.PPIReader):
    def leave_only(self, inlist):
        newp = ppir() 
        for key in self.piv.keys():
            if key not in inlist:
                continue
            newp.piv[key] = self.piv[key]
            newp.piv[key] = list(set(newp.piv[key]) & set(inlist))
        return newp
    def push(self, major, partner):
        self.piv[major] = list(set(self.piv.get(major, []) + [partner]))
class PPIStratify(object):
    BLASTs = 'blasts'
    Interactions = 'interactions.csv'
    def __init__(self, BLASTs = None, Interactions = None):
        if BLASTs:
            self.BLASTs = BLASTs
        if Interactions:
            self.Interactions = Interactions
        #build ppi object and check path to blasts is a dir
        self.ppi = ppir()
        self.ppi.read(self.Interactions)
        assert(os.path.isdir(self.BLASTs))
        self.run()
    def run(self):
        self.recordmap = {}
        for p in self.ppi.participant():
            try:
                self.test_and_set_fasta(p)
            except Exception as e:
                self.recordmap[p] = []
                print str(e)
                print 'skip search %s' %p
        print 'evaluation complete'
    def Tier3(self):
        return copy.deepcopy(self.ppi)
    def Tier2(self):
        res = self.recordmap.items()
        res = filter(lambda pair: len(pair[1]), res)
        res = map(lambda k: k[0], res)
        return self.ppi.leave_only(res)
    def Tier1(self): 
        newp = ppir()
        for mp in self.ppi.majors():
            partners = filter(lambda k:
                self.detectInteraction(mp,k),self.ppi[mp])
            if len(partners) == 0:
                continue
            newp.update({mp: partners})
        return newp
                
    def discoverPPI(self):
        #reset hit before startint to search
        self.recordmap = {}
        for mp in self.ppi.majors():
            for p in self.ppi[mp]:
                try:
                    self.detectInteraction(mp, p)
                except KeyboardInterrupt as kie:
                    print str(kie)
                    print 'user stop program when running %s, %s' %(mp, p)
                    raise
                except Exception as e:
                    print str(e)
                    print 'skip search %s v.s. %s' %(mp, p)
        print 'done'
    def log(self, msg, destfile = 'log'):
        with open(destfile, 'a') as f:
            f.write(msg)
            f.write('\n')
        
    def detectInteraction(self, fas1, fas2):
        self.test_and_set_fasta(fas1)
        self.test_and_set_fasta(fas2)
        hit1 = pd.DataFrame(self.recordmap[fas1])
        hit2 = pd.DataFrame(self.recordmap[fas2])
        if hit1.shape[0] * hit2.shape[0] == 0:
            return False
        hit12 = pd.merge(hit1, hit2, on='accession')
        hit12 = hit12.loc[hit12['chain_x'] != hit12['chain_y'],:]
        if hit12.shape[0] > 0:
            print '%s v.s. %s' %(fas1, fas2)
            print '\t' + '\n\t'.join(hit12['accession'])
            map(lambda x: self.log(','.join([fas1, fas2, x])), hit12)
            return True
        else:
            return False
    def test_and_set_fasta(self, fas):
        if fas in self.recordmap.keys():
            return
        path = os.path.join(
            self.BLASTs,
            '{}.xml'.format(fas))
        blast = br.BLASTReader(path)
        self.recordmap[fas] = blast.Identicals(0.9)
