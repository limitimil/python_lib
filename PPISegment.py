import BLASTReader as br
import PPIReader 

from shutil import copyfile
import csv, os, re, sys
import pandas as pd
from Bio.PDB import PDBList
def help():
    print """
        >>> ppi = PPISegment.PPISegment(BLASTs = 'blasts', Interactions=
        ... 'interactions.csv')
        >>> tab = ppi.make_table()
            .
            .
            .
        >>> PPISegment.transform(tab, pdbrepo = 'pdbs')
        chainA chainB PDBid codeA codeB path locA locB
            .
            .
            .
        >>> 
    """
class BLASTReader(br.BLASTReader):
    def IdenticalRegions(self, threshold):
        hits = self.root.findall('.//Hit')
        ret = []
        for h in hits:
            identity = h.find('.//Hsp_identity').text
            alignment = h.find('.//Hsp_align-len').text
            if threshold > float(identity) / float(alignment):
                continue
            hit = h.find('.//Hit_accession').text.split('_')
            region = (h.find('.//Hsp_query-from').text, 
                h.find('.//Hsp_query-to').text)
            ret.append({'accession':hit[0],'chain':hit[1], 'from': region[0],
                'to': region[1]})
        return ret
class PPISegment(object):
    BLASTs = 'blasts'
    Interactions = 'interactions.csv'
    def __init__(self, BLASTs = None, Interactions = None):
        if BLASTs:
            self.BLASTs = BLASTs
        if Interactions:
            self.Interactions = Interactions
        #build ppi object and check path to blasts is a dir
        self.ppi = PPIReader.PPIReader()
        self.ppi.read(self.Interactions)
        assert(os.path.isdir(self.BLASTs))
    def make_table(self):
        self.recordmap = {}
        ret = pd.DataFrame()
        for mp in self.ppi.majors():
            for p in self.ppi[mp]:
                try:
                    tab = self._detect(mp, p)
                    if tab.shape[0] == 0:
                        continue
                    tab['A'] = mp
                    tab['B'] = p
                    ret = ret.append(tab)
                except KeyboardInterrupt as kie:
                    print str(kie)
                    print 'user stop program when running %s, %s' %(mp, p)
                    raise
                except Exception as e:
                    print str(e)
                    print 'skip search %s v.s. %s' %(mp, p)
        return ret
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
    def _detect(self, fas1, fas2):
        self.test_and_set_fasta(fas1)
        self.test_and_set_fasta(fas2)
        hit1 = pd.DataFrame(self.recordmap[fas1])
        hit2 = pd.DataFrame(self.recordmap[fas2])
        if hit1.shape[0] == 0 or hit2.shape[0] == 0:
            return pd.DataFrame()
        hit12 = pd.merge(hit1, hit2, on='accession', suffixes=('A','B'))
        hit12 = hit12.loc[hit12['chainA'] != hit12['chainB'],:]
        return hit12
        
        
    def detectInteraction(self, fas1, fas2):
        hit12 = self._detect(fas1, fas2)
        if hit12.shape[0] > 0:
            print hit12
            return True
        else:
            return False
    def test_and_set_fasta(self, fas):
        if fas in self.recordmap.keys():
            return
        path = os.path.join(
            self.BLASTs,
            '{}.xml'.format(fas))
        blast = BLASTReader(path)
        self.recordmap[fas] = blast.IdenticalRegions(0.9)
def transform(tab, pdbrepo = 'pdbs'):
    def trans1(p, Id):
        path = os.path.join(
            os.path.dirname(p),
            '{}.pdb'.format(Id)
            )
        try:
            copyfile(p, path)
        except:
            sys.stderr.write('error when downloading structure {}\n'.format(Id))
            return None
        return os.path.relpath(path, os.getcwd())
    def func1(row):
        path = pl.retrieve_pdb_file(row['accession'], file_format='pdb')
        return trans1(path, row['accession'])
    def func2(row):
        ret = [row['fromA'] + '-' + row['toA'],
            row['fromB'] + '-' + row['toB']
            ]
        return pd.Series(ret)
    ret = tab[['A','B','accession','chainA','chainB']]
    #prepare for path/ locA/ locB
    pl = PDBList(pdb = 'pdbs')
    pdbs = tab.apply(func1, 1)
    ret= pd.concat([ret, pdbs], axis= 1)
    locs = tab.apply(func2, 1)
    ret= pd.concat([ret, locs], axis= 1)
    ret.columns = ['chainA', 'chainB', 'PDBid', 'codeA', 'codeB','path',
        'locA', 'locB']
    return ret

class TemplateHandler(object):
    pdbs = 'pdbs'
#    indir = 'pdbs'
#    outdir = 'pdbs2'
#    def __init__(self, indir=None, outdir=None):
#        if indir:
#            self.indir = indir
#        if outdir:
#            self.outdir = outdir
#        assert(os.path.isdir(self.indir))
#        assert(os.path.isdir(self.outdir))
#    def run(self):
#        for folder in os.listdir(self.indir):
#            target = os.path.join(self.indir, folder)
#            for ent in os.listdir(target):
#                self.single(os.path.join(folder, ent))
#        print 'done!'
#        return
#    def single(self, path):
#        fn = os.path.basename(path)[3:7].upper()+'.pdb'
#        shutil.copyfile(
#            os.path.join(self.indir, path),
#            os.path.join(self.outdir, fn))
#        return
    def transfer(self, path, indir, outdir):
        fn = os.path.basename(path)
        copyfile(
            os.path.join(path),
            os.path.join(outdir, fn))
        return 
    def format1(self, tab, pdbrepo):
        return transform(tab, pdbrepo = pdbrepo)
    def format2(self, tab, pdbrepo):
        tab = transform(tab, pdbrepo = self.pdbs)
        for i, row in tab.iterrows():
            if pd.isnull(row['path']):
                continue
            self.transfer(row['path'], self.pdbs, pdbrepo)
            row['path'] = os.path.basename(row['path'])
        return tab
