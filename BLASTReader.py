import xml.etree.ElementTree as ET
def blast_read(filename, mode = 'pdbaa'):
    if mode in _modelist :
        reader = _modelist[mode] #
    else:
        raise Exception('unexpect reader mode {}'.format(mode))
    return reader(filename)
class PDBAA(object):
    def __init__ (self, filename):
        self.tree = ET.parse(filename)
        self.root = self.tree.getroot()
    def Hits(self):#Hit accession
        res = self.root.findall('.//Hit_accession')
        res = map(lambda x: x.text, res)
        res = map(lambda x: x.split('_'), res)
        res = map(lambda x: {'accession':x[0],'chain':x[1]}, res)
        return res
    def Identicals(self, thredshold):
        hits = self.root.findall('.//Hit')
        ret = []
        for h in hits:
            identity = h.find('.//Hsp_identity').text
            alignment = h.find('.//Hsp_align-len').text
            if thredshold > float(identity) / float(alignment):
                continue
            hit = h.find('.//Hit_accession').text.split('_')
            ret.append({'accession':hit[0],'chain':hit[1]})
        return ret
class BLASTReader(PDBAA):
    pass
class REFSEQ(BLASTReader):
    def Matches(self, thredshold):
        hits = self.root.findall('.//Hit')
        query_len = int(self.root.find('.//BlastOutput_query-len').text)
        hits = filter(lambda k:
            int(k.find('.//Hsp_hit-from').text) == 1 and\
            int(k.find('.//Hsp_hit-to').text) == query_len,
            hits
        )
        ret = []
        for h in hits:
            identity = h.find('.//Hsp_identity').text
            alignment = h.find('.//Hsp_align-len').text
            if thredshold > float(identity) / float(alignment):
                continue
            ret.append({'accession': h.find('.//Hit_accession').text})
        return ret
_modelist = {
    'pdbaa': PDBAA,
    'refseq_protein': REFSEQ
}
