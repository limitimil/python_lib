import xml.etree.ElementTree as ET
class BLASTReader(object):
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

                
