
class BidirectedOverlapGraph(object):
    def __init__(self, alignments, maxOverhang=25):
        # alignments is a vector of 9-tuples:
        #     (seqV, seqH, rc, begV, endV, begH, endH, lenV, lenH)

        self.graph = {}
        self.maxOverhang = maxOverhang

        for seqV,seqH,rc,begV,endV,begH,endH,lenV,lenH in alignments:
            if not (seqV-1) in self.graph:
                self.graph[seqV-1] = {}
            e = {
              'transpose' :  0,
                    'dir' : -1,
                 'suffix' : -1,
                 'prefix' : -1,
                     'rc' : rc,
                   'begV' : begV, 'endV' : endV,
                   'begH' : begH, 'endH' : endH,
                   'lenV' : lenV, 'lenH' : lenH
                }

            self.graph[seqV-1][seqH-1] = e

            if e['begV'] > 0 and e['begH'] < maxOverhang:
                e['suffix'] = e['lenH'] - e['endH']
                e['prefix'] = e['begV']
                e['dir'] = 1 if not e['rc'] else 3
            elif e['begH'] > 0 and e['begV'] < maxOverhang:
                e['suffix'] = e['begH']
                e['prefix'] = e['lenV'] - e['endV']
                e['dir'] = 2 if not e['rc'] else 0

            eT = {
              'transpose' :  1,
                    'dir' : e['dir'] if e['dir'] in {-1,0,3} else (3 - e['dir']),
                 'suffix' : e['prefix'],
                 'prefix' : e['suffix'],
                     'rc' : rc,
                   'begV' : begV, 'endV' : endV,
                   'begH' : begH, 'endH' : endH,
                   'lenV' : lenV, 'lenH' : lenH
                 }

            if (seqH-1) not in self.graph:
                self.graph[seqH-1] = {}
            self.graph[seqH-1][seqV-1] = eT

    def tofile(self, filename):
        with open(filename, "w") as f:
            f.write("\t".join(["row","col","transpose","dir","suffix","prefix","rc","begV","endV","begH","endH","lenV","lenH"]) + "\n")
            for row in self.graph:
                for col in self.graph[row]:
                    e = self.graph[row][col]
                    f.write(("{}\t"*12 + "{}\n").format(row+1,col+1,e['transpose'], e['dir'], e['suffix'], e['prefix'], e['rc'], e['begV'], e['endV'], e['begH'], e['endH'], e['lenV'], e['lenH']))

    @classmethod
    def read(cls, filename, maxOverhang=25):
        def fileread_generator():
            with open(filename, "r") as f:
                next(f)
                next(f)
                for line in f.readlines():
                    yield [int(v) for v in line.rstrip().split()]
        return cls(fileread_generator(), maxOverhang)


