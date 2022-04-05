import sys

class BidirectedOverlapGraph(object):
    def __init__(self, alignments, maxOverhang=25):
        # alignments is a vector of 9-tuples:
        #     (seqV, seqH, rc, begV, endV, begH, endH, lenV, lenH)

        self.graph = {}
        self.maxOverhang = maxOverhang

        for seqV,seqH,rc,begV,endV,begH,endH,lenV,lenH in alignments:
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

            if not e['rc']:
                if e['begV'] > e['begH']:
                    e['dir'] = 1
                    e['suffix'] = e['lenH'] - e['endH']
                    e['prefix'] = e['begV']
                else:
                    e['dir'] = 2
                    e['suffix'] = e['begH']
                    e['prefix'] = e['lenV'] - e['endV']
            else:
                if e['begV'] > 0 and e['begH'] > 0 and e['lenV']-e['endV'] == 0 and e['lenH']-e['endH'] == 0:
                    e['dir'] = 0
                    e['suffix'] = e['begH']
                    e['prefix'] = e['begV']
                else:
                    e['dir'] = 3
                    e['suffix'] = e['lenH'] - e['endH']
                    e['prefix'] = e['lenV'] - e['endV']

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

            if not (seqV-1) in self.graph:
                self.graph[seqV-1] = {}

            if (seqH-1) not in self.graph:
                self.graph[seqH-1] = {}

            self.graph[seqV-1][seqH-1] = e
            self.graph[seqH-1][seqV-1] = eT

    def write_gml(self, filename):
        with open(filename, "w") as f:
            f.write("graph\n[\n\tdirected 1\n")
            for i in self.graph:
                f.write("\tnode\n\t[\n\t\tid " + str(i) + "\n\t]\n")
            for seqV in self.graph:
                for seqH in self.graph[seqV]:
                    e = self.graph[seqV][seqH]
                    f.write("")
                    f.write("\tedge\n\t[\n\t\tsource {}\n\t\ttarget {}\n\t\tdir {}\n\t\tsuffix {}\n\t]\n".format(seqV, seqH, e['dir'], e['suffix']))
            f.write("]");


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

    @classmethod
    def direct_to_file(cls, filename, output_prefix, maxOverhang=25):
        G = cls.read(filename, maxOverhang)
        G.tofile(output_prefix + ".txt")
        G.write_gml(output_prefix + ".gml")


if __name__ == "__main__":
    if len(sys.argv) > 3:
        print("running with maxOverhang={}".format(int(sys.argv[3])))
        BidirectedOverlapGraph.direct_to_file(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        BidirectedOverlapGraph.direct_to_file(sys.argv[1], sys.argv[2])
