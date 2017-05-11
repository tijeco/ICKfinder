from itertools import groupby
from itertools import (takewhile,repeat)

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


import sys
def getOptionValue(option):
    optionPos = [i for i, j in enumerate(sys.argv) if j == option][0]
    optionValue = sys.argv[optionPos + 1]
    return optionValue
if "-i" in sys.argv:
    fileName = getOptionValue("-i")
else:
    print "\nplease specify input file name using -i <file_name> \n"
    sys.exit()

sequence_iterator = fasta_iter(fileName)


cysDict = {}
with open("Id.CysMotif.txt","w") as out:
    out.write("ID\tCysMotif\n")
    for ff in sequence_iterator:
        headerStr, seq = ff

        firstC = False
        cysMotif = ""
        number = 0
        for i in seq:

            if i == "C":
                if firstC == False:
                    firstC=True
                    number = 0
                    cysMotif+="C"
                else:
                    cysMotif+="-"+str(number)+"-C"
                    number = 0
            else:
                number+=1
        if cysMotif not in cysDict:
            cysDict[cysMotif] = 1
        else:
            cysDict[cysMotif] +=1
        writeMe = headerStr.split()[0]+"\t"+ cysMotif  +'\n'
        out.write(writeMe)
with open("CysMotifAbundance.txt","w") as out:
    out.write("Abundance\tCysMotif\n")
    for i in cysDict.keys():
        out.write(str(cysDict[i])+"\t"+i+"\n")
