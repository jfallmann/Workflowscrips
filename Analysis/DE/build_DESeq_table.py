import sys, argparse
import gzip

#import os
parser = argparse.ArgumentParser(description='build DESeq2 table')
parser.add_argument("-l", "--list", dest='my_list', required=True, type=argparse.FileType('r'), help="toDo")
parser.add_argument("-n", "--sample_name", action="store_true", help=" provide -n if sample names instead of group names should be used for header" )
parser.add_argument("-o", "--order", dest='order', required=False, default='NA',help="if wanted the order of conditions can be given as comma separated list" )
parser.add_argument("-c", "--cutoff", dest='cutoff', default=int(0) ,help="cutoff for minimum count" )

args = parser.parse_args()

my_groups={}

class Sample_list(object):

    group_name=""
    replicate_names=[]
    replicate_paths=[]
    # the class constructor
    def __init__(self, group_name):
        self.group_name = group_name
        self.replicate_names = []
        self.replicate_paths = []


def make_sample_list(group_name):
    sample_list=Sample_list(group_name)
    return sample_list

list_size=0
for line in args.my_list:
    list_size+=1
    columns = line.strip().split(' ')
    if columns[0] in my_groups:
        my_groups[columns[0]].replicate_paths.append(columns[2])
        my_groups[columns[0]].replicate_names.append(columns[1])
    else:
        my_groups[columns[0]]=make_sample_list(columns[0])
        my_groups[columns[0]].replicate_names.append(columns[1])
        my_groups[columns[0]].replicate_paths.append(columns[2])

myMatrix = []
myMatrix.append([])
myMatrix[0].append("names")
sample_counter=0

conds = []
if not 'NA' in args.order:
   conds = args.order.split(',')
else:
    conds = [x for x in my_groups]

for gruppies in conds:

    condition_index=-1

    for replicates in my_groups[gruppies].replicate_paths:
        print >> sys.stderr,'Processing '+str(replicates)
        condition_index +=1
        sample_counter+=1
        if (args.sample_name):
            myMatrix[0].append(my_groups[gruppies].replicate_names[condition_index])
        else:
            myMatrix[0].append(my_groups[gruppies].group_name)
        if '.gz' in replicates:
            myInput = gzip.open(replicates,'r')
        else:
            myInput = open(replicates,'r')
        lineNumber=0
        for line in myInput:
            columns = line.strip().split('\t')
            if columns[0] != "name" and columns[1]!="count":
                lineNumber+=1
                if sample_counter==1:
                    newListi=[]
                    myMatrix.append(newListi)
                    myMatrix[lineNumber].append("l_"+str(lineNumber)+"_"+str(columns[0]))
                myMatrix[lineNumber].append(columns[1])

line = "\t".join(myMatrix[0])
sys.stdout.write(str(line)+"\n")

for z in range(1,len(myMatrix)):
    zeilen = myMatrix[z]
    willprint = False
    line = str(zeilen[0])+"\t"
    for x in range(1,len(zeilen)):
        line = line + str(zeilen[x])+"\t"
        if (int(zeilen[x]) >= int(args.cutoff)):
            willprint = True
    if willprint is True:
        sys.stdout.write(str(line)+"\n")
