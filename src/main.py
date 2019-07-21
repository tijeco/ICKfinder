# import tensorflow as tf
import numpy as np
import pandas as pd
from random import shuffle
import random
import os
import sys
from itertools import groupby
import argparse
import json

import src
import src.create as cr
import src.classify as cl
# import src/create as cr
# import toxify.fifteenmer as fm
# import toxify.protfactor as pf
# import toxify.seq2window as sw
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
def fa2dict(fa):
    fasta_dict = {}
    with open(fa) as f:
        header = ""
        for line in f:
            if line.strip() != '':
                if line[0] == ">":
                    header = line.strip().strip(">").split()[0]
                    fasta_dict[header] = ""
                else:
                    fasta_dict[header] += line.strip()
    return fasta_dict

class ParseCommands(object):

    def __init__(self):

        parser = argparse.ArgumentParser(
            description='Pretends to be git',
            usage='''cluck <command> [<args>]
The most commonly used cluck commands are:

    describe
    create
    classify
    train
    predict''')

        parser.add_argument("command", help="Subcommand to run")
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)
        self.args = parser.parse_args(sys.argv[1:2])
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def describe(self):
        parser = argparse.ArgumentParser(
            description="Provides description of cysteine motifs in dataset")
        args = parser.parse_args(sys.argv[2:])

        print("Running cluck describe")
        self.args = args
        return(self.args)

    def classify(self):
        parser = argparse.ArgumentParser(
            description="Provides annotation of cysteine motifs in dataset")
        parser.add_argument('fasta')
        parser.add_argument("-motif", type=str2bool, nargs='?',
                        const=True, default=True,
                        help="Activate nice mode.")
        parser.add_argument('-out',type = str,default = "cluck.csv")
        args = parser.parse_args(sys.argv[2:])

        self.args = args
        print("Running cluck classify",self.args.fasta)
        fasta_file = self.args.fasta
        fasta_dict = fa2dict(fasta_file)
        cluck_json = os.path.abspath(src.__file__).replace("__init__.py","cluck.json")

        fasta_annotate = cl.annotateFastaDict(fasta_dict,cluck_json,motif=self.args.motif)
        fasta_annotate.to_csv(self.args.out,index=False)
        print(self.args.motif)
        # print(os.path.abspath(src.__file__).replace("__init__.py","models/max_len_50"))

        return(self.args)

    def create(self):
        parser = argparse.ArgumentParser(
            description="Creates cysteine motif json file")
        parser.add_argument('fasta')
        parser.add_argument('-json',type = str,default = "cluck.json")
        args = parser.parse_args(sys.argv[2:])

        self.args = args
        print("Running cluck create",self.args.fasta)
        fasta_file = self.args.fasta
        output_json = self.args.json

        faDict = fa2dict(fasta_file)
        faICK = cr.faDict2json(faDict)
        with open(output_json, 'w') as out:
            json.dump(faICK, out)

        return(self.args)
def main():
    cluck_args = ParseCommands().args
    print(cluck_args)
