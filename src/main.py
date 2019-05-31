# import tensorflow as tf
import numpy as np
import pandas as pd
from random import shuffle
import random
import os
import sys
from itertools import groupby
import argparse
# import toxify
# import toxify.fifteenmer as fm
# import toxify.protfactor as pf
# import toxify.seq2window as sw

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

    def create(self):
        parser = argparse.ArgumentParser(
            description="Creates cysteine motif json file")
        parser.add_argument('fasta')
        args = parser.parse_args(sys.argv[2:])
        print("Running cluck create")
        self.args = args
        return(self.args)
def main():
    cluck_args = ParseCommands().args
    print(cluck_args)
