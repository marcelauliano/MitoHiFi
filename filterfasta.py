#! /usr/bin/env python

'''
    Copyright 2014-2019 Felix Heeger
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    '''



import sys, re

try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False
from random import sample
from optparse import OptionParser

from Bio import SeqIO

def filterFasta(inStream, outPath, minLength=None, idList=None, 
                random=None, fastq=False, regex=False, neg=False, 
                log=sys.stderr):
    if fastq:
        format = "fastq"
    else:
        format = "fasta"
    if random:
        sampleRandom(inStream, outPath, format, random, log)
    else:
        filterLengthIdList(inStream, outPath, format, minLength, idList, 
                           regex, neg, log)
    

def filterLengthIdList(inStream, outPath, format, minLength=None, 
                       idList=None, regex=False, neg=False, log=sys.stderr):
    log.write("".join(("isNeg? ", str(neg))))
    if not idList is None:
        if regex:
            idRes = [re.compile(x) for x in idList]
        else:
            #use a dict to have random acces in O(1)
            idDict = dict(zip(idList, [None]*len(idList)))
    try:
        out = open(outPath, "w")
    except TypeError: 
        #if the outPath paramerter is not a path but already a stream
        out = outPath
    if log:
        log.write("Filtering and writing...\n")
    try:
        l = 0
        n = 0
        for rec in SeqIO.parse(inStream, format):
            l+=1
            if log and l%1000==0:
                log.write("\r%i records done" % l)
            write = True
            if not minLength is None and len(rec)<minLength: 
                #if a min length was set skip the record if is to short
                write = False
            elif not idList is None:
                if regex:
                    if all([r.match(rec.id) is None for r in idRes]):
                        #id the record id does not match any of the given REs
                        # i.e. all results of matching are None
                        write = False
                elif rec.id not in idDict:
                    #if a ID list was given skip the record if it is not in it
                    write = False
            if (not neg and write) or (neg and not write) :
                out.write(rec.format(format))
                n+=1
        if log:
            log.write("\n%i out of %i sequences remained after filtering.\n"
                      % (n,l))
    finally:
        out.close()   

def sampleRandom(inStream, outPath, format, number, log):
    try:
        try:
            out = open(outPath, "w")
        except TypeError:
            out = outPath
        if log:
            log.write("Counting sequences in file...\n")
        l = 0
        for rec in SeqIO.parse(inStream, format):
            l+=1
        inStream.seek(0, 0)
        if number>l:
            raise ValueError("Sample size(%i) is bigger than number of "
                             "sequences in input file(%i)." % (number, l))
        if log:
            log.write("Sampling %i from %i sequences...\n" % (number, l))
        n = 0 #record number
        i = 0 #index of next sampled record
        ls = sample(range(0,l), number)
        ls.sort()
        for rec in SeqIO.parse(inStream, format):
            if n == ls[i]: #if this record is the nextsampled
                out.write(rec.format(format))
                i+=1
                if i>=len(ls):
                    #stop iteration over the input file if all
                    # sampled records were written
                    break
            n+=1    
    finally:
        out.close()   

if __name__ == "__main__":

    usage = "usage: %prog [options] input.fasta [output.fasta]"

    parser = OptionParser(usage)
    
    parser.add_option("-q", "--quite",
                       action="store_true", dest="quiet", default=False, 
                       help="do not print status messages to the screen",)
    parser.add_option("-u", "--fastq",
                       action="store_true", dest="fastq",
                       default=False, help="input file is fastq",)
    if gzImported:
        parser.add_option("-z", "--gzip",
                           action="store_true", dest="gzip",
                           default=False, help="input file is gzipped",)
    parser.add_option("-l", "--min-length",
                       action="store", type="int", dest="minLength",
                       default=None, 
                       help="write only sequence with lengths at least X",
                       metavar="X")
    parser.add_option("-i", "--id-list",
                       action="store", type="string", dest="idList",
                       default=None, 
                       help="write only sequence with an ID from this list. "
                            "List can be comma separated string of IDs or a "
                            "path to a file with a line separated list of IDs",
                       metavar="X")
    parser.add_option("-r", "--random",
                      action="store", type="int", dest="random",
                      default=None, 
                      help="randomly sample X sequence from input file",
                      metavar="X")
    parser.add_option("-e", "--regexp",
                       action="store_true", dest="regexp",
                       default=False, 
                       help="use regular expression instead of exact "
                            "matching for IDs",)
      
    parser.add_option("-n", "--negative",
                       action="store_true", dest="neg",
                       default=False, 
                       help="do exactly the opposite of what would normally "
                            "be done",)
    (options, args) = parser.parse_args()
    
    if (options.idList and options.random):
        parser.error("Options -i and -r are mutually exclusive.")
    if (options.minLength and options.random):
        parser.error("Options -l and -r are mutually exclusive.")
    if (options.regexp and not options.idList):
        parser.error("Options -e can only be used with -i.")
    if (options.random and options.neg):
        parser.error("Negative mode does not work with random mode.")
    
    if options.quiet:
        log = None
    else:
        log = sys.stderr
        
    if len(args) < 1:
        if gzImported and options.gzip:
            parser.error("Pipe mode (no input file argument) does not work together with -z (gzipped input).")
        log.write("NOTE: Running in pipe mode. Witing for input from stdin.\n")
        log.write("Will be writing to stdout.\n")
        out = sys.stdout
    elif len(args) == 1:
        #if no output file was given write to std out
        log.write("Will be writing to stdout.\n")
        out = sys.stdout
    else:
        out = args[1]
    
    if len(args)>2:
        if log:
            log.write("Additional arguments will be ignored!\n")
    if options.neg:
        if log:
            log.write("NOTE: Running in negative mode.\n")
    idList = None
    if options.idList:
        idList = []
        try:
            iListFile = open(options.idList)
            for line in iListFile:
                idList.append(line.strip())
            iListFile.close()
        except IOError:
            pass
            if log:
                log.write("ID list parameter is not a valid path. Assume it to "
                          "be comma separated string.\n")
            idList = options.idList.strip().split(",")
    if log:
        if options.regexp:
            log.write("Using list of Regular Expression on IDs to filter.\n")
        elif idList:
            log.write("Using list of IDs to filter.\n")
                    
    if len(args) == 0:
        inStream = sys.stdin
    elif gzImported and options.gzip:
        inStream = gzip.open(args[0], "r")
    else:
        inStream = open(args[0], "r")
    
    try:
        filterFasta(inStream, out, options.minLength, idList, options.random, 
                    options.fastq, options.regexp, options.neg, log=log)
    finally:
        if len(args) > 0:
            inStream.close()
