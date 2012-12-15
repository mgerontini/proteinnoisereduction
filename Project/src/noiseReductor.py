#! /usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio import Alphabet
from Bio import Entrez, SeqIO  
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import collections
import os
import glob
import sys
import time
import itertools
from os.path import basename



'''
Method which calculates the percentage of the
unique characters in a column
'''
def uniqueCharsPercentage(column):
    column_string = ''.join(column)
    my_dict = collections.Counter(column_string)
    counter = 0
    for value in my_dict.itervalues():
       # print value
        if value == 1 or value == 2:
            counter+=1
    perc = counter * (100/len(column_string))
    return perc


def uniqueOccurancePercentage(column):
    column_string = ''.join(column)
    my_dict = collections.Counter(column_string)
    counter = 0
    for value in my_dict.itervalues():
       # print value
        if value == 1 or value == 2:
            counter+=1
    perc = counter * (100/len(column_string))
    return perc



'''
Method which calculates the percentage of the
indels(insertions & deletions) in a column
Type of indels( Deletion, Insertion) 
'''
def indelsPercentage(column):
    column_string = ''.join(column)
    my_dict = collections.Counter(column_string)
    count_del = 0
    count_let = 0
    #print my_dict
    for key,value in my_dict.iteritems():
        if key == '-':
            count_del = value
      
    #check what type of indel is
        perc = count_del * (100/len(column_string))
        #print "delete indel"
        
    
    return perc

def main(folder):
    #read file
    path = '../data/'
    
   # print os.listdir('../data')  
   # for infile in glob.glob( os.path.join(path, 'test/*/*.msl') ):

    for infile in glob.glob( os.path.join(path, folder+'/*/*.msl') ):
        print "**current file is: " + infile + "**"
        size = os.path.getsize(infile)
        if size == 0:
            sys.stderr.write("The file '"+infile+"' is empty, the program will continue with the next non empty file.\n")
            time.sleep(3)
            continue
       
        sequences =  []
        indices = []
        columns = list()
        
        input_handle = open(infile, "rU")
        record_iter = SeqIO.parse(input_handle, 'fasta')
         

        #check if it is not fasta the biopython returns null iterator
        record_iter, walk2 = itertools.tee(record_iter)
        sentinel = object()
        if next(walk2, sentinel) is sentinel:
            sys.stderr.write('The file '+infile+' is not in fasta format.\n The program will continue with the next file.\n')
            continue
        
        
        
        
        
        
        for seq_record in record_iter:  
            sequences.append(seq_record)
            
        len_columns = len(list(sequences[0].seq))
        
        if len_columns<20:
            sys.stderr.write("The file "+infile+" contains very short sequences.\n The program will continue with the next file.\n")
            continue
        
        for col_ind in range(0,len_columns):
            column = list() 
            for record in sequences:
                record_char = list(record.seq)
                column.append(record_char[col_ind])
            columns.append(column)
            #print column
        index = 0;
        for col in columns:
            perc = uniqueCharsPercentage(col)
            #print str(perc) +" unique"
            if perc >50:
                indices.append(index)
            perc = indelsPercentage(col)
            #print str(perc) +" indels"
            if perc >=50:
                indices.append(index)
            index +=1
        print str(indices) +" Indexes with noise"
        
        final_output = []
        newfile = ""
        if  folder != 'test_data':
            newfile = infile.replace(folder,'new')
        else:
                newfile = "../data/test_data/own/new_"+basename(infile)
        output_handle = open(newfile, "w");
        for seqs in sequences:
            chars = list(seqs.seq)
            
            for ind in indices:
                chars[ind]=''
            #print ''.join(chars)
            new_column = ''.join(chars)
            if len(new_column)==0:
                sys.stderr.write("The sequence has no columns left.\n Program will exit.")
                sys.exit(-1)
            record = SeqRecord(Seq(''.join(chars),
                   IUPAC.protein),
                   id=seqs.id, name=seqs.name,
                   description=seqs.description)
            final_output.append(record)
        SeqIO.write(final_output, output_handle, "fasta")
        output_handle.close()
    
    return 0
folder = 'old'
sts = main(folder)
if sts ==0:
    sys.stdout.write("The  program finished the calculations.\n Program exists.\n")
#test = ['S','-','V','V','L','A','K','V','N','K','Q','P','H','H','L','E']
#
#c = uniqueCharsPercentage(test)
#print str(c)+'% :counter'
#c = indelsPercentage(test)
#print str(c)+'% :counter'