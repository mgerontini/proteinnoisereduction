import subprocess
from subprocess import PIPE,STDOUT,Popen
import dendropy
import os
import glob
import sys


def calc(infile, reference_tree):
	string = 'fastprot ' + infile + ' | fnj -O newick'
	print string
	p = subprocess.Popen([string], shell=True, stdout=subprocess.PIPE)

	f2 = open(reference_tree, 'r')
	s = ''
	for line in f2:
		s = s + line

	tree1 = dendropy.Tree.get_from_string(p.communicate()[0], 'newick')
	tree2 = dendropy.Tree.get_from_string(s, 'newick')
	#print tree1
	#print tree2

	val1 = tree1.symmetric_difference(tree2)
	#val1 = tree1.robinson_foulds_distance(tree2)
	f2.close()
	return val1



def main():
	path = '../data/'
	f1 = open('../results/Result2.txt' , 'w')
   	#print os.listdir('../data')  
    	for infile in glob.glob( os.path.join(path, '*/*/*.msl') ):
		string = infile.split('/')
		reference_tree = '../data/old/' + string[3] + '/' + string[3] + '.tree'
		#print reference_tree
		val1 = calc(infile, reference_tree)
		#print val1, val2
		f1.write(string[3] + ' ' + string[4] +' ' + string[2]  + ' ' +  str(val1) + '\n' )
	
	f1.close()

main()



