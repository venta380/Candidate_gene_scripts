"""
###python scrypt to obtain concensus sequence from readcount.
###this script will provide concensus based on the ocuuring nucleotides in a specific position and convert them to IUPAC CODES.
##dependencies: python2
##usage: python2 getfasta2.py [readcount.txt] > output.fasta
##example: python getfasta2.py test.readcount > test.fasta
##note: i being a basic script can't deel with indels
#####################################
#####        Venkat Talla       #####
#####   Venkat.Talla@ebc.uu.se  #####
#####        Backstrom Lab      #####
#####################################
"""


import sys
import string


infile = sys.argv[1]
limit = int(sys.argv[2])
table = [line.strip().split() for line in open(infile)]
seq=""
i=int(table[0][1])
if i != 1:
	sequence=str("N"*int(i-2))
else:
	sequence=""

with open(infile) as f:
	for line in f:
		bases=[]
		pos = int(line.split()[1])
		if pos != i:
			info=line.split()[5:]
			indel=line.split()[9:]
			for item in indel:
				base6=[]
				max_base_indel={}
				if item.split(':')[0][0] != "N":
					if item.split(':')[0][0] == '+':
						base6.append(int(item.split(':')[1]))
						base_indel=max(base6) 
						max_base_indel[int(item.split(':')[1])]=item.split(':')[0][1:]
						base_indel_to_append=max_base_indel[base_indel]
					if item.split(':')[0][0] == '-':
						base6.append(int(item.split(':')[1]))
						base_indel=max(base6)
						max_base_indel[int(item.split(':')[1])]=item.split(':')[0][1:]
						base_indel_to_append=""
				else:
					base_indel=0
					bases=[]
 			base1=int(info[0].split(':')[1])
			base2=int(info[1].split(':')[1])
			base3=int(info[2].split(':')[1])
			base4=int(info[3].split(':')[1])
			base5=int(info[4].split(':')[1])
			if base1 != 0:
				bases.append(info[0].split(':')[0])
			if base2 != 0:
				bases.append(info[1].split(':')[0])
			if base3 != 0:
				bases.append(info[2].split(':')[0])
			if base4 != 0:
				bases.append(info[3].split(':')[0])
			if base5 != 0:
				bases.append(info[4].split(':')[0])
			if "A" in bases: base_in_the_pos = "A"
			if "C" in bases: base_in_the_pos = "C"
			if "G" in bases: base_in_the_pos = "G"
			if "T" in bases: base_in_the_pos = "T"
			if "A" in bases and "G" in bases: base_in_the_pos = "R"
			if "C" in bases and "T" in bases: base_in_the_pos = "Y"
			if "G" in bases and "C" in bases: base_in_the_pos = "S"
			if "A" in bases and "T" in bases: base_in_the_pos = "W"
			if "G" in bases and "T" in bases: base_in_the_pos = "K"
			if "A" in bases and "C" in bases: base_in_the_pos = "M"
			if "C" in bases and "G" in bases and "T" in bases: base_in_the_pos = "B"
			if "A" in bases and "G" in bases and "T" in bases: base_in_the_pos = "D"
			if "A" in bases and "C" in bases and "T" in bases: base_in_the_pos = "H"
			if "A" in bases and "G" in bases and "C" in bases: base_in_the_pos = "V"
			if "A" in bases and "G" in bases and "C" in bases and "T" in bases: base_in_the_pos = "N"
			sequence += str("N"*((pos-prev_positon)-1))+base_in_the_pos
			i = pos+1
		else:
			info=line.split()[5:]
			indel=line.split()[9:]
			for item in indel:
				base6=[]
				max_base_indel={}
				if item.split(':')[0][0] != "N":
					if item.split(':')[0][0] == '+':
						base6.append(int(item.split(':')[1]))
						base_indel=max(base6) 
						max_base_indel[int(item.split(':')[1])]=item.split(':')[0][1:]
						base_indel_to_append=max_base_indel[base_indel]
					if item.split(':')[0][0] == '-':
						base6.append(int(item.split(':')[1]))
						base_indel=max(base6)
						max_base_indel[int(item.split(':')[1])]=item.split(':')[0][1:]
						base_indel_to_append=""
				else:
					base_indel=0
 			base1=int(info[0].split(':')[1])
			base2=int(info[1].split(':')[1])
			base3=int(info[2].split(':')[1])
			base4=int(info[3].split(':')[1])
			base5=int(info[4].split(':')[1])
			if base1 != 0:
				bases.append(info[0].split(':')[0])
			if base2 != 0:
				bases.append(info[1].split(':')[0])
			if base3 != 0:
				bases.append(info[2].split(':')[0])
			if base4 != 0:
				bases.append(info[3].split(':')[0])
			if base5 != 0:
				bases.append(info[4].split(':')[0])
			if "A" in bases: base_in_the_pos = "A"
			if "C" in bases: base_in_the_pos = "C"
			if "G" in bases: base_in_the_pos = "G"
			if "T" in bases: base_in_the_pos = "T"
			if "A" in bases and "G" in bases: base_in_the_pos = "R"
			if "C" in bases and "T" in bases: base_in_the_pos = "Y"
			if "G" in bases and "C" in bases: base_in_the_pos = "S"
			if "A" in bases and "T" in bases: base_in_the_pos = "W"
			if "G" in bases and "T" in bases: base_in_the_pos = "K"
			if "A" in bases and "C" in bases: base_in_the_pos = "M"
			if "C" in bases and "G" in bases and "T" in bases: base_in_the_pos = "B"
			if "A" in bases and "G" in bases and "T" in bases: base_in_the_pos = "D"
			if "A" in bases and "C" in bases and "T" in bases: base_in_the_pos = "H"
			if "A" in bases and "G" in bases and "C" in bases: base_in_the_pos = "V"
			if "A" in bases and "G" in bases and "C" in bases and "T" in bases: base_in_the_pos = "N"
			sequence += base_in_the_pos
			prev_positon = pos
			i+=1
			
			


if len(sequence)<limit:
        end=limit-len(sequence)
        #print ">" + line.split()[0]
        print sequence+(end*"N")
else:
		#print ">" + line.split()[0]
		print sequence

