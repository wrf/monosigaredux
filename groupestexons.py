#!/usr/bin/env python

# v1.0 created 2015-06-16
#
# groupestexons.py use EST information to group exons into plausible transcripts
#
# and using information from:
# https://www.sanger.ac.uk/resources/software/gff/spec.html
# http://www.sequenceontology.org/gff3.shtml

'''
GROUPESTEXONS.PY v1 2015-07-01

  ## GENERAL PRINCIPLES
  assume EST information is always right, or more right than de novo genes
  assume that different ESTs with the same exons are the same transcript
  alternate exon use can be defined by the EST
  ignore antisense transcription
  within a forward reverse pair, any number of exons can be missing

groupestexons.py -e ests.gmap.gff -g gene_models.gff

  specify the delimiter for the ESTs and reads with -n
  default is ., such as for seq123.1 and seq123.2

  ## GETTING GFF RESULTS WITH GMAP

  GMAP should use -f 3 (format 3 gff) as output
  such as:

gmap -d monbr1 Monbr1_cDNA_reads.fasta -f 3 > Monbr1_cdna_reads.gmap.gff

  use multithread option -t for faster results
'''

# TEST PARAMETERS

# TEST CASES
# scaffold_17 167000 XYM20158.x1 forward and reverse overlap, one blast hit as only evidence
# scaffold_17 130200 fgenesh2 170022 fake protein, contains exons going the wrong direction
#    real transcripts are bounded on both ends by ESTs
# scaffold_17 125050 short mapping, appx 20bp to provide false fusion
# scaffold_17 67000 run-on protein containing probably several genes and no RNAseq
# scaffold_17 62000 long EST antisense of several exons, bridging different blast hits
# scaffold_17 43000 two long multiexon genes, one with a single blast hit, other with no evidence

#
import sys
import os
import re
import argparse
import time
from collections import defaultdict,namedtuple
from itertools import groupby

### GFF FORMAT
# GFF3 format is a tabular of:
#
# seqid source type start end score strand phase attributes
#
# in the context of blast results, this is converted as:
# seqid is the subject ID, so sseqid
# source is BLAST
# type could be mRNA, exon, CDS, but in this case HSP as the BLAST result
# start must always be less or equal to end, so if BLAST hits the reverse 
#   complement, these must be switched
# score is the bitscore
# strand is either + or - for forward and reverse hits
# phase is intron phase used only for "CDS" type, otherwise is "."
# attributes is a list of tag=value; pairs usually including query ID

def exon_nums_to_set(exondict):
	for k,v in exondict.iteritems():
		exondict[k] = set(v)

def get_read_bounds(exondict, fnc, ind):
	# returns dictionary of mins for forward reads and maxs for reverse reads
	bounds_d = {}
	for k,v in exondict.iteritems():
		# operates like min( position in tuple for x in interval)
		bounds_d[k] = fnc(x[ind] for x in v)
	return bounds_d

def numbers_to_ranges(numberset):
	# makes groups of value-positions, where values in a row end up as the same group
	for i,l in groupby(enumerate(sorted(numberset)), lambda (x,y): y-x):
		l = list(l)
		# returns generator of first and last positions in the group
		yield l[0][1], l[-1][1]

def make_exon_ranges(exondict):
	exoncounter = 0
	rangeDict = {}
	for k,v in exondict.iteritems():
		exonranges = list(numbers_to_ranges(v))
		rangeDict[k] = exonranges
	# dictionary of lists, where lists contain tuples of exon boundaries
	return rangeDict

def within_exon_range(targetexon, lower, upper):
	return lower <= targetexon <= upper

def get_strand_consensus(strandlist, estscaf):
	numexons = len(strandlist)
	strandconsensus = float(sum(strandlist))/len(strandlist)
	if strandconsensus > 0:
		return "+"
	elif strandconsensus < 0:
		return "-"
	else:
		print >> sys.stderr, "WARNING EST {1} on {2} with {0} exons has exons in both directions".format(numexons, *estscaf)
		return "."

def get_best_part_exon(refgenedict, exon, verbose, exnum, pos, antipos):
	'''for non exact matches, find closest part of an exon'''
	partmatch = namedtuple("partialmatch", "mtype, gexc, extension")
	pm = partmatch(mtype='n', gexc=None, extension=0)
	for gr in refgenedict:
		opm = partmatch(mtype='n', gexc=None, extension=0)
		if exon[pos]==gr[pos]: # pos is either 0 or 1, for start or end of exon or gr
			opm = partmatch(mtype='p', gexc=gr, extension=exon[pos]-gr[antipos])
			if verbose:
				print >> sys.stderr, "Exon {} partial {} at {} {}".format(exnum, exon, gr, refgenedict[gr] )
		elif exon[0]>gr[0] and exon[1]<gr[1]: # for embedded exon
			opm = partmatch(mtype='m', gexc=gr, extension=exon[pos]-gr[antipos])
			if verbose:
				print >> sys.stderr, "Exon {} embedded {} at {} {}".format(exnum, exon, gr, refgenedict[gr] )
		if opm.gexc: # prioritize matched exons, then partial, then embedded
			if pm.gexc:
				if pm.mtype=='m' and opm.mtype=='p': # take partial over embedded
					pm = partmatch(*opm)
				elif pm.mtype=='p' and opm.mtype=='p': # take longer partials
					if pm.extension < opm.extension:
						pm = partmatch(*opm)
			else: # if no pm is found yet, take the first match
				pm = partmatch(*opm)
	# return tuple of best matching exon, or None
	return pm.gexc

def find_matching_exons(exbyreaddict, read, refgenedict, verbose, doreverse):
	'''for each EST exon find any matching gene exons'''
	# each kept exon is in format of (exon to call line, new attribute ID, possible new coords)
	keptexons = []
	lastexon = len(exbyreaddict[read])-1
	# check if any exons are in this direction
	if exbyreaddict[read]:
		readmin = min(e[0] for e in exbyreaddict[read])-100
		readmax = max(e[1] for e in exbyreaddict[read])+100
	else: # otherwise return the empty list
		return keptexons
	# then check if exons are in range
	withinrange = False
	for exon in refgenedict:
		# check if there are exons before the end of the read or after the beginning 
		# or if the read is within an exon
		if readmax >= exon[0] >= readmin or readmax >= exon[1] >= readmin or (readmin >= exon[0] and readmax <= exon[1]):
			withinrange=True
	novelexons = 0
	if verbose:
		print >> sys.stderr, "Checking read {} with exons 0,{}".format(read, lastexon)
	# iterate through exons in sorted order of the direction of the sequence (+/-)
	for i,exon in enumerate(sorted(exbyreaddict[read],key=lambda x: x[0], reverse=doreverse) ):
		if withinrange:
			refex = refgenedict.get(exon, False)
			if refex: # if the exon interval is a perfect match
				if verbose:
					print >> sys.stderr, "Exon {} found {} {}".format(i, exon, refex)
				keepexon = (exon, read, exon)
				keptexons.append(keepexon)
			else:
				if verbose:
					print >> sys.stderr, "No exact match for {}".format(exon)
				#if i==order[0]: # first exon, 5prime
				if i<lastexon:
					if verbose:
						print >> sys.stderr, "Using middle exon {} at {}".format(i, exon)
					keepexon = (exon, read, exon)
					keptexons.append(keepexon)
				# case of only one exon where order is the same
				elif i==lastexon: # last exon, 3prime
					if doreverse:
						bestmatch = get_best_part_exon(refgenedict, exon, verbose, exnum=i, pos=1, antipos=0)
						if bestmatch:
							keepexon = (bestmatch, read, (bestmatch[0], exon[1]) ) # reverse strand
						else:
							if verbose:
								print >> sys.stderr, "Using terminal exon {} at {}".format(i, exon)
							keepexon = (exon, read, exon)
					else:
						bestmatch = get_best_part_exon(refgenedict, exon, verbose, exnum=i, pos=0, antipos=1)
						if bestmatch:
							keepexon = (bestmatch, read, (exon[0], bestmatch[1]) ) # forward strand
						else:
							if verbose:
								print >> sys.stderr, "Using terminal exon {} at {}".format(i, exon)
							keepexon = (exon, read, exon)
					if verbose:
						print >> sys.stderr, keepexon
					keptexons.append(keepexon)
		else:
			# for case where no exons are within the boundaries of the read
			# suggesting a novel transcript
			if verbose:
				print >> sys.stderr, "No exons within range, keeping {}".format(exon)
			keepexon = (exon, read, exon)
			keptexons.append(keepexon)
			novelexons += 1
	if verbose and novelexons==len(exbyreaddict[read]):
		print >> sys.stderr, "NOVEL transcript for {}".format(read)
	return keptexons

def print_updated_gff(matchingexons, ref_lines, scaffold, transstrand, program):
	# generate transcript line
	transstart = str(min(me[2][0] for me in matchingexons) )
	transend = str(max(me[2][1] for me in matchingexons) )
	transattrs = 'gene_id "{0}"; transcript_id "{0}.1";'.format(matchingexons[0][1])
	translist = [scaffold, program, "transcript", transstart, transend, "100", transstrand, ".", transattrs]
	print >> sys.stdout, "\t".join(translist)
	# need dummy line in case 
	dummylist = "\t".join([scaffold, "", "exon", "", "", ".", transstrand, ".", ""])
	exoncount = 0
	# me is a tuple of matching ref exon, read for updating attributes, and possible new positions
	# set is used to remove duplicate exons from forward and reverse
	#for me in sorted(set(matchingexons), key=lambda x: x[2][0]):
	for me in matchingexons: # should be presorted
		exoncount += 1
		outsplits = ref_lines.get(me[0],dummylist).split("\t")
		outsplits[1] = program
		outsplits[3] = str(me[2][0]) # from tuple like (1160, 1498) should be 1160
		outsplits[4] = str(me[2][1]) # and 1498
		outsplits[6] = transstrand
		attributes = '{} exon_number "{}";'.format(transattrs, exoncount)
		outsplits[8] = attributes
		outline = "\t".join(outsplits)
		print >> sys.stdout, outline

def get_fasta_dict(fastafile):
	seqdict = defaultdict(str)
	for line in open(fastafile, 'r'):
		line = line.rstrip()
		if line[0] == ">":
			seqid = line[1:]
		else:
			seqdict[seqid] += line
	return seqdict

def find_polya(seqstring, polyalen):
	seqtail = seqstring[-polyalen:]
	return float(seqtail.count("A"))/polyalen >= 0.8

def get_est_directions(fastafile, polyalen):
	directdict = {}
	seqdict = get_fasta_dict(fastafile)
	for k,v in seqdict.iteritems():
		if find_polya(v, polyalen):
			directdict[k] = True
	return directdict

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-e','--ests', help="EST gmap results file")
	parser.add_argument('-f','--fasta', help="optional fasta file of raw ESTs for polyA detection")
	parser.add_argument('-g','--genes', help="gene models that will be clustered or corrected")
	parser.add_argument('-b','--blast', help="optional tabular blast for direction evidence")
	parser.add_argument('-a','--poly-a', type=int, default=25, help="terminal bases to check for polyA [25]")
	parser.add_argument('-l','--map-length', type=int, default=100, help="minimum mapping length in bases [100]")
	parser.add_argument('-i','--intron-distance', type=int, default=5000, help="max distance for introns [5000]")
	parser.add_argument('-G','--gene-length', type=int, default=50000, help="max length allowed for genes [50000]")
	parser.add_argument('-T','--trans-length', type=int, default=20000, help="max length allowed for transcripts [20000]")
	parser.add_argument('-r','--read-limit', type=int, default=10, help="max allowed reads per EST group [10]")
	parser.add_argument('-n','--number-split', default=".", help="delimited for EST numbers [.]")
	parser.add_argument('-p','--program', help="program for 2nd column in output [EST]", default="EST")
	parser.add_argument('-o','--overlap', help="keep overlapping exons", action="store_false")
	parser.add_argument('-E','--extract', help="filter by a specific EST for debugging")
	parser.add_argument('-S','--select', help="filter by a specific scaffold for debugging")
	parser.add_argument('-v','--verbose', help="verbose output", action="store_true")
	args = parser.parse_args(argv)

	# counter for number of lines, and inferred exons
	linecounter = 0
	exoncounter = 0
	# counter for number of hits that are filtered
	badhits = 0

	starttime = time.time()

	# PART 1
	#
	# make this a dictionary of lists, containing a tuple of exon start and stop positions
	print >> sys.stderr, "Starting exon parsing on %s" % (args.genes), time.asctime()
	gene_ref_exons = defaultdict(dict) # keys are scaffolds, then tuples
	gene_ref_lines = defaultdict(dict) # keys are scaffolds, then lists from lsplits
	for line in open(args.genes,'r'):
		linecounter += 1
		line = line.rstrip()
		if line and not line[0]=="#":
			lsplits = line.split("\t")
			if lsplits[2]=="exon":
				exoncounter += 1
				scaffold = lsplits[0]
				strand = lsplits[6]
				exonbound = ( int(lsplits[3]),int(lsplits[4]) )
				if strand=='+':
					gene_ref_exons[scaffold][exonbound] = 1
				else: # assume '-' strand
					gene_ref_exons[scaffold][exonbound] = -1
				# keep track of all lists for later output
				gene_ref_lines[scaffold][exonbound] = line
			elif lsplits[2]=="CDS":
				pass ### TODO maybe this is needed
	print >> sys.stderr, "Parsed %d lines" % (linecounter), time.asctime()
	print >> sys.stderr, "Counted %d putative exons" % (exoncounter), time.asctime()

	# OPTIONAL
	#
	if args.blast:
		blasthitcount = 0
		blasthitdict = defaultdict(dict)
		for line in open(args.blast, 'r'):
			blasthitcount += 1
			qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.rstrip().split("\t")
			# convert strings of start and end to integers for calculations
			isend = int(send)
			isstart = int(sstart)
			iqend = int(qend)
			iqstart = int(qstart)

			# is checks are passed, accumulate that hit an info into dict
			hittuple = (iqstart, iqend, isstart, isend)
			# must use setdefault since get(qseqid, []) returns None, and cannot append
			blasthitdict[sseqid].setdefault(qseqid, []).append(hittuple)

	# OPTIONAL
	#
	if args.fasta:
		fastaestcount = 0
		print >> sys.stderr, "Searching ESTs for poly-A tails in {}".format(args.fasta), time.asctime()
		directiondict = get_est_directions(args.fasta, args.poly_a)
		print >> sys.stderr, "Found direction for {} ESTs".format(len(directiondict) ), time.asctime()

	# PART II
	#
	# read EST mappings
	plusexons, minusexons = 0,0	
	plusstrand, minusstrand = 0,0
	# counter for potentially overlapping hits or tandem duplicates
	tandemdups = 0
	hitDictCounter = defaultdict(int)

	#print >> sys.stderr, "Reading database sequences from %s" % (args.db), time.asctime()
	#dbdict = SeqIO.to_dict(SeqIO.parse(args.db, 'fasta'))
	#print >> sys.stderr, "Counted %d db sequences" % (len(dbdict)), time.asctime()

	#est_exons = defaultdict(list)
	estexoncount = 0
	readlengths = defaultdict(int) # sum of exons by readname
	estreadsbyscaffold = defaultdict( lambda: defaultdict(list) ) # list of all possible readnames for each estname
	exonsbyscaffold = defaultdict( lambda: defaultdict(list) ) # by scaffold then readname
	estbyread = {}
	forwardreads = defaultdict(list)
	reversereads = defaultdict(list)
	#scafreadcount = defaultdict( lambda: defaultdict(int) )

	readidre = "ID=([\w\d.|:]+);" # must have . | and : for different EST types
	print >> sys.stderr, "Reading EST sequences from %s" % (args.ests), time.asctime()
	for line in open(args.ests, 'r'):
		line = line.rstrip()
		if line and not line[0]=="#":
			lsplits = line.split("\t")
			attrs = lsplits[8]
			# from ID=XYM14183.y1.path1;Name=XYM14183.y1;Target=XYM14183.y1 88 130;Gap=M43
			# or ID=jgi|JGI_XYM10825.rev|.path1;Name=jgi|JGI_XYM10825.rev|;Target=jgi|JGI_XYM10825.rev|
			# or ID=3720288:1.path1;Name=3720288:1;Target=3720288:1 12 684;Gap=M673
			readid = re.search(readidre,attrs).group(1) # should be XYM14183.y1.path1
			readname, readpath = readid.rsplit(".",1) # should be XYM14183.y1 and path1
			estname = readname.split(args.number_split)[0] # should be XYM14183

			if not readpath=="path1": # this takes only the first path
				continue
			if args.extract: # added for debugging specific ESTs
				if not estname==args.extract:
					continue

 			scaffold = lsplits[0]
			strand = lsplits[6]
			estbyread[readname] = estname
			estreadsbyscaffold[scaffold][estname].append(readname)
			exonbound = ( int(lsplits[3]),int(lsplits[4]) )
			exonsbyscaffold[scaffold][readname].append(exonbound)
			#est_exons[scaffold].extend(range(*exonbound) ) # ests are treated strand-non-specific
			readlengths[readname] += exonbound[1]-exonbound[0]+1
			#scafreadcount[scaffold][readname] += 1
			if strand=='+':
				forwardreads[readname].append(exonbound)
			else: # assume '-' strand
				reversereads[readname].append(exonbound)
			estexoncount += 1
	print >> sys.stderr, "Counted %d forward and %d reverse reads" % (len(forwardreads), len(reversereads) ), time.asctime()
	print >> sys.stderr, "Counted %d ESTs with %d bases" % (len(readlengths), sum(readlengths.values()) ), time.asctime()
	print >> sys.stderr, "Counted %d putative exons" % (estexoncount), time.asctime()

	#forwardbounds = get_read_bounds(forwardreads, fnc=min, ind=0)
	#reversebounds = get_read_bounds(reversereads, fnc=max, ind=1)

	matchcount = 0
	transcriptcount = 0
	print >> sys.stderr, "Finding predicted exons that match ESTs", time.asctime()
	for sc,ebr in exonsbyscaffold.iteritems(): # scaffold, dicts of exons by read
		if args.select: # added for debugging of scaffolds
			if not sc==args.select:
				continue

		readgroups = defaultdict(list) # should be list by final transcript boundaries
		readonscaffold = ebr.keys() # take all reads for that scaffold
		estonscaffold = list(set(estbyread[r] for r in readonscaffold)) # take all ESTs for reads on that scaffold
		for est in estonscaffold:
			forwardmatches = []
			reversematches = []
			matchingexons = []
			if args.verbose:
				print >> sys.stderr, "Sorting {} on scaffold {}".format(est, sc)
			if args.extract: # added for debugging specific ESTs
				if not est==args.extract:
					continue
			readsperest = list(set(estreadsbyscaffold[sc][est]))
			if len(readsperest) > args.read_limit:
				print >> sys.stderr, "WARNING counted {} reads from {}, skipping".format(len(readsperest), est)
				continue
			for read in readsperest:
				if args.verbose:
					print >> sys.stderr, "Sorting exons from {}".format(read)
				# matches are in format of (exon, read, new values for position)
				forwardmatches.extend(find_matching_exons(forwardreads, read, gene_ref_exons[sc], args.verbose, doreverse=False) )
				reversematches.extend(find_matching_exons(reversereads, read, gene_ref_exons[sc], args.verbose, doreverse=True) )
			try:
				forward_end = max(me[2][1] for me in forwardmatches)
			except ValueError: # for empty list max fails, so no forward read
				forward_end = 0
			try:
				reverse_end = min(me[2][0] for me in reversematches)
			except ValueError: # for empty list min fails, meaning no reverse read
				reverse_end = 0
			matchingexons.extend(forwardmatches)
			matchingexons.extend(reversematches)

			# add exons in between two ends of reads
			### TODO add optional boundaries for cases with only one read sense
			if forward_end and reverse_end:
				if reverse_end > forward_end: # only if reverse_end is greater than forward_end
					if args.verbose:
						print >> sys.stderr, "Checking internal exons between {} and {}".format(forward_end, reverse_end)
					betweenmatches = []
					for exon in gene_ref_exons[sc]:
						if exon[0]>forward_end and exon[1]<reverse_end:
							### TODO also check against intron distance
							if args.verbose:
								print >> sys.stderr, "De novo internal exon found {}".format(exon)
							midmatch = (exon, est, exon)
							betweenmatches.append(midmatch)
					matchingexons.extend(betweenmatches)

			if matchingexons:
				sortedmatches = sorted(set(matchingexons), key=lambda x: x[2][0])
				# remove exons that span introns of other exons
				if args.overlap:
					poplist = []
					for i,mep in enumerate(zip(sortedmatches, sortedmatches[1:] ) ):
					#	if mep[0][2][1] >= mep[1][2][0]:
						if mep[0][2][1] >= mep[1][2][1]: # if end of one match is after end of another match
							poplist.append(i)
					for pi in reversed(poplist):
						if args.verbose:
							print >> sys.stderr, "Removing overlapping exon {} for {}".format(pi, est)
						sortedmatches.pop(pi)

				# before adding to counts, do sanity check on transcript
				genelength = max(me[2][1] for me in sortedmatches) - min(me[2][0] for me in sortedmatches) + 1
				if genelength > args.gene_length:
					print >> sys.stderr, "WARNING GENE {} on {} is too long: {} with {} exons".format(est, sc, genelength, len(sortedmatches) )
					continue # should skip to next est
				# before adding to counts, do sanity check on transcript
				translength = sum(me[2][1]-me[2][0] for me in sortedmatches)
				if translength > args.trans_length:
					print >> sys.stderr, "WARNING EST {} on {} is too long: {} with {} exons".format(est, sc, translength, len(sortedmatches) )
					continue # should skip to next est

				transcriptcount += 1
				matchcount += len(sortedmatches)

				# determine transcript direction, using poly-A, BLAST or de novo exons
				if args.fasta and directiondict:
					for read in readsperest:
						if directiondict.get(read,False) and forwardreads.get(read,False):
							transstrand = "+"
							if args.verbose:
								print >> sys.stderr, "Found poly-A for {} strand on {}".format(transstrand, est)
							break
						elif directiondict.get(read,False) and reversereads.get(read,False):
							transstrand = "-"
							if args.verbose:
								print >> sys.stderr, "Found poly-A for {} strand on {}".format(transstrand, est)
							break
					else:
						transstrand = get_strand_consensus([gene_ref_exons[sc].get(x[0],1) for x in sortedmatches ], (est,sc) )
				else:
					# forces '+' if there is no gene model
					### TODO incorporate evidence from blast here
					transstrand = get_strand_consensus([gene_ref_exons[sc].get(x[0],1) for x in sortedmatches ], (est,sc) )
				print_updated_gff(sortedmatches, gene_ref_lines[sc], sc, transstrand, args.program)
	print >> sys.stderr, "Counted %d matching exons or partial exons" % (matchcount), time.asctime()
	print >> sys.stderr, "Printed %d transcripts" % (transcriptcount), time.asctime()
	print >> sys.stderr, "# Processed completed in %.1f minutes" % ((time.time()-starttime)/60)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
