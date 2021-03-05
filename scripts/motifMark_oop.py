#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------------
# Motif Marking tool
# Author: Cameron Watson
# Last Updated: 9 Feb 2021
#--------------------------------------------------------------------------------------------------------------

import argparse
import re
import cairo
import numpy as np

#--------------------------------------------------------------------------------------------------------------
# USER INPUT
#--------------------------------------------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(description = "A tool for marking specified motifs in FASTA files. This program returns \
        an image (SVG) of genes with exons and introns denoted and specified motifs marked.",
    add_help = True)

    parser.add_argument("-f", "--fasta", 
    help = "A FASTA file containing sequences with exons capitalized and all other nucleotides in lower-case", required = True)

    parser.add_argument("-m", "--motif",
    help = "A plain-text file with one motif (<=10 base nucleotide sequence) per line", required = True)

    return parser.parse_args()

args = get_args()

#--------------------------------------------------------------------------------------------------------------
# CLASSES
#--------------------------------------------------------------------------------------------------------------

class Gene:

    '''
    Gene object contains length of the full gene, drawGene draws a thin line based on the gene length
    '''

    __slots__ = ['start', 'length', 'width']

    def __init__(self, record):
        self.start = 0.01 # margin
        self.length = len(record[1])
        self.width = 0.008 # intron width

    def drawGene(self, context, y_coord, surface_width):

        '''
        Draws gene length on surface context, located by y_coord and scaled to the surface width
        '''

        context.set_source_rgb(0,0,0)
        context.set_line_width(self.width)
        context.move_to(self.start, y_coord)
        context.line_to((self.length / surface_width) + 0.01, y_coord)
        context.stroke()

class Exon:

    '''
    Exon object contains all exons for a gene, drawExons draws a thick line based on coords of each exon
    Exon class finds all exons because regex search tool makes it easier than making an object for each
    individual exon
    '''

    __slots__ = ['coordinates', 'width']

    def __init__(self, record):
        self.coordinates = []
        self.findExons(record)
        self.width = 0.05

    def findExons(self, record):

        '''
        Finds all capitalized segments in sequence 
        '''

        exon_matches = re.finditer("[A-Z]+", record[1])
        # extracting coordinates (first char index: last char index) from re.iter objects
        for i in exon_matches:
            self.coordinates.append(i.span())
        return self

    def drawExons(self, context, y_coord, surface_width):

        '''
        Draws all exons on one gene, scaled to width of surface
        '''
        context.set_line_width(self.width)
        for ex in self.coordinates:
            start, stop = ex[0], ex[1] # exon coordinates, to be scaled by image width
            context.move_to((start/surface_width + 0.01), y_coord)
            context.line_to((stop/surface_width + 0.01), y_coord)
            context.stroke()

class Motif:

    __slots__ = ['motif', 'reMotif', 'coords', 'length', 'width', 'color']

    def __init__(self, motif, counter):

        # best approximation of the Darjeeling limited color pallete
        colors = [(.008,.47,.69),
                (.68,.22,0.09),
                (.117,.75,.803),
                (0.82,.61,.18),
                (.29,.6,0.49),
                (.35,.7,.9),
                (0,.6,.5)]

        self.motif = motif
        self.reMotif = ""
        self.expandMotif()
        self.length = len(self.motif)
        self.coords = []
        self.color = colors[counter]

    def expandMotif(self):

        '''
        Generates a Regex motif based on original motif and IUPAC dictionary
        '''

            # RNA/DNA indifferent dictionary for IUPAC degenerate bases
        iupac_degens = { 
            "w" : "[atuATU]", "b" : "[cgtuCGTU]", "r" : "[agAG]", "t" : "[tuTU]",
            "s" : "[cgCG]", "d" : "[agtuAGTU]", "n" : "[acgtuACGTU]", "u" : "[tuTU]",
            "m" : "[acAC]", "h" : "[actuACTU]", "y" : "[ctuCTU]", "a" : "[aA]",
            "k" : "[gtuGTU]", "v" : "[acgACG]", "z" : "-", "c" : "[cC]", "g" : "[gG]"
        }

        temp_motif = self.motif.lower()
        for i in temp_motif:
            self.reMotif += iupac_degens[i] # rebuild input motif with IUPAC regex from dict
        return self

    def locateMotif(self, record):
        '''
        Finds all instances of one motif object in a sequence
        '''
        seq = record[1]
        counter = 0
        for nuc in range(len(seq)): # iterates through each base in seq
            # prevent indexing past length of seq
            if (len(seq) - counter) > (self.length -1): 
                if re.fullmatch(self.reMotif, seq[nuc:(nuc+self.length)]):
                    # (start pos, stop pos)
                    self.coords.append((nuc,(nuc+self.length)))
            counter += 1
        return self

    def drawMotif(self, context, y_coord, surface_width):

        '''
        Draws all instances of motif on gene
        '''
        context.set_source_rgba(self.color[0], self.color[1], self.color[2], 0.7)
        # places motifs on gene same way as exons
        for coord in self.coords:
            start, stop = coord[0], coord[1]
            context.move_to((start/surface_width + 0.01), y_coord)
            context.line_to((stop/surface_width + 0.01), y_coord)
            context.stroke()

class FastaHeader:

    __slots__ = ['start', 'text']

    def __init__(self, record):
        self.start = 0.01
        self.text = record[0]
        self.parseHeader() # calls function to update self.text

    def parseHeader(self):

        '''
        Gets rid of fasta carrot and unnecessary header info
        '''
        longName = re.split(">", self.text)[1]
        self.text = re.split(" ", longName)[0]
        return self

    def drawHeader(self, context, space):

        '''
        Writes gene label slightly offset up from gene drawing
        '''
        context.set_source_rgb(0,0,0)
        context.move_to(self.start, space - 0.1)
        context.select_font_face('Calibri')
        context.set_font_size(.02)
        context.show_text(self.text)
        return None


class GeneGroup:

    __slots__ = ['rank', 'record', 'gene', 'exon', 
                'motifs', 'header', 'motifList']

    def __init__(self, record, record_count, motif_list):

        '''
        defines instance of class GeneGroup, some attributes dependent on other classes
        Dependencies: FastaHeader, Gene, Exon, Motif
        '''
        self.record = record
        self.rank = record_count
        self.header = FastaHeader(self.record)
        self.gene = Gene(self.record)
        self.exon = Exon(self.record)
        self.motifList = motif_list
        self.motifs = []
        self.callMotifs()

    def callMotifs(self):

        '''
        Calling Motif class locateMotif method for each motif in motif list
        '''
        num = 0
        for i in self.motifList:
            temp_motif = Motif(i, num)
            temp_motif = temp_motif.locateMotif(self.record)
            self.motifs.append(temp_motif)
            num += 1
        return self

    def draw(self, context, spacer, width):

        '''
        Draw everything
        '''
        surface_space = self.rank * spacer
        self.header.drawHeader(context, surface_space)
        self.gene.drawGene(context, surface_space, width)
        self.exon.drawExons(context, surface_space, width)
        jitter = 0 # offset for legend lines
        for motif in self.motifs:
            motif.drawMotif(context, surface_space, width)
            # putting legend generator here for now, will functionalize later
            color = motif.color
            context.set_source_rgb(color[0],color[1],color[2])
            context.set_font_size(.02)
            context.move_to(0.8, 0.1+jitter)
            context.show_text(motif.motif)
            jitter += 0.03



#--------------------------------------------------------------------------------------------------------------
# MAIN
#--------------------------------------------------------------------------------------------------------------

# open input motif file, store motifs in list

input_motifs = []
with open(args.motif, "r") as fh:
    for line in fh:
        input_motifs.append(line.strip("\n"))

print("parsed input motifs")

# iterate through input fasta file

fasta_fh = open(args.fasta, "r")

print("fasta opened")

ln = 0
all_geneGroups = []

for line in fasta_fh:
    if ln == 0: # captures first line
        current_record = [line.strip("\n"),""]
        rec_ctr = 1
    elif re.match("^>", line): # captures all subsequent header lines
        current_geneGroup = GeneGroup(current_record, rec_ctr, input_motifs)
        all_geneGroups.append(current_geneGroup)
        current_record = [line.strip("\n"),""] # reset and restart for next record
        rec_ctr += 1
    else: # sequence lines
        current_record[1] += line.strip("\n")
    ln += 1

# handle EOF, close input fasta

current_geneGroup = GeneGroup(current_record, rec_ctr, input_motifs)

all_geneGroups.append(current_geneGroup)

fasta_fh.close()

print("fasta closed")

# output naming

basename = re.split("/", args.fasta)[-1]

prefix = re.split("\.", basename)[0]

outName = "./" + prefix + ".svg"

# initializing drawing surface, dimensions based on longest gene

print("initializing drawing surface")

longest_gene = 0

for ggroup in all_geneGroups: # find longest gene
    if ggroup.gene.length > longest_gene:
        longest_gene = ggroup.gene.length

if longest_gene == 0:
    raise ValueError("Check sequence lengths")

width = longest_gene * 1.1
height = width * 0.75

lnSpacing = 1 / (rec_ctr + 1)

surface = cairo.SVGSurface(outName, width, height)

context = cairo.Context(surface)

context.scale(width, height) # scale so everything is 0-1

# drawing all gene groups on surface

print("drawing all groups")

for ggroup in all_geneGroups:
    ggroup.draw(context, lnSpacing, width)

surface.finish()