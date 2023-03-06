#!/usr/bin/env python

#Import methods
import argparse
import cairo
import re

## ARGPARSE ###################################################################################
#create input arguments for input samfile, input umi file, output file and he for help, using "h" or "help" causes error
def get_args():
    parser = argparse.ArgumentParser(description = "client input program")
    parser.add_argument("-f", "--file", help = "input fasta filename", type=str, required = True)
    parser.add_argument("-m", "--motifs", help = "input motifs filename", type=str, required = True)
    return parser.parse_args()

#make arg input name easier to type
clinput = get_args()
f = clinput.file
m = clinput.motifs


##### OBJECT ORIENTED ################################################################################
class GeneClass: # create gene class
    def __init__(self, sequence, header): #define class attributes
        self.sequence = sequence
        self.header = header
        self.motif_list = []
        

    def find_exon(self): #create method to find exon
        exon = re.search("[A-Z]+", self.sequence)
        self.exon_location = exon.span()


class MotifClass:
    def __init__(self, motif):
        self.motif = motif.upper()
        self.corr_seq = self.motif
        self.DegenCorrector()

    def FindMotif(self, sequenceObject): #method to search sequence and find motifs
            motifLocation = re.finditer(self.regexMotif, sequenceObject.sequence, re.IGNORECASE)
            
            matched_motifs = []
            for match in motifLocation: #loop through iterable object to store the start location each time
                matched_motifs.append((match.start(), match.start()+len(self.motif)))
            
            return(matched_motifs)
    
    def DegenCorrector(self):
        degen_dict = {"U":"[TU]", "W":"[AT]", "S":"[CG]", "M":"[AC]", "K":"GT", "R":"[AG]", "Y":"[CT]", "B":"[CGT]", "D":"[AGT]", "H":"ACG", "V":"ACG", "N":"ACGT"}
        for key, value in degen_dict.items():
            self.corr_seq = self.corr_seq.upper().replace(key, value)
         
        self.regexMotif = "(?=("+self.corr_seq+"))" #regex to lookahead of our self.corr_seq when it is called


        


###### READ IN FASTA #################################################################################################
#create empty dictionary to store header & sequence
genes = {}

#open fasta file and loop through it stripping the new line so we have a 1 line fasta sequence
with open(f, "rt") as fh:
    for line in fh:
        line = line.strip('\n')
        if line.startswith('>'):
            headline = line.strip('>')
            genes[headline] = ''
        else:
            genes[headline] += line

#Get longest sequence length so we have a max length for pycairo plot
# loop through gene dictionary and find longest sequence

longest_seq = ""

for key, value in genes.items():
    if len(genes[key]) > len(longest_seq):
        longest_seq = genes[key]
    


# iterate through genes dictionary and store in GeneClass
for key, value in genes.items():
    genes[key] = GeneClass(value, key)

    genes[key].find_exon() #at key in dictionary run find_exon method

for key in genes:
    sequenceObject = genes[key]

###### READ IN MOTIFS ######################################################################################################
motif_dict = {}
motif_names_list = []

#loop through motifs and create a dictionary with the motif as the key and a motif object as the value
with open(m, "rt") as fh:
    for line in fh:
        line = line.strip('\n')
        motif_dict[line] = MotifClass(line)
        motif_names_list.append(line)

#create a surface height variable based off of how many sequences we have in our dataset
surf_height = 200
for key, value in genes.items():
    surf_height += 100

#create surface layer
surface = cairo.SVGSurface("MotfiMark.svg", len(longest_seq)+100, surf_height)
#create context layer
context = cairo.Context(surface)

#set contect background to white
context.set_source_rgb(1,1,1)
context.paint()

context.set_line_width(4)

# draw exon legend
context.set_source_rgb(.2,.2,.2)
context.rectangle(50, 120, 50, 10)
context.fill()

# write exon label
context.move_to(50, 115)
context.show_text("Exon")


# create variables that we are going to use to modularize our graph based off in put motif and fasta files
ystart = 200
xstart = 50
recx = 50
context.set_source_rgb(0,0,0)
label_x = 50
seq_label_y = 180


#create list of tuples for different colors
#DRE, for some reason it skips a value??
colors = [(1,0,0), (0,0.5,0), (0,0,1), (1,0,1), (1,1,0)]

#create list of different x-axis values for rectangle
#DRE not using this right now
recstart = 50
counter = 0
name_counter = 0

#loop through genes dictionary and output a line at each one

for key, seq in genes.items():
    #draw our full sequences
    context.set_line_width(4)
    context.set_source_rgb(.5,.5,.5)
    context.move_to(xstart, ystart)
    context.line_to(xstart + len(genes[key].sequence), ystart)
    context.stroke()

    #draw our exons on that sequence 
    context.set_line_width(10)
    context.set_source_rgb(.2,.2,.2)
    context.move_to(xstart + genes[key].exon_location[0], ystart)
    context.line_to(xstart + genes[key].exon_location[1], ystart)
    context.stroke()

    #label our sequences
    context.move_to(xstart, seq_label_y)
    context.show_text(genes[key].header.split(' ')[0])
    seq_label_y += 100

    #DRE, ask for help with lookahead to get overlapping motifs
    #loop through our motifs to plot our motifs, then loops to plot legend and labels
    for i, (motifkey , motifval) in enumerate(motif_dict.items()):
        plot_list = motifval.FindMotif(seq)
        


        for line in plot_list:
            context.set_line_width(20)
            context.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
            context.move_to(xstart + line[0], ystart)
            context.line_to(xstart + line[1], ystart)
            context.stroke()

        #now draw our rectangles for legend, use a counter to stop writing legend squares once we've worked through every motif
        if counter < len(motif_names_list):
            context.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
            context.rectangle(recstart, 60, 25, 25)
            context.fill()

            #loop through our list of input motifs and write the motif name as a label, increase a counter until it is the length of the file
            for motif in motif_names_list:
                if name_counter < len(motif_names_list):
                    context.set_source_rgb(0, 0, 0)
                    context.move_to(label_x, 55)
                    context.show_text(motif)
                    label_x += 150
                    name_counter += 1
            counter += 1
        
        recstart += 150

    ystart += 100
 

#output to png with same title as input fasta file but changing suffix to png
surface.write_to_png(f[:-5] + "png")

#write out drawing
surface.finish()