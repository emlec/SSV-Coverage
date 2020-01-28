#!/usr/bin/env python3
# -*- coding: utf-8 -*

"""
@package    SSV-Coverage
@brief      Contain the template of the empty configuration file for SSV-Coverage
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Mathieu Bolteau - 2018
* <mathieu.bolteau1@gmail.com>
* <mathieu.bolteau@univ-nantes.fr>
* [Github] (https://github.com/mablt)
* [INSERM UMR 1089] (http://umr1089.univ-nantes.fr/)
"""


def write_conf_file(number_bam):
    """
    Create an empty configuration file for the program SSV-Coverage. By default, the file contains one
    bam section. 
    """

    # Initialize variables
    counter = 1
    bam_to_add = str()

    # Write parameters for reference informations
    while counter <= number_bam - 1:
        bam_to_add += """
[Bam{}]

sample_name: 
bam_path:
line_color:
line_width:\n\n""".format(str(counter + 1))
        counter += 1

    # Write the empty configuration file
    with open("Conf.txt", 'w') as f:
        f.write("""

#########################################################################################
#                                                                                       #
#                          SSV-COVERAGE CONFIGURATION FILE                              #
#                                                                                       #
#########################################################################################

#########################################################################################
#                                 INSTRUCTIONS                                          #
#########################################################################################

#   Values can by modify according to the samples, but the file template must remain
#   unchanged, otherwise the program will not be able to load default values.
#   
#   File path should be indicated as absolute path preferably and should not contain
#   blank spaces.  
                                                                                       


#########################################################################################
#                                 PARAMETERS                                            #
#########################################################################################

##################
#   REFERENCE    #
##################

#   These section correspond to the definition of the reference. It is defined by a path
#   to a fasta file that can contain only one sequence. SSV-Coverage program is 
#   responsible for determining the name and the length of the reference. 

[Ref]

#   Path to the reference fasta file that can contain multiple sequences. (String)
#   Fasta file must be in fa.gz format.
ref_path : 


##################
#   BAM FILES    #
##################

#   Each Bam file must be defined by its name (sample_came) and its path (bam_path).
#   For each sample, it is possible to specify the line color and the line width of the curve. 
#   The number of samples is limited to 10.

[Bam1]

#   Name of the sample. (String)
sample_name : 

#   Path to the bam file. (String)
bam_path : 

#   Line color of the sample curve. Color examples : blue, green, red, cyan, magenta, yellow, black, white.
#   If leaved blank 'black' will be used by default. (String)
line_color : 

#   Line width of the sample curve. If leaved blank, '1' will be used by default. (Int)
line_width : 

"""
+ bam_to_add +
"""
##################
#   DEPTH FILES    #
################## 

[Depth]

#   Depth files can be created from BAM files using bedtools genomecov (bedtools) or pysam. 
#   Write either 'bedtools' or 'pysam'. 
#   If leaved blank, bedtools genomecov will be used by default. (String) 
depth_program : 

#####################
#   GRAPH OPTIONS   #
#####################

#   This section contains the options for the graphic design. 
#   If necessary, a zoom can be made. Only one graphic can be created.

[Graph]

#   Title of the graphic. 
#   If leaved blank, 'Coverage along the reference [Name of reference]' will be used by default. (String)
title : 

#   Width of the graphic window. If leaved blank, '6' will be used by default. (String)   
width : 

#   Length of the graphic window. If leaved blank, '8' will be used by default. (String)
length : 

#   Representation of the sequencing coverage after normalization by the total number of sequenced bases.
#   Write 'True' if normalization is required, otherwise 'False'. 
#   If leaved blank, the sequencing coverage will be normalized by default. (Boolean)
depth_normalization: 

#   Scale of the graphic. Scale can be 'log' or 'linear'. 
#   If leaved blank, 'linear' will be used by default. (String)
scale : 

#   x-axis name. If leaved blank, 'Position' will be used by default. (String)
xlabel : 

#   y-axis name. If leaved blank, 'Normalized coverage of sequencing per base' will be used by default. (String)
ylabel :



[Zoom]

#   The zoom concern the x-axis. It is possible to visualize only a region of the DNA sequence.

#   Write either 'True' or 'False' to visualize or not only a region of the DNA sequence.
#   If False, the following options of the section won't be parsed. (Boolean) 
zoom : False

#   Starting position of the region. (Int)
min_position : 

#   Ending position of the region. (Int)
max_position : 

######################
#   OUTPUT OPTIONS   #
######################

#   By default, a svg file containing the graphic is record in the current folder. 
#   It is possible to create a supplementary out file in addition of the svg.

[Output]

#   Name of the output file. If leaved blank, 'out' will be used by default. (String)
output_name : 

[SupplementaryOutput]

#   Convert the svg file to another format (optional). 

#   Write either 'True' or 'False' to create or not a supplementary output file. 
#   If False, the following options of the section won't be parsed. (Boolean)
supplementary_output: 

#   Name of the supplementary output file. If leaved blank 'out' will be used by default. (String)
supplementary_output_prefix :

#   Available format of the supplementary output file: 'png', 'pdf', 'tif', 'ps', 'eps'. (String)
supplementary_output_format : """)
