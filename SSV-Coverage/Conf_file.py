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


def write_conf_file (number_bam):
    """
    Create an empty configuration file for SSV-Coverage. By default, the file contains one 
    bam section. 
    """

    # Initialize variables
    counter = 1
    bam_to_add = str()

    # Write parameters for reference informations
    while counter<=number_bam-1:
        bam_to_add += """
[Bam{}]

sample_name : 
bam_path :
line_color :
line_width :\n\n""".format(str(counter+1))
        counter+=1

    # Write the empty configuration file
    with open ("Conf.txt", 'w') as f:
        f.write ("""

#########################################################################################
#                                                                                       #
#                          SSV-COVERAGE CONFIGURATION FILE                              #
#                                                                                       #
#########################################################################################

#########################################################################################
#                                 INSTRUCTIONS                                          #
#########################################################################################

#   Values can by customized with users values, but the file template must remain
#   unchanged, otherwise the program will not be able to load default values.
#   Please follow the following recommendations :
#   
#   - File path should be indicated as absolute path preferably and should not contain
#     blank spaces  
                                                                                       


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
ref_path : 


##################
#   BAM FILES    #
##################

#   These sections correspond to the definitions bam files. They are defined by a name, 
#   a path of the bam files. For each sample, it is possible to customize the curve 
#   which represent this sample. For that, the curve color and the curve width can be 
#   personalize. The number of sample is limited at 10 maximum.

[Bam1]

#   Name of the sample. (String)
sample_name : 

#   Path to the bam file. (String)
bam_path : 

#   Curve color of the sample for the graph. The possible colors are : blue, green, 
#   red, cyan, magenta, yellow, black, white. If leaved blank 'black' will be used as a 
#   scale.(String)
line_color : 

#   Curve width of the sample for the graph. If leaved blank '1' will be used as a scale.
#   (Int)
line_width : 

"""
+bam_to_add+
"""
#####################
#   GRAPH OPTIONS   #
#####################

#   This section groupe options to parameters the graph construction. If necessary a 
#   zoom can be create. When it is flag, only a graph with his parameters is 
#   generate. 

[Graph]

#   Global title of the graph. If leaved blank, 'Coverage along the reference [Name 
#   of reference]' will be used as a scale. (String)
title : 

#   Width of the graph window. If leaved blank, '6' will be used as a scale. (String)   
width : 

#   Length of the graph window. If leaved blank, '8' will be used as a scale. (String)
length : 

#   This value correspond to the scale of the graph. It will be 'log' or 'linear' 
#   only. If leaved blank 'linear' will be used as a scale. (String)
scale : 

#   Name of the x axe title. If leaved blank 'Position' will be used as a scale. 
#   (String)
xlabel : 

#   Name of the y axe title. If leaved blank 'Normalized coverage of sequencing per 
#   base' will be used as a scale. (String)
ylabel :


[Zoom]

#   The zoom concern the x axe. It is possible to select a part of the positions.

#   Flag to activate or deactivate the zoom. If False the following options
#   of the section won't be parsed (Boolean) 
zoom : 

#   Number of the start position for the zoom. (Int)
min_position : 

#   Number of the stop position for the zoom. (Int)
max_position : 

######################
#   OUTPUT OPTIONS   #
######################

#   By default, a svg file containing the graph is record in thecurrent folder. It is 
#   possible to create a supplementary out file in addition of the svg.

[Output]

#   Name of the out file. If leaved blank 'out' will be used as a scale. (String)
output_name : 

[SupplementaryOutput]
#   Create an other out file with the format if it is necessaray.

#   Flag to activate or desactivate the creation of supplementary out file. If False,
#   the following options of the section won't be parsed (Boolean)
supplementary_output: 

#   Name of the supplementary out file. If leaved blank 'out' will be used as a scale. 
#   (String)
supplementary_output_prefix :

#   The possible format of the supplementary out file are : 'png', 'pdf', 'tif',
#   'ps', 'eps'. (String)
supplementary_output_format : """)

