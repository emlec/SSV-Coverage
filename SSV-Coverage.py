#!/usr/bin/env python3
# -*- coding: utf-8 -*

"""
@package    SSV-Coverage
@brief      Main file of the program
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Mathieu Bolteau - 2018
* <mathieu.bolteau1@gmail.com>
* <mathieu.bolteau@univ-nantes.fr>
* [Github] (https://github.com/mablt)
* [INSERM UMR 1089] (http://umr1089.univ-nantes.fr/)
"""

try:
    # Standard library imports
    import sys
    import optparse  # Mandatory package
    import matplotlib  # Mandatory package

    matplotlib.use('svg')  # Mandatory for the hardcopy backend
    import matplotlib.pyplot as plt  # Mandatory package
    from os import path, listdir  # Mandatory package
    from time import time
    from datetime import datetime
    import configparser  # Mandatory package
    import pysam  # Mandatory package
    from Bio import SeqIO  # Mandatory package

    # Local imports
    from Conf_file import write_conf_file  # if not imported = no creation of configuration file

except ImportError as E:
    print(E)
    print("Please verify your dependencies. See Readme for more informations\n")
    exit()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Main(object):
    """
    Create a graphic representing the coverage of sequencing from bam files.
    """
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    # ~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "SSV-Coverage.py 0.1.0"
    USAGE = "SSV-Coverage.py -c Conf.txt [-i -h]"

    # ~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init(self):
        """
        Init class method for instantiation from command line. Parse arguments.
        """

        # Define parser usage, options
        parser = optparse.OptionParser(usage=self.USAGE, version=self.VERSION)
        parser.add_option('-i', dest="bam_number", type="int",
                          help="Generate an empty configuration file adapted to the number of bam files to create coverage graph and exit. [Mandatory]")
        parser.add_option('-c', dest="conf_file",
                          help="Path to the configuration file [Mandatory]")

        # Parse arguments
        options, args = parser.parse_args()

        # If the program is run without option, show the help and exit
        if len(sys.argv[1:]) == 0:
            parser.print_help()
            exit()

        return Main(options.bam_number, options.conf_file)

    # ~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, bam_number=None, conf_file=None):

        self.start_time = time()

        # Create a empty configuration file if needed with the number of bam files necessary
        if bam_number:
            if bam_number > 1:
                print("Create an empty configuration file, appropriate for {} bam files, in the current folder.".format(
                    bam_number))
            else:
                print("Create an empty configuration file, appropriate for {} bam file, in the current folder.".format(
                    bam_number))
            write_conf_file(bam_number)
            exit()

        # Define a configuration file parser object and load the configuration file
        self.conf_file = conf_file
        self.config = configparser.RawConfigParser(allow_no_value=True)
        self.config.read(conf_file)
        if not self.config.sections():
            print("Empty or invalid configuration file. See Readme for more informations\n")
            exit()

        # Declaration of global values
        self.dict = {
            "sample_name": list(),
            "bam_path": list(),
            "line_color": list(),
            "line_width": list(),
            "x_list": list(),
            "y_list": list(),
            "x_max_list": list(),
            "y_max_list": list()}
        self.title_write = True
        self.max_position_write = True
        self.reference_name = ""
        self.supplementary_output = False
        self.created_files = list()

        try:

            self.reference_path = self.config.get("Ref", "ref_path")
            if self.reference_path:
                self.reference()

            for i in range(1, 10):
                bam_num = "Bam" + str(i)
                # If there are no more bam files, stop the data recovery
                if not self.config.has_section(bam_num):
                    break

                self.dict["sample_name"].append(self.config.get(bam_num, "sample_name"))

                self.dict["bam_path"].append(self.config.get(bam_num, "bam_path"))

                # If line color is not indicate, black is used
                if not self.config.get(bam_num, "line_color"):
                    self.dict["line_color"].append("black")
                else:
                    self.dict["line_color"].append(self.config.get(bam_num, "line_color"))

                # If line width is not indicate, 1 is used
                if not self.config.get(bam_num, "line_width"):
                    self.dict["line_width"].append(1)
                else:
                    self.dict["line_width"].append(self.config.get(bam_num, "line_width"))

            self.scale = self.config.get("Graph", "scale")
            if not self.config.get("Graph", "scale"):
                self.scale = "linear"

            # if (self.scale != "linear" or self.scale != "log"):
            #    print ("The 'scale' option of the configuration file is not in the correct format")
            #    print ("Please report to the descriptions in the configuration file\n")
            #    exit()

            self.xlabel = self.config.get("Graph", "xlabel")
            if not self.xlabel:
                self.xlabel = "Position"

            if not self.config.get("Graph", "ylabel"):
                self.ylabel = "Normalized coverage of sequencing per base"
            else:
                self.ylabel = self.config.get("Graph", "ylabel")

            if not self.config.get("Graph", "title"):
                self.title_write = False
            else:
                self.title = self.config.get("Graph", "title")

            if not self.config.get("Graph", "length"):
                self.length = 8
            else:
                self.length = int(self.config.get("Graph", "length"))

            if not self.config.get("Graph", "width"):
                self.width = 6
            else:
                self.width = int(self.config.get("Graph", "width"))

            if not self.config.get("Output", "output_name"):
                self.output_name = "out"
            else:
                self.output_name = self.config.get("Output", "output_name")

            #   If zoom is flag
            if self.config.getboolean("Zoom", "zoom"):

                # Test if the min position if superior to max position
                if int(self.config.get("Zoom", "min_position")) > int(self.config.get("Zoom", "max_position")):
                    print("One of the value in the 'Zoom' section of the configuration file is not correct")
                    print("Please report to the descriptions in the configuration file\n")
                    exit()

                # If the min position is inferior to 1, means is inferior to the first position (1)
                if int(self.config.get("Zoom", "min_position")) < 1:
                    print("One of the value in the 'Zoom' section of the configuration file is not correct")
                    print("Please report to the descriptions in the configuration file\n")
                    exit()
                else:
                    self.min_position = int(self.config.get("Zoom", "min_position"))

                # If the max position is superior to the length reference (maximum of possible position)
                if int(self.config.get("Zoom", "max_position")) > self.reference_length:
                    print("One of the value in the 'Zoom' section of the configuration file is not correct")
                    print("Please report to the descriptions in the configuration file\n")
                    exit()
                else:
                    self.max_position = int(self.config.get("Zoom", "max_position"))

            else:
                self.min_position = 0
                # Indicate that the zoom is not given
                self.max_position_write = False

            # If supplementary output is flag
            if self.config.getboolean("SupplementaryOutput", "supplementary_output"):
                self.supplementary_output = True

                if not self.config.get("SupplementaryOutput", "supplementary_output_prefix"):
                    self.supplementary_output_prefix = "out"
                else:
                    self.supplementary_output_prefix = self.config.get("SupplementaryOutput",
                                                                       "supplementary_output_prefix")

                self.supplementary_output_format = self.config.get("SupplementaryOutput", "supplementary_output_format")
                if self.supplementary_output_format not in ('png', 'pdf', 'tif', 'ps', 'eps'):
                    print(
                        "The 'supplementary_output_format' option of the configuration file is not in the correct format")
                    print("Please report to the descriptions in the configuration file\n")
                    exit()



        except configparser.NoOptionError as E:
            print(E)
            print("An option is missing in the configuration file")
            print("Please report to the descriptions in the configuration file\n")
            exit()
        except configparser.NoSectionError as E:
            print(E)
            print("An section is missing in the configuration file")
            print("Please report to the descriptions in the configuration file\n")
            exit()
        except ValueError as E:
            print(E)
            print("One of the value in the configuration file is not in the correct format")
            print("Please report to the descriptions in the configuration file\n")
            exit()

        for i in range(0, 10):
            # Try if the bam_path exist in the configuration file
            try:
                self.dict["bam_path"][i]
            # If not exist, stop the execution on the bam_path
            except IndexError:
                break
            # If exist, execute
            else:

                print("\n\n====== SAMPLE {} ======\n".format(self.dict["sample_name"][i]))
                # Recovery of the current directory
                directory = path.dirname(self.dict["bam_path"][i])
                if directory == "":
                    directory = "./"
                list_dir = list()
                # Recovery of the list of files in the current directory
                list_dir = listdir(directory)

                # Recovery of the base name of the bam path and substitute the extension by '.depth'
                bed_name = path.basename(path.splitext(self.dict["bam_path"][i])[0] + ".bed")

                # If depth file does not exist, it will be create
                if not bed_name in list_dir:
                    print("\n\t=== MAKE BED FILE ===")
                    self.make_bed(self.dict["bam_path"][i])

                bed_file = open(path.splitext(self.dict["bam_path"][i])[0] + ".bed", "r")

                print("\n\t=== READING ===")
                self.x_list = self.read(bed_file)
                # Recover the total number of reads
                reads = self.number_reads(path.splitext(self.dict["bam_path"][i])[0])

                print("\n\t=== NORMALIZATION ===")
                self.dict["y_list"][i] = self.normalize(self.dict["y_list"][i], reads)

        # Determine the y max for the x axis of the graph
        self.y_max = max(self.dict["y_max_list"])
        # Create a title for the graph if is not given
        if not self.title_write:
            self.title = "Coverage along the reference {}".format(self.reference_name)

        # If a max position is not parse, determine the max position for the x axis of the graph
        if not self.max_position_write:
            self.max_position = max(self.dict["x_max_list"])

        print("\n\t=== CREATION OF THE GRAPH ===")
        # Make the graph
        self.coverage()

        self.total_time = round(time() - self.start_time, 1)
        self.make_report()

        print("\n\n====== DONE ======")
        print("\tTotal execution: {}s".format(self.total_time))

        # ~~~~~~~PUBLIC METHODS~~~~~~~#

    def reference(self):
        """
        Recover the name and the length of the reference.
        """
        with open(self.reference_path) as f:
            for ref in SeqIO.parse(f, "fasta"):
                self.reference_name = ref.id
                self.reference_length = len(ref.seq)

    def make_bed(self, bam_path):
        """"
        Make a bed file from bam path in argument.
        """

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Generate index file if not present. Index file is required by pileup.
            if bam.has_index():
                print("\t\tBam index file is present.")
            else:
                print("\t\tCreate a bam index file...")
                pysam.index(bam_path)
                bai = bam_path + ".bai"
                self.created_files.append(bai)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            print("\t\tCreate a bed file...")
            bed_name = path.splitext(bam_path)[0] + ".bed"
            # Create a bed file with the same name of the bam file
            with open(bed_name, "w+") as bed_file:

                self.created_files.append(bed_name)
                for seq_dictionnary in bam.header['SQ']:
                    read_exist = False
                    # Mandatory : The file must contain only alignments against one reference
                    for read in bam.fetch(reference=self.reference_name):
                        if read:
                            read_exist = True
                            break

                # Make a bed file using only the reads aligned to a specific reference
                if read_exist:
                    print("\t\t... using reads aligned to the reference {}.".format(self.reference_name))
                    pos_prec = -1
                    for pileupcolumn in bam.pileup(self.reference_name, max_depth=100000000):

                        # pileupcolumn.pos + 1 == base position in the reference
                        base_position = pileupcolumn.pos + 1

                        # Write the missing positions since pileup do not report uncovered regions
                        if pileupcolumn.pos != pos_prec + 1:
                            # Note: i is a 0-based coordinate
                            for i in range(pos_prec + 1, pileupcolumn.pos):
                                bed_file.write("{}\t{}\t{}\t{}\n".format(self.reference_name, i, i + 1, 0))

                        # Then, write the current position
                        bed_file.write("{}\t{}\t{}\t{}\n".format(self.reference_name, base_position - 1, base_position,
                                                                 pileupcolumn.n))
                        pos_prec = pileupcolumn.pos

                    # Write the last entries for the sequence if needed
                    if pos_prec + 1 < self.reference_length:
                        for i in range(pos_prec + 1, self.reference_length - 1):
                            bed_file.write("{}\t{}\t{}\t{}\n".format(self.reference_name, i, i + 1, 0))

    def number_reads(self, bam_path):
        """
        Determine the number of the reads for a sample.
        """
        reads = 0
        y_list = list()
        # with bam_path as bed:
        with open(bam_path + ".bed", "r") as bed:

            for line in bed:
                y = line.split()[3]
                y_list.append(y)

            # Transform the list content in int
            y_list = list(map(int, y_list))

            for i in y_list:
                reads += i

            return reads

    def read(self, file_name):
        """
        Read bam file done in argument and return lists of positions and coverages.
        """
        x_list = list()
        y_list = list()

        # If zoom is needed, read only the depth of the position of zoom
        if self.max_position_write:
            line_list = file_name.readlines()

            # Add -1 at because list start at O and bed file start at 1
            # No -1 for max_position because range do not include the stop
            for i in range(self.min_position - 1, self.max_position):
                # Recovery of the position
                x = line_list[i].split()[2]
                x_list.append(x)
                # Recovery of the depth
                y = line_list[i].split()[3]
                y_list.append(y)

        # If not zoom, read all of the file
        else:
            for line in file_name:
                # Recovery of the position
                x = line.split()[2]
                x_list.append(x)
                # Recovery of the depth
                y = line.split()[3]
                y_list.append(y)

        # Transform the list content in int
        x_list = list(map(int, x_list))
        y_list = list(map(int, y_list))

        # Add the maximum of the x list in the dictionnary to configure the x axe of the graph
        self.dict["x_max_list"].append(max(x_list))

        # Add x and y list to the global dictionnary
        self.dict["x_list"].append(x_list)
        self.dict["y_list"].append(y_list)

    def normalize(self, depth_list, reads):
        """
        Normalization of the depth by division by the total number of sequenced nucleotides.
        """

        # At each position, the depth is divided by the total number of sequenced nucleotides.
        depth_list[:] = [x / reads for x in depth_list]
        self.dict["y_max_list"].append(max(depth_list))

        return depth_list

    def coverage(self):
        """
        Construction of the coverage graph.
        """
        plt.figure(figsize=(self.length, self.width))

        # For each sample, execute
        for i in range(0, 10):
            # If x list exist,
            try:
                self.dict["x_list"][i]
            # If there is no more sample, break
            except IndexError:
                break
            # If the sample exist, create the graph
            else:
                if self.scale == "log":
                    # Make curve for all samples            
                    plt.semilogy(self.dict["x_list"][i],
                                 self.dict["y_list"][i],
                                 linewidth=self.dict["line_width"][i],
                                 label=self.dict["sample_name"][i],
                                 color=self.dict["line_color"][i])

                    # Add the axis with logarithmic scale on the y axe
                    plt.axis([self.min_position, self.max_position, 1e-10, 2 * self.y_max])

                if self.scale == "linear":
                    # Make curve for all samples
                    plt.plot(self.dict["x_list"][i],
                             self.dict["y_list"][i],
                             linewidth=self.dict["line_width"][i],
                             label=self.dict["sample_name"][i],
                             color=self.dict["line_color"][i])

                    # Add the axis with linear scale on the y axe
                    plt.axis([self.min_position, self.max_position, 0, 1.2 * self.y_max])

        # Add titles
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)

        plt.legend()
        print("\t\tCreate the '{}.svg' file...".format(self.output_name))

        # Save the graph in SVG file (always)
        plt.savefig(self.output_name + ".svg", dpi=200)
        self.created_files.append(self.output_name + ".svg")

        # If a supplementary out is choose
        if self.supplementary_output:
            supplementary_output_name = self.supplementary_output_prefix + "." + self.supplementary_output_format
            print("\t\tCreate the '{}' file...".format(supplementary_output_name))
            self.created_files.append(supplementary_output_name)
            plt.savefig(supplementary_output_name)

    def make_report(self):
        """
        Make a report file with the files created during the run.
        """
        name = "report.txt"
        with open(name, 'w') as report_file:
            report_file.write("\n====== {} ======\n".format(self.VERSION))
            report_file.write("\n{}\n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))

            # Configuration file name
            report_file.write("\n=== Configuration File Name ===\n")
            report_file.write("\t{}\n".format(self.conf_file))
            # Created files
            report_file.write("\n=== Created Files ===\n")
            for file in self.created_files:
                report_file.write("\t{}\n".format(file))

            # Packages versions
            # report_file.write("\n=== Packages versions ===\n")
            ### Maybe try modulefinder package to report the packages versions


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':
    Main.class_init()
