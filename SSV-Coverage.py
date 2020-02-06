#!/usr/bin/env python3
# -*- coding: utf-8 -*

"""
@package    SSV-Coverage
@brief      Main file of the program
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@creation of the first version:
@author     Mathieu Bolteau - 2018
* <mathieu.bolteau1@gmail.com>
* [Github] (https://github.com/mablt)
@new developments
@author     Emilie Lecomte
* <emilie.lecomte@univ-nantes.fr>
* [Github] (https://github.com/emlec)
* [INSERM UMR 1089] (http://umr1089.univ-nantes.fr/)
"""

try:
    # Standard library imports
    import sys  # Mandatory package
    import os  # Mandatory package
    import argparse  # Mandatory package
    import matplotlib  # Mandatory package
    import subprocess  # Required if depth file created using bedtools genomecov
    import pathlib  # Mandatory package
    import gzip  # Mandatory package

    matplotlib.use('svg')  # Mandatory for the hardcopy backend
    import matplotlib.pyplot as plt  # Mandatory package
    from os import path, listdir  # Mandatory package
    import time  # Mandatory package
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
class Main:
    """
    Create a graphic representing the coverage of sequencing along a reference from bam files.
    """
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    # ~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "SSV-Coverage.py 0.2.0"
    USAGE = "SSV-Coverage.py -c Conf.txt [-i -h]"

    # ~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init(self):
        """
        Init class method for instantiation from command line. Parse arguments.
        """

        # Define parser usage, options
        parser = argparse.ArgumentParser(description='Coverage graph creator from BAM files.', usage=self.USAGE)
        parser.add_argument('--version', action='version', version=self.VERSION,
                            help="show program's version number and exit")
        parser.add_argument('-i', dest='bam_number', type=int,
                            help="Generate an empty configuration file adapted to the \
                            number of bam files to create coverage graph and exit. [Mandatory]")
        parser.add_argument('-c', dest='conf_file',
                            help="Path to the configuration file [Mandatory]")

        # Parse arguments
        options = parser.parse_args()

        # If the program is run without option, show the help and exit
        if len(sys.argv[1:]) == 0:
            parser.print_help()
            exit()

        return Main(options.bam_number, options.conf_file)

    # ~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, bam_number=None, conf_file=None):
        """
        Parse the configuration file.
        """

        # Returns the time as a floating point number expressed in seconds
        self.start_time = time.time()

        # Create an empty configuration file if needed with the number of bam files necessary
        if bam_number:
            print(f"Create an empty configuration file, appropriate for {bam_number} bam files, in the current folder.")
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
            "depth_program": list(),
            "max_depth_pysam": int,
            "x_list": list(),
            "y_list": list(),
            "x_max_list": list(),
            "y_max_list": list(),
            "y_max_list_norm": list()
        }

        self.title_write = True
        self.max_position_write = True
        self.reference_name = ""
        self.supplementary_output = False
        self.created_files = list()

        try:
            # Parse [Ref] section
            self.reference_path = self.config.get("Ref", "ref_path")
            if not ".gz" in pathlib.Path(self.reference_path).suffixes:
                print("The reference must be in fa.gz format.")
                exit()
            else:
                if self.reference_path:
                    self.reference()

            # Parse [Bam] section
            for i in range(1, 10):
                bam_num = "Bam" + str(i)
                # If there are no more bam files, stop the data recovery
                if not self.config.has_section(bam_num):
                    break

                self.dict["sample_name"].append(self.config.get(bam_num, "sample_name"))

                self.dict["bam_path"].append(self.config.get(bam_num, "bam_path"))

                # If line color is not indicated, black is used
                if not self.config.get(bam_num, "line_color"):
                    self.dict["line_color"].append("black")
                else:
                    self.dict["line_color"].append(self.config.get(bam_num, "line_color"))

                # If line width is not indicated, 1 is used
                if not self.config.get(bam_num, "line_width"):
                    self.dict["line_width"].append(1)
                else:
                    self.dict["line_width"].append(self.config.get(bam_num, "line_width"))

            # Parse [depth] section
            self.depth_program = self.config.get("Depth", "depth_program")
            if not self.depth_program:
                self.depth_program = "bedtools"

            self.max_depth_pysam = self.config.get("Depth", "max_depth_pysam")
            if not self.max_depth_pysam:
                self.max_depth_pysam = 100000000

            # Parse [Graph] section
            if not self.config.get("Graph", "title"):
                self.title_write = False
            else:
                self.title = self.config.get("Graph", "title")

            if not self.config.get("Graph", "width"):
                self.width = 6
            else:
                self.width = int(self.config.get("Graph", "width"))

            if not self.config.get("Graph", "length"):
                self.length = 8
            else:
                self.length = int(self.config.get("Graph", "length"))

            if self.config.get("Graph", "depth_normalization"):
                self.depth_normalization = self.config.get("Graph", "depth_normalization")
            else:
                self.depth_normalization = 'True'

            self.scale = self.config.get("Graph", "scale")
            if not self.config.get("Graph", "scale"):
                self.scale = "linear"

            if self.scale != "linear":
                if self.scale != "log":
                    print("The 'scale' option of the configuration file is not in the correct format")
                    print("Please report to the descriptions in the configuration file\n")
                    exit()

            self.xlabel = self.config.get("Graph", "xlabel")
            if not self.xlabel:
                self.xlabel = "Position"

            if not self.config.get("Graph", "ylabel"):
                if self.depth_normalization == 'True':
                    self.ylabel = "Normalized coverage of sequencing per base"
                elif self.depth_normalization == 'False':
                    self.ylabel = "Coverage of sequencing per base"
            else:
                self.ylabel = self.config.get("Graph", "ylabel")

            # Parse [Zoom] section

            #   If zoom is required
            if self.config.getboolean("Zoom", "zoom"):

                # If the min position is superior to the max position
                if int(self.config.get("Zoom", "min_position")) > int(self.config.get("Zoom", "max_position")):
                    print("One of the value in the 'Zoom' section of the configuration file is not correct")
                    print("Please report to the descriptions in the configuration file\n")
                    exit()

                # If the min position is inferior to 1
                if int(self.config.get("Zoom", "min_position")) < 1:
                    print("One of the value in the 'Zoom' section of the configuration file is not correct")
                    print("Please report to the descriptions in the configuration file\n")
                    exit()
                else:
                    self.min_position = int(self.config.get("Zoom", "min_position"))

                # If the max position is superior to the length reference
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

            # Parse [Output] section

            if not self.config.get("Output", "output_name"):
                self.output_name = "out"
            else:
                self.output_name = self.config.get("Output", "output_name")

            # If supplementary output file is required
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

        # Workflow
        for i in range(0, 10):
            # Try if the bam_path exist in the configuration file
            try:
                self.dict["bam_path"][i]
            # If not exist, stop the execution on the bam_path
            except IndexError:
                break
            # If exist, execute
            else:
                print(f"\n\n====== SAMPLE {self.dict['sample_name'][i]} ======\n")
                # Recovery of the current directory
                directory = path.dirname(self.dict["bam_path"][i])
                if directory == "":
                    directory = "./"
                # List containing the names of the files in the directory.
                list_dir = listdir(directory)

                # depth file naming from bam file name
                depth_name = path.basename(
                    path.splitext(self.dict["bam_path"][i])[0] + "_" + self.depth_program + ".depth")

                # If depth file doesn't exist, it will be created
                if not depth_name in list_dir:
                    print("\n\t=== MAKE DEPTH FILE ===")
                    self.make_depth(self.dict["bam_path"][i])

                depth_file = open(path.splitext(self.dict["bam_path"][i])[0] + "_" + self.depth_program + ".depth", "r")

                print("\n\t=== READING ===")
                self.x_list = self.read(depth_file)

                # Count the total number of sequenced bases
                all_bases = self.all_bases(
                    path.splitext(self.dict["bam_path"][i])[0] + "_" + self.depth_program + ".depth")

                if self.depth_normalization == "True":
                    print("\n\t=== NORMALIZATION ===")
                    self.dict["y_list"][i] = self.normalize(self.dict["y_list"][i], all_bases)
                    self.dict["y_max_list_norm"].append(max(self.dict["y_list"][i]))

        # Determine the y max for the x axis of the graph
        if self.depth_normalization == "True":
            self.y_max = max(self.dict["y_max_list_norm"])
        else:
            self.y_max = max(self.dict["y_max_list"])
        # Create a title for the graph if is not given
        if not self.title_write:
            self.title = f"Coverage along the reference {self.reference_name}"
        # If a max position is not parse, determine the max position for the x axis of the graph
        if not self.max_position_write:
            self.max_position = max(self.dict["x_max_list"])

        print("\n\t=== CREATION OF THE GRAPH ===")
        # Make the graph
        self.coverage()
        self.total_time = round(time.time() - self.start_time, 1)
        self.make_report()

        print("\n\n====== DONE ======")
        print(f"\tTotal execution: {self.total_time}s")

        # ~~~~~~~PUBLIC METHODS~~~~~~~#

    def reference(self):
        """
        Recover the name and the length of the reference.
        """
        with gzip.open(self.reference_path, "rt") as f:
            for ref in SeqIO.parse(f, "fasta"):
                self.reference_name = ref.id
                self.reference_length = len(ref.seq)

    def make_depth(self, bam_path):
        """"
        Make a depth file from bam path passed as argument.
        """
        # Generate index file if not present. Index file is required by pileup.
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if bam.has_index():
                print("\t\tBam index file is present.")
            else:
                print("\t\tCreate a bam index file...")
                pysam.index(bam_path)
                bai = bam_path + ".bai"
                self.created_files.append(bai)

        # Depth file creation ...
        print(f"\t\tCreate a depth file using {self.depth_program} ...")
        print(f"\t\t... and reads aligned to the reference {self.reference_name}.")
        # Create a depth file with the name of the depth program in addition to the name of the bam file
        depth_name = path.splitext(bam_path)[0] + "_" + self.depth_program + ".depth"
        self.created_files.append(depth_name)

        with open("depth.temp", "a") as depth_temp_file:

            # ... using bedtools
            if self.depth_program == "bedtools":
                # Report the depth at each genome position with 0-based coordinates, excluding zero positions.
                subprocess.run(["bedtools", "genomecov", "-dz", "-ibam", bam_path], stdout=depth_temp_file,
                               check=True)


            # ... using pysam
            elif self.depth_program == "pysam":
                with pysam.AlignmentFile(bam_path, "rb") as bam:
                    # pileupcolumn.pos + 1 == base position in the reference
                    pos_prec = -1
                    for pileupcolumn in bam.pileup(self.reference_name, max_depth=self.max_depth_pysam):
                        # Write the missing positions since pileup do not report uncovered regions
                        if pileupcolumn.pos != pos_prec + 1:
                            # Note: i is a 0-based coordinate
                            for i in range(pos_prec + 1, pileupcolumn.pos):
                                depth_temp_file.write(f"{self.reference_name}\t{i - 1}\t{0}\n")
                        # Then, write the current position
                        depth_temp_file.write(
                            f"{self.reference_name}\t{pileupcolumn.pos}\t{pileupcolumn.n}\n")
                        pos_prec = pileupcolumn.pos
                        # Write the last entries for the sequence if needed
                    if pos_prec + 1 < self.reference_length:
                        for i in range(pos_prec + 1, self.reference_length - 1):
                            depth_temp_file.write(f"{self.reference_name}\t{i - 1}\t{0}\n")

        # Mofify the depth file including a column with 1-based coordinates and adding 0-depth if necessary (bedtools output)
        with open("depth.temp", "r") as depth_temp_file:
            with open(depth_name, "a") as depth_file:
                # current position refers to 0-based coordinates
                current_position = 0
                for line in depth_temp_file:
                    columns = line.split("\t")
                    if columns[1] != str(current_position):
                        print(columns[1])
                        print(current_position)
                        depth_file.write(self.reference_name + "\t" + str(current_position) + "\t" + str(
                            current_position + 1) + "\t" + "0" + "\n")
                        current_position = current_position + 1
                    if columns[1] == str(current_position):
                        if current_position < self.reference_length:
                            # Convert the type of the depth number from string to float with 0 digit, because of exponential number incompatibility
                            # that can be present in the bedtools output
                            columns[2] = f'{float(columns[2]):.0f}'
                            # Convert the type of the depth number from float to string in order to include the line return
                            columns[2] = str(columns[2])+"\n"
                            columns.insert(2, str(current_position + 1))
                            depth_file.write("\t".join(columns))
                            current_position = current_position + 1

        os.remove("depth.temp")

        # Check the depth file
        number_of_lines = len(open(depth_name).readlines())
        if number_of_lines == self.reference_length:
            print("\t\tThe depth file is created with success.")
        else:
            print("\t\tBe careful, the depth file is truncated.")
            exit()

    def all_bases(self, depth_path):
        """
        Determine the number of reads for a sample.
        """
        all_bases = 0
        y_list = list()
        # with bam_path as depth:
        with open(depth_path, "r") as depth:
            for line in depth:
                y = line.split()[3]
                y_list.append(y)

            # Transform the list content in int
            y_list = list(map(int, y_list))

            for i in y_list:
                all_bases += i

            return all_bases

    def read(self, file_name):
        """
        Read bam file done in argument and return the lists of positions and coverages.
        """
        x_list = list()
        y_list = list()

        # If zoom is needed, read only the depth of the position of zoom
        if self.max_position_write:
            line_list = file_name.readlines()

            # Add -1 because list start at O and depth file start at 1
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
        self.dict["y_max_list"].append(max(y_list))

    def normalize(self, depth_list, all_bases):
        """
        Normalization of the depth by division by the total number of sequenced nucleotides.
        """

        # At each position, the depth is divided by the total number of sequenced nucleotides.
        depth_list[:] = [base_number / all_bases for base_number in depth_list]
        return depth_list

    def coverage(self):
        """
        Creation of the graphic of coverage.
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
        print(f"\t\tCreate the '{self.output_name}.svg' file...")

        # Save the graph in SVG file (always)
        plt.savefig(self.output_name + ".svg", dpi=200)
        self.created_files.append(self.output_name + ".svg")

        # If a supplementary out is choose
        if self.supplementary_output:
            supplementary_output_name = self.supplementary_output_prefix + "." + self.supplementary_output_format
            print(f"\t\tCreate the '{supplementary_output_name}' file...")
            self.created_files.append(supplementary_output_name)
            plt.savefig(supplementary_output_name)

    def make_report(self):
        """
        Make a report file with the files created during the run.
        """
        name = "SSV-Coverage_report_" + self.depth_program + ".txt"
        with open(name, 'w') as report_file:
            report_file.write(f"\n====== {self.VERSION} ======\n")
            report_file.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

            # Configuration file name
            report_file.write("\n=== Configuration File Name ===\n")
            report_file.write(f"\t{self.conf_file}\n")
            # Created files
            report_file.write("\n=== Created Files ===\n")
            for file in self.created_files:
                report_file.write(f"\t{file}\n")

            # Packages versions
            # report_file.write("\n=== Packages versions ===\n")
            ### Maybe try modulefinder package to report the packages versions


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':
    Main.class_init()
