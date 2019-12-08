from datetime import datetime
import exceptions
import pandas as pd
import os

class HLAPipeline:
    """This class keeps track of all the files/parameters for this complex pipeline"""
    
    def __init__(self):
        """Initializes HLA pipeline class and sets the default number of threads to 1"""
        self.num_threads = "1"

    def time_stamp(self):
        now = datetime.now()
        return(now.strftime("%m/%d/%Y, %H:%M:%S"))

    def write_log(self, message):
        """Writes message to output log in case user wants to see the steps executed"""
        with open(self.output_dir + "log.txt", "a+") as log:
            log.write(self.time_stamp() + ":  " + message + "\n")

    def string_name(self):
        """Returns the name of the output directory, which should be the sample name"""
        return self.output_dir.split("/")[-2]

    def set_output_dir(self, output_dir):
        """Sets path to output directory according to user input, throws error if it cannot make the directory"""
        try:
            os.makedirs(output_dir, exist_ok=True)
            if output_dir[-1] != "/":
                output_dir += "/"
            self.output_dir = output_dir
        except:
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=output_dir)
            raise exceptions.IncorrectPathError(message)
    
    def set_dir_genome(self, dir_genome):
        """Sets path to genome directory according to user input"""
        if not os.path.isdir(dir_genome):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=dir_genome)
            raise exceptions.IncorrectPathError(message)
        else:
            if dir_genome[-1] != "/":
                dir_genome += "/"
            self.dir_genome = dir_genome

    def set_num_threads(self, num_threads):
        """Sets number of threads according to user input"""
        try:
            int_min = int(num_threads)
            if int_min > 0:
                self.num_threads = num_threads
            else:
                message = ("Error: the following value for number of threads <{num_threads}> "
                    "is not correct").format(num_threads=num_threads)
                raise exceptions.WrongArgumentError(message)
        except:
            message = ("Error: the following value for number of threads <{num_threads}> "
                    "is not correct").format(num_threads=num_threads)
            raise exceptions.WrongArgumentError(message)
    
    def set_path_correct_HLA(self, correct_HLA):
        """Sets path to correct HLA typing CSV file according to user input"""
        if not os.path.isfile(correct_HLA) or not correct_HLA.endswith(".csv"):
            message = ("Error: the following file path <{path}> "
                    "is incorrect").format(path=correct_HLA)
            raise exceptions.IncorrectPathError(message)
        else:
            self.path_correct_HLA = correct_HLA

    def set_path_WGS(self, input_WGS):
        """Sets path to WGS file according to user input"""
        if not os.path.isfile(input_WGS):
            message = ("Error: the following file path <{path}> "
                    "is incorrect").format(path=input_WGS)
            raise exceptions.IncorrectPathError(message)
        else:
            self.path_WGS = input_WGS
    
    def set_path_WES(self, input_WES):
        """Sets path to WES file according to user input"""
        if not os.path.isfile(input_WES):
            message = ("Error: the following file path <{path}> "
                    "is incorrect").format(path=input_WES)
            raise exceptions.IncorrectPathError(message)
        else:
            self.path_WES = input_WES
    
    def set_path_RNAseq(self, input_RNAseq):
        """Sets path to RNAseq file according to user input"""
        list_RNAseq = input_RNAseq.split(",")
        if len(list_RNAseq) != 2:
            message = ("Error: wrong number of RNAseq files").format(path=input_RNAseq)
            raise exceptions.IncorrectPathError(message)
        else:
            for path_RNA in list_RNAseq:
                if not os.path.isfile(path_RNA):
                    message = ("Error: the following file path <{path}> "
                            "is incorrect").format(path=path_RNA)
                    raise exceptions.IncorrectPathError(message)
            self.path_RNAseq1 = list_RNAseq[0]
            self.path_RNAseq2 = list_RNAseq[1]
    
    def set_RNAseq_min_length(self, RNAseq_min_length):
        """Sets minimum length for RNAseq sequences according to user input"""
        try:
            int_min = int(RNAseq_min_length)
            if int_min > 0:
                self.RNAseq_min_length = RNAseq_min_length
            else:
                message = ("Error: the following value for minimum RNAseq length <{min_len}> "
                    "is not correct").format(min_len=RNAseq_min_length)
                raise exceptions.WrongArgumentError(message)
        except:
            message = ("Error: the following value for minimum RNAseq length <{min_len}> "
                    "is not correct").format(min_len=RNAseq_min_length)
            raise exceptions.WrongArgumentError(message)
    
    def set_population(self, population):
        """Sets population origin of sample according to user input"""
        if population in ["asian_pacific_islander", "black", "caucasian", "hispanic", "native_american"]:
            self.population = population
        else:
            message = ("Error: the following population <{population}> "
                    "is not correct").format(population=population)
            raise exceptions.WrongArgumentError(message)
    
    def set_dir_picard(self, dir_picard):
        """Sets path to picard directory according to user input"""
        if not os.path.isdir(dir_picard):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=dir_picard)
            raise exceptions.IncorrectPathError(message)
        else:
            if dir_picard[-1] != "/":
                dir_picard += "/"
            self.dir_picard = dir_picard

    def set_dir_HLA_HD(self, dir_HLA_HD):
        """Sets path to HLA-HD directory according to user input"""
        if not os.path.isdir(dir_HLA_HD):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=dir_HLA_HD)
            raise exceptions.IncorrectPathError(message)
        else:
            if dir_HLA_HD[-1] != "/":
                dir_HLA_HD += "/"
            self.dir_HLA_HD = dir_HLA_HD
    
    def set_dir_seq2HLA(self, dir_seq2HLA):
        """Sets path to seq2HLA directory according to user input"""
        if not os.path.isdir(dir_seq2HLA):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=dir_seq2HLA)
            raise exceptions.IncorrectPathError(message)
        else:
            if dir_seq2HLA[-1] != "/":
                dir_seq2HLA += "/"
            self.dir_seq2HLA = dir_seq2HLA
    
    def set_dir_Kourami(self, dir_Kourami):
        """Sets path to Kourami directory according to user input"""
        if not os.path.isdir(dir_Kourami):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=dir_Kourami)
            raise exceptions.IncorrectPathError(message)
        else:
            if dir_Kourami[-1] != "/":
                dir_Kourami += "/"
            self.dir_Kourami = dir_Kourami
    
    def set_dir_arcasHLA(self, dir_arcasHLA):
        """Sets path to arcasHLA directory according to user input"""
        if not os.path.isdir(dir_arcasHLA):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=dir_arcasHLA)
            raise exceptions.IncorrectPathError(message)
        else:
            if dir_arcasHLA[-1] != "/":
                dir_arcasHLA += "/"
            self.dir_arcasHLA = dir_arcasHLA
    
    def set_dir_HLA_LA(self, dir_HLA_LA):
        """Sets path to HLA-LA directory according to user input"""
        if not os.path.isdir(dir_HLA_LA):
            message = ("Error: the following directory path <{path}> "
                    "is not correct").format(path=dir_HLA_LA)
            raise exceptions.IncorrectPathError(message)
        else:
            if dir_HLA_LA[-1] != "/":
                dir_HLA_LA += "/"
            self.dir_HLA_LA = dir_HLA_LA