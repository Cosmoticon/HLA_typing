import sys
import getopt
import exceptions
import pipeline
import analysis
import runprogram

def main():
    """
    The main function: 
    - parses through all the command line arguments
    - creates the hla_pipeline class
    - verifies there are no missing arguments
    -runs the program and analysis
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                "g:e:r:G:l:p:o:t:H:s:K:a:P:L:c:h",
                                ["input_WGS=",
                                    "input_WES=",
                                    "input_RNAseq=",
                                    "genome_dir=",
                                    "RNAseq_min_length=",
                                    "population=",
                                    "output_dir=",
                                    "num_threads=",
                                    "dir_HLA_HD=",
                                    "dir_seq2HLA=",
                                    "dir_Kourami=",
                                    "dir_arcasHLA",
                                    "dir_picard=",
                                    "dir_HLA_LA=",
                                    "correct_HLA",
                                    "help"
                                    ]
                                )
    except getopt.GetoptError as e:
        print(e)
        sys.exit(2)
    if len(args) > 0:
        message = "Error: non-paired arguments are not allowed."
        raise exceptions.WrongArgumentError(message)
    
    HLA_pipeline = pipeline.HLAPipeline()
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            description()
            sys.exit()
        elif opt in ("-g", "--input_WGS"):
            HLA_pipeline.set_path_WGS(arg)
        elif opt in ("-e", "--input_WES"):
            HLA_pipeline.set_path_WES(arg)
        elif opt in ("-r", "--input_RNAseq"):
            HLA_pipeline.set_path_RNAseq(arg)
        elif opt in ("-G", "--dir_genome"):
            HLA_pipeline.set_dir_genome(arg)
        elif opt in ("-l", "--RNAseq_min_length"):
            HLA_pipeline.set_RNAseq_min_length(arg)
        elif opt in ("-p", "--population"):
            HLA_pipeline.set_population(arg)
        elif opt in ("-o", "--output_dir"):
            HLA_pipeline.set_output_dir(arg)
        elif opt in ("-t", "--num_threads"):
            HLA_pipeline.set_num_threads(arg)
        elif opt in ("-H", "--dir_HLA_HD"):
            HLA_pipeline.set_dir_HLA_HD(arg)
        elif opt in ("-s", "--dir_seq2HLA"):
            HLA_pipeline.set_dir_seq2HLA(arg)
        elif opt in ("-K", "--dir_Kourami"):
            HLA_pipeline.set_dir_Kourami(arg)
        elif opt in ("-a", "--dir_arcasHLA"):
            HLA_pipeline.set_dir_arcasHLA(arg)
        elif opt in ("-P", "--dir_picard"):
            HLA_pipeline.set_dir_picard(arg)
        elif opt in ("-L", "--dir_HLA_LA"):
            HLA_pipeline.set_dir_HLA_LA(arg)
        elif opt in ("-c", "--correct_HLA"):
            HLA_pipeline.set_path_correct_HLA(arg)
        else:
            message = "Error: {opt} is not a valid option".format(opt=opt)
            raise exceptions.WrongArgumentError(message)
    for pipeline_attr in ["input_WGS",
                            "input_WES",
                            "input_RNAseq",
                            "RNAseq_min_length",
                            "population",
                            "output_dir",
                            "dir_genome",
                            "dir_HLA_HD",
                            "dir_arcasHLA",
                            "dir_seq2HLA",
                            "dir_Kourami",
                            "dir_picard",
                            "dir_HLA_LA"]:
        if not hasattr(HLA_pipeline, pipeline_attr):
            message = ("Error: you must indicate --{attr}.").format(attr=pipeline_attr)
            raise exceptions.MissingArgumentError(message)
    
    for path in (HLA_pipeline.path_WGS, HLA_pipeline.path_WES):
        runprogram.Kourami(HLA_pipeline, path)
        runprogram.HLA_LA(HLA_pipeline, path)
    runprogram.seq2HLA(HLA_pipeline)
    runprogram.HLA_HD(HLA_pipeline)
    runprogram.arcasHLA(HLA_pipeline)
    analysis.extract_results(HLA_pipeline)
    analysis.calculate_accuracy(HLA_pipeline)

def description():
    manual = ("\nOPTIONS:\n\t -h or --help : display the manual\n"
        "\t -g or --input_WGS : specify the path to the whole genome sequence file\n"
        "\t -e or --input_WES : specify the path to the whole exome sequence file \n"
        "\t -r or --input_RNAseq : specify the path to both RNAseq sequence files (e.g. path/to/RNAseq1,path/to/RNAseq2)\n"
        "\t -l or --RNAseq_min_length : specify RNAseq sequence length threshold for mapping\n"
        "\t -G or --genome_dir : specify the path to directory containing all the genome files\n"
        "\t -p or --population : specify the population origin of the sample\n"
        "\t -o or --output_dir : specify the output directory for all the algorithms\n"
        "\t -t or --num_threads : specify the number of threads to be used for parallel computing\n"
        "\t -H or --dir_HLA_HD : specify the directory containing all the HLA-HD program files\n"
        "\t -s or --dir_seq2HLA : specify the directory containing all the seq2HLA program files\n"
        "\t -K or --dir_Kourami : specify the directory containing all the Kourami program files\n"
        "\t -a or --dir_arcasHLA : specify the directory containing all the arcasHLA program files\n"
        "\t -L or --dir_HLA_LA : specify the directory containing all the HLA-LA program files\n"
        "\t -P or --dir_picard : specify the directory containing all the picard program files\n"
        "\t -c or --correct_HLA: specify the CSV file containing the correct HLA typing information from PCR-Sanger sequencing\n")
    print(manual)

if __name__ == "__main__":
    main()