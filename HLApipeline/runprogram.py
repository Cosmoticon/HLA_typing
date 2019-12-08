import os
import subprocess
import sys
import exceptions
import pipeline


def HLA_LA(hla_pipeline, path_file):
    """Runs the HLA-LA algorithm, and indexes/converts the input file if necessary"""
    if path_file.endswith("cram"):
        bam_file = CRAM_to_BAM(hla_pipeline, path_file)
    else:
        bam_file = hla_pipeline.path_WGS
    if not os.path.isfile(bam_file.replace(".bam", ".bai")):
        SAM_index(hla_pipeline, bam_file)
    HLA_LA_script = ("{HLA_LA_dir}HLA-LA.pl "
                    "--BAM {bam} "
                    "--graph PRG_MHC_GRCh38_withIMGT "
                    "--sampleID {id} "
                    "--maxThreads {threads} "
                    "--picard_sam2fastq_bin {picard_dir}picard.jar "
                    "--workingDir {out}HLA_LA/").format(HLA_LA_dir=hla_pipeline.dir_HLA_LA,
                                                bam=bam_file,
                                                id=hla_pipeline.string_name,
                                                threads=hla_pipeline.num_threads,
                                                picard_dir=hla_pipeline.dir_picard,
                                                out=hla_pipeline.output_dir)
    hla_pipeline.write_log("Running HLA-LA for {file}".format(file=bam_file))
    HLA_LA_sub = subprocess.Popen(HLA_LA_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = HLA_LA_sub.communicate()
    print(outs)
    print(errs)
    if HLA_LA_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()

def HLA_HD(hla_pipeline):
    """Runs the HLA-HD algorithm"""
    HLA_HD_script = ("{HLA_HD_dir}/hlahd.sh "
                    "-t {threads} -m {min_len} "
                    "-c 0.95 -f {HLA_HD_dir}freq_data "
                    "{fastq1} {fastq2} {HLA_HD_dir}HLA_gene.split.txt {HLA_HD_dir}dictionary "
                    "{id} {out}HLA_HD").format(HLA_HD_dir=hla_pipeline.dir_HLA_HD,
                                                min_len=hla_pipeline.RNAseq_min_length,
                                                id=hla_pipeline.string_name,
                                                fastq1=hla_pipeline.path_RNAseq1,
                                                fastq2=hla_pipeline.path_RNAseq2, 
                                                threads=hla_pipeline.num_threads,
                                                out=hla_pipeline.output_dir)
    hla_pipeline.write_log("Running HLA-HD for {file}".format(id=hla_pipeline.path_RNAseq1+","+hla_pipeline.path_RNAseq2))
    HLA_HD_sub = subprocess.Popen(HLA_HD_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = HLA_HD_sub.communicate()
    print(outs)
    print(errs)
    if HLA_HD_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()

def seq2HLA(hla_pipeline):
    """Runs the seq2HLA algorithm"""
    seq2HLA_script = ("{seq2HLA_dir}/seq2HLA.py "
                        "-1 {fastq1} -2 {fastq2} "
                        "-r {out}seq2HLA/{id} -p {threads} ").format(seq2HLA_dir=hla_pipeline.dir_seq2HLA,
                                                fastq1=hla_pipeline.path_RNAseq1,
                                                fastq2=hla_pipeline.path_RNAseq2, 
                                                threads=hla_pipeline.num_threads,
                                                out=hla_pipeline.output_dir,
                                                id=hla_pipeline.string_name)
    hla_pipeline.write_log("Running seq2HLA for {file}".format(id=hla_pipeline.path_RNAseq1+","+hla_pipeline.path_RNAseq2))
    seq2HLA_sub = subprocess.Popen(seq2HLA_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = seq2HLA_sub.communicate()
    print(outs)
    print(errs)
    if seq2HLA_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()

def arcasHLA(hla_pipeline):
    """
    Runs the arcasHLA algorithm, first extracting the relevant sequences 
    from the bam file, and then running the HLA typing algorithm
    """
    bamfile = RNAseq_align(hla_pipeline)
    arcasHLA_extract_script = ("arcasHLA extract --paired True --unmapped True "
                                "--outdir {out}arcasHLA/ "
                                "--threads {threads} {bam} ").format(bam=bamfile, 
                                                threads=hla_pipeline.num_threads,
                                                out=hla_pipeline.output_dir)
    hla_pipeline.write_log("Running arcasHLA extract for {file}".format(id=hla_pipeline.path_RNAseq1+","+hla_pipeline.path_RNAseq2))
    arcasHLA_extract_sub = subprocess.Popen(arcasHLA_extract_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = arcasHLA_extract_sub.communicate()
    print(outs)
    print(errs)
    if arcasHLA_extract_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()
    arcasHLA_genotype_script = ("arcasHLA genotype --genes A,B,C,DQB1,DRB1 --population {pop} "
                                "--min_count 75 --tolerance 10e-7 --max_iterations 1000 "
                                "--outdir {out}arcasHLA/ --threads {threads} "
                                "{out}arcasHLA/{id}Aligned.out.1.fq.gz  {out}arcasHLA/{id}Aligned.out.2.fq.gz ").format(pop=hla_pipeline.population, 
                                                                                                    threads=hla_pipeline.num_threads,
                                                                                                    out=hla_pipeline.output_dir,
                                                                                                    id=hla_pipeline.string_name)
    hla_pipeline.write_log("Running arcasHLA genotype for {file}".format(id=hla_pipeline.path_RNAseq1+","+hla_pipeline.path_RNAseq2))
    arcasHLA_genotype_sub = subprocess.Popen(arcasHLA_genotype_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = arcasHLA_genotype_sub.communicate()
    print(outs)
    print(errs)
    if arcasHLA_genotype_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()

def Kourami(hla_pipeline, path_file):
    """
    Runs the Kourami algorithm, first extracting the relevant sequences 
    from the bam file, and then running the HLA typing algorithm
    """
    if path_file.endswith("cram"):
        bam_file = CRAM_to_BAM(hla_pipeline.path_WGS, path_file)
    else:
        bam_file = path_file
    bam_file = RNAseq_align(hla_pipeline)
    Kourami_extract_script = ("{kourami}/scripts/alignAndExtract_hs38DH.sh {id} {bam} ").format(bam=bam_file, 
                                                                                            kourami=hla_pipeline.dir_kourami,
                                                                                            id=hla_pipeline.string_name)
    hla_pipeline.write_log("Running Kourami extract for {file}".format(id=bam_file))
    Kourami_extract_sub = subprocess.Popen(Kourami_extract_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = Kourami_extract_sub.communicate()
    print(outs) 
    print(errs)
    if Kourami_extract_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()
    Kourami_run_script = ("java -jar {kourami}/Kourami.jar --msaDirectory {kourami}resources/hs38NoAltDH.fa "
                            "--outputfilePrefix {out}Kourami/{id} {bam_dir}{id}_on_KouramiPanel.bam").format(kourami=hla_pipeline.dir_Kourami, 
                                                                                                    bam_dir=bam_file.rsplit(",", 1)[0]+"/",
                                                                                                    id=hla_pipeline.string_name,
                                                                                                    out=hla_pipeline.output_dir)
    hla_pipeline.write_log("Running Kourami for {file}".format(id=bam_file))
    Kourami_run_sub = subprocess.Popen(Kourami_run_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = Kourami_run_sub.communicate()
    print(outs)
    print(errs)
    if Kourami_run_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()

def SAM_index(hla_pipeline, bam_file):
    """Indexes the WGS/WES bam file"""
    SAM_ind_script = ("samtools index "
                    "{bam}").format(bam=bam_file)
    hla_pipeline.write_log("Generating index file for {f}".format(f=bam_file))
    SAM_ind_sub = subprocess.Popen(SAM_ind_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = SAM_ind_sub.communicate()
    print(outs)
    print(errs)
    if SAM_ind_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()

def CRAM_to_BAM(hla_pipeline, cram_file):
    """Converts .cram to .bam file conversion"""
    SAM_view_script = ("samtools view -b -o {output}"
                    "{cram}").format(output=cram_file.replace(".cram", ".bam"),
                        cram=cram_file)
    hla_pipeline.write_log("Converting cram to bam for {f}".format(f=cram_file))
    SAM_view_sub = subprocess.Popen(SAM_view_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = SAM_view_sub.communicate()
    print(outs)
    print(errs)
    if SAM_view_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()
    return cram_file.replace(".cram", ".bam")

def RNAseq_align(hla_pipeline):
    """Runs the STAR alignment of RNAseq data to the reference genome"""
    align_prefix = hla_pipeline.path_RNAseq1.rsplit('/', 1)[0] + "/" + hla_pipeline.string_name
    if not os.path.isfile(align_prefix + "Aligned.out.bam"):
        STAR_align_script = ("STAR --genomeDir {g_dir} "
                            "--readFilesIn {r1} {r2} "
                            "--outSAMstrandField intronMotif "
                            "--outSAMattrIHstart 0 "
                            "--alignSoftClipAtReferenceEnds No "
                            "--outFilterIntronMotifs RemoveNoncanonical "
                            "--outSAMtype BAM SortedByCoordinate "
                            "--outFileNamePrefix {align_dir} "
                            "--alignEndsType EndToEnd "
                            "--outSAMunmapped Within").format(g_dir=hla_pipeline.dir_genome,
                                r1=hla_pipeline.path_RNAseq1,
                                r2=hla_pipeline.path_RNAseq2,
                                align_dir=align_prefix)
        hla_pipeline.write_log("Alignment of {prefix}".format(prefix=align_prefix))
        STAR_align_sub = subprocess.Popen(STAR_align_script.split(),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        universal_newlines=True)
        outs, errs = STAR_align_sub.communicate()
        print(outs)
        print(errs)
        if STAR_align_sub.returncode != 0:
            print("An error occurred, the program did not run to completion.")
            sys.exit()
    return align_prefix + "Aligned.out.bam"