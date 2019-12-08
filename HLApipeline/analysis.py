import pandas as pd
import pipeline
import json


def extract_results(hla_pipeline):
    """Extracts results from all the algorithms, trims the alleles, and outputs into CSV file"""
    dict_HLA_HD = {"algorithm": "HLA_HD"}
    with open("{out}{id}/result/{id}_final.result.txt", "r") as HLA_HD:
        results_HLA_HD = HLA_HD.readlines()
        for result in results_HLA_HD:
            if len(result.split()) == 3:
                if result.split()[0] in ["A", "B", "C", "DQB1", "DRB1"]:
                    dict_HLA_HD[result.split()[0] + "-1"] = result.split()[1].rsplit(":", 1)[0]
                    dict_HLA_HD[result.split()[0] + "-2"] = result.split()[2].rsplit(":", 1)[0]
    dict_arcasHLA = {"algorithm": "arcasHLA"}
    with open(hla_pipeline.output_dir + "arcasHLA/" + 
                hla_pipeline.string_name + "Aligned.out.genotype.json", "r") as json_file:
        dict_json = json.load(json_file)
        for allele in dict_json:
            dict_arcasHLA[allele + "-1"] = dict_json[allele][0]
            dict_arcasHLA[allele + "-2"] = dict_json[allele][1]
    dict_seq2HLA = {"algorithm": "seq2HLA"}
    with open(("{out}seq2HLA/{id}-ClassII.HLAgenotype4digits").format(out=hla_pipeline.output_dir,
                                                                        id=hla_pipeline.string_name), "r") as seq2HLA_class2:
        results_seq2HLA_class2 = seq2HLA_class2.readlines()
        for result in results_seq2HLA_class2:
            if len(result.split()) == 3:
                if result.split()[0] in ["DQB", "DRB"]:
                    dict_seq2HLA[result.split()[0] + "1-1"] = result.split()[1]
                    dict_seq2HLA[result.split()[0] + "1-2"] = result.split()[2]
    with open(("{out}seq2HLA/{id}-ClassI.HLAgenotype4digits").format(out=hla_pipeline.output_dir,
                                                                        id=hla_pipeline.string_name), "r") as seq2HLA_class1:
        results_seq2HLA_class1 = seq2HLA_class1.readlines()
        for result in results_seq2HLA_class1:
            if len(result.split()) == 3:
                if result.split()[0] in ["A", "B", "C"]:
                    dict_seq2HLA[result.split()[0] + "-1"] = result.split()[1]
                    dict_seq2HLA[result.split()[0] + "-2"] = result.split()[2]
    dict_kourami = {"algorithm": "Kourami"}
    with open(hla_pipeline.output_dir + "Kourami/" + hla_pipeline.string_name + ".result") as kourami:
        results_kourami = kourami.readlines()
        for result in results_kourami:
            if result.split()[0].split("*")[0] in ["A", "B", "C", "DQB1", "DRB1"]:
                if result.split()[0].split("*")[0] + "-1" in dict_kourami:
                    dict_kourami[result.split()[0].split("*")[0] + "-2"] = result.split()[0].rsplit(":", 1)[0]
                else:
                    dict_kourami[result.split()[0].split("*")[0] + "-1"] = result.split()[0].rsplit(":", 1)[0]
    dict_HLA_LA = {"algorithm": "HLA_LA"}
    with open(hla_pipeline.output_dir + "HLA_LA/" + hla_pipeline.string_name + "/hla/R1_bestguess_G.txt", "r") as HLA_LA:
        results_HLA_LA = HLA_LA.readlines()
        for result in results_HLA_LA:
            if result.split()[0].split("-")[1] in ["A", "B", "C", "DQB1", "DRB1"]:
                dict_HLA_LA[result.split()[0].split("-")[1] + "-" + result.split()[1]] = result.split()[2].rsplit(":", 1)[0]
    total_results = pd.DataFrame([dict_HLA_LA, dict_HLA_HD, dict_kourami, dict_seq2HLA, dict_arcasHLA])
    total_results.to_csv(hla_pipeline.output_dir + hla_pipeline.string_name + "_total_results.csv")

def calculate_accuracy(hla_pipeline):
    """Calculates the number of correct and incorrect alleles for each algorithm using the 1000 Genomes data"""
    correct_HLA = pd.read_csv(hla_pipeline.path_correct_hla)
    results = pd.read_csv(hla_pipeline.output_dir + "total_results.csv")
    correct_values = correct_HLA.loc[correct_HLA['id'] == hla_pipeline.string_name]
    list_accuracy_dict = []
    for row in results.iterrows():
        count_incorrect = 0
        count_correct = 0
        for column in ["A", "B", "C", "DQB1", "DRB1"]:
            if results[column + "-1"].split("*")[1] == correct_values[column]:
                if results[column + "-2"].split("*")[1] == correct_values[column + ".1"]:
                    count_correct += 2
                else:
                    count_correct += 1
            elif results[column + "-1"].split("*")[1] == correct_values[column + ".1"]:
                if results[column + "-2"].split("*")[1] == correct_values[column]:
                    count_correct += 2
                else:
                    count_correct += 1
            elif (results[column + "-2"].split("*")[1] == correct_values[column] or 
                results[column + "-2"].split("*")[1] == correct_values[column + ".1"]):
                count_correct += 1
        list_accuracy_dict.append({"algorithm": row["agorithm"], "incorrect": count_incorrect, "correct": count_correct})
    accuracy_results = pd.DataFramce(list_accuracy_dict)
    accuracy_results.to_csv(hla_pipeline.output_dir + hla_pipeline.string_name + "_accuracy_results.csv")