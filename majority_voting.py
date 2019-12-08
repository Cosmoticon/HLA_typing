import random
import pandas as pd
import json
import sys
from collections import Counter

def main():
    """
    Extracts the HLA typing alleles from the CSV files, 
    runs the custom majority voting algorithm, and saves the result in a CSV output file
    """
    file_path = sys.argv[1]
    results = pd.read_csv(file_path)
    dict_majority = {"algorithm": "majority_voting"}
    for column in ["A", "B", "C", "DQB1", "DRB1"]:
        total_alleles1 = results[column + "-1"].tolist()
        total_alleles2 = results[column + "-2"].tolist()
        paired_alleles = tuple(zip(total_alleles1, total_alleles2))
        allele1 = find_best(paired_alleles)
        for tup in paired_alleles:
            if allele1 in tup:
                tup.remove(allele1)
        allele2 = find_best(paired_alleles)
        dict_majority[column + "-1"] = allele1
        dict_majority[column + "-2"] = allele2
    results.append(dict_majority)
    results.to_csv(file_path.replace(".csv", "maj_voting.csv"))
        
def find_best(paired_values):
    """
    Finds the best allele to output, given the list of allele pairs
    (If there is a tie, it is broken randomly)
    """
    total_values = list(sum(paired_values, ()))
    count_values = Counter(total_values).most_common()
    max_count = count_values[0][1]
    list_max = []
    for tuple_count in count_values:
        if tuple_count[1] == max_count:
            list_max.append(tuple_count[0])
    while len(list_max) > 0:
        random_max = random.choice(list_max)
        if total_values.count(random_max) >= (len(total_values) / 2):
            return random_max
        else:
            count_tuples = 0
            for tup in paired_values:
                if random_max in tup:
                    count_tuples
            if count_tuples >= (len(paired_values) / 2):
                return random_max
    return ""

if __name__ == "__main__":
    main()