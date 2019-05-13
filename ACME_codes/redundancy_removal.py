def redundancy_removal(data_dict):
    '''
    Removes the redundant data from the training set
    '''
    for allele in sorted(data_dict.keys()):
        allele_data = data_dict[allele]
        unique_9mers = []
        nonredundant_data = []
        overlap = 0
        for pep_info in allele_data:
            redundant = False
            seq = pep_info[-1]
            for i in range(len(seq) - 9 + 1):
                _9mer = seq[i:i + 9]
                if _9mer not in unique_9mers:
                    unique_9mers.append(_9mer)
                else:
                    redundant = True
            if not redundant or len(seq) == 8:
                nonredundant_data.append(pep_info)
        data_dict[allele] = nonredundant_data
    
    return data_dict