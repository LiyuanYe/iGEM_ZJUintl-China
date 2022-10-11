# mouse
    # A in vitro hep & KP & 152329 (3) √
    # B in vivo mouse organ atlas: marrow, liver, lung (3) df √
    # C in vivo aging 100906 (1) √

# human
    # D protein serum chip-seq limma (1) √
    # E protein BM (1) √
    # F rna-seq BM (1) √

# rat
    # G BM (1) √

# read in
from uniprot_web_scraping import *
def read_markers(filepath):
    with open(filepath, 'r') as infile:
        lines = infile.readlines()
        potential_markers = [line.split("\n")[0] for line in lines]
        infile.close()
    return potential_markers

    # mouse
A = read_markers("D:\\SRTP2021\\potential_markers\\python_selection\\A.csv")
B = read_markers("D:\\SRTP2021\\potential_markers\\python_selection\\B.csv")
C = read_markers("D:\\SRTP2021\\potential_markers\\python_selection\\C.csv")
    # human
# D = read_markers("D:\\SRTP2021\\potential_markers\\python_selection\\D.csv")
E = read_markers("D:\\SRTP2021\\potential_markers\\python_selection\\EX.csv")
F = read_markers("D:\\SRTP2021\\potential_markers\\python_selection\\F.csv")
    # rat
G = read_markers("D:\\SRTP2021\\potential_markers\\python_selection\\G.csv")

# preprocess
def no_dup_na(data):
    no_dup_na = []
    for gene in data:
        if gene not in no_dup_na and gene != "NA" and gene != 'None':
            no_dup_na.append(gene.lower())
    return no_dup_na

A = no_dup_na(A)
B = no_dup_na(B)
C = no_dup_na(C)
# D = no_dup_na(D)
E = no_dup_na(E)
F = no_dup_na(F)
G = no_dup_na(G)

# count occurrences
res_dict = {}
data_merged = A + B + C + E + F + G
for gene in data_merged:
    if gene.lower() not in res_dict:
        res_dict[gene.lower()] = 1
    else:
        res_dict[gene.lower()] += 1
# max = 4
# max(res_dict.values())

def get_gene_by_n_occurrence(res_dict, n):
    return [gene for gene in res_dict if res_dict[gene] == n]
four_times_genes = get_gene_by_n_occurrence(res_dict, 4)
three_times_genes = get_gene_by_n_occurrence(res_dict, 3)
twice_genes = get_gene_by_n_occurrence(res_dict, 2)
once_genes = get_gene_by_n_occurrence(res_dict, 1)

# find membrane proteins
four_mem = main_function(object=four_times_genes, mode='object') # 0/3
three_mem = main_function(object=three_times_genes, mode='object') # /37
twice_mem = main_function(object=twice_genes, mode='object') # /274

# which three datasets?
def which_datasets(markers):
    dict_datasets = {}
    for gene in markers:
        dict_datasets[gene] = []
        if gene in A:
            dict_datasets[gene].append('A')
        if gene in B:
            dict_datasets[gene].append('B')
        if gene in C:
            dict_datasets[gene].append('C')
        if gene in E:
            dict_datasets[gene].append('E')
        if gene in F:
            dict_datasets[gene].append('F')
        if gene in G:
            dict_datasets[gene].append('G')
    return dict_datasets

gene_datasets = which_datasets(twice_genes)
for gene in gene_datasets:
    print(gene, ": ", gene_datasets[gene])

# which datasets of 10 paper markers
paper_markers = read_markers('D:\\SRTP2021\\potential_markers\\paper_10_markers.txt')
paper_markers_lower = [gene.lower() for gene in paper_markers]
gene_datasets = which_datasets(paper_markers_lower)
for gene in gene_datasets:
    print(gene, ": ", gene_datasets[gene])

# intersect best 500 mem scores
marrow = read_markers("D:\\SRTP2021\\12\\scores\\df_marrow_best_500_memscore.csv")
for i in range(len(marrow)):
    marrow[i] = marrow[i].split(',')[0]

lung = read_markers("D:\\SRTP2021\\12\\scores\\df_lung_best_500_memscore.csv")
for i in range(len(lung)):
    lung[i] = lung[i].split(',')[0]

liver = read_markers("D:\\SRTP2021\\12\\scores\\df_liver_best_500_memscore.csv")
for i in range(len(liver)):
    liver[i] = liver[i].split(',')[0]

three_organs = marrow + lung + liver

cnt_dict = {}
for gene in three_organs:
    if gene not in cnt_dict:
        cnt_dict[gene] = 1
    else:
        cnt_dict[gene] += 1

intersect = [gene for gene in cnt_dict if cnt_dict[gene] == 3]

with open("D:\\SRTP2021\\12\\scores\\three_organs_mem.csv", 'w') as out:
    for gene in intersect:
        out.write(gene+"\n")


# Heart Selection
b_heart = read_markers('D:/SRTP2021/potential_markers/python_selection/heart/B_heart_limma.csv')
g_aorta = read_markers('D:/SRTP2021/potential_markers/python_selection/heart/G_heart_converted.csv')
monkey_AA = read_markers('D:/SRTP2021/potential_markers/python_selection/heart/monkey_AA_genes.csv')
monkey_CA = read_markers('D:/SRTP2021/potential_markers/python_selection/heart/monkey_CA_genes.csv')

b_heart = no_dup_na(b_heart)
g_aorta = no_dup_na(g_aorta)
monkey_AA = no_dup_na(monkey_AA)
monkey_CA = no_dup_na(monkey_CA)

res_dict = {}
data_merged = b_heart + g_aorta + monkey_AA + monkey_CA
for gene in data_merged:
    if gene.lower() not in res_dict:
        res_dict[gene.lower()] = 1
    else:
        res_dict[gene.lower()] += 1

four_times_genes = get_gene_by_n_occurrence(res_dict, 4)
three_times_genes = get_gene_by_n_occurrence(res_dict, 3)
twice_genes = get_gene_by_n_occurrence(res_dict, 2)
once_genes = get_gene_by_n_occurrence(res_dict, 1)

with open('D:/SRTP2021/potential_markers/python_selection/heart/all_four_three_markers.csv', 'w') as out:
    for gene in four_times_genes:
        out.write(gene+'\n')
    for gene in three_times_genes:
        out.write(gene+'\n')
    out.close()
