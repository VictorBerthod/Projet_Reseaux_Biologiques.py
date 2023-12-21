import numpy as np

filename = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\sample.txt"


#1.2.1 Question structure 1
def read_interaction_file_dict(filename):
    interaction_dict = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
        first_line = lines[0].strip()
        if first_line.isdigit():
            lines = lines[1:]
        for line in lines:
            proteins = line.strip().split('\t')
            protein1, protein2 = proteins
            if protein1 not in interaction_dict:
                interaction_dict[protein1] = []
            if protein2 not in interaction_dict:
                interaction_dict[protein2] = []
            interaction_dict[protein1].append(protein2)
            interaction_dict[protein2].append(protein1)
    return interaction_dict


#1.2.2 Question structure 2
def read_interaction_file_list(filename):
    interaction_list = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        first_line = lines[0].strip()
        if first_line.isdigit():
            lines = lines[1:]
        for line in lines:
            proteins = line.strip().split('\t')
            interaction_list.append(tuple(proteins))
            interaction_list.append((proteins[1], proteins[0])) #Pour avoir les deux sens (à voir si besoin)
    return interaction_list

# print(read_interaction_file_list(filename))


#1.2.3 Question structure 3
def read_interaction_file_mat(filename):
    adjacency_matrix = None
    nodes = set()
    with open(filename, 'r') as file:
        lines = file.readlines()
        first_line = lines[0].strip()
        if first_line.isdigit():
            lines = lines[1:]
        for line in lines:
            proteins = line.strip().split('\t')
            node1, node2 = proteins
            nodes.add(node1)
            nodes.add(node2)
    nodes_list = list(nodes)
    num_nodes = len(nodes_list)

    # Initialise la matrice d'adjacence avec des zéros
    adjacency_matrix = np.zeros((num_nodes, num_nodes), dtype=int)
    with open(filename, 'r') as file:
        lines = file.readlines()
        first_line = lines[0].strip()
        if first_line.isdigit():
            lines = lines[1:]
        for line in lines:
            proteins = line.strip().split('\t')
            node1, node2 = proteins
            index1 = nodes_list.index(node1)
            index2 = nodes_list.index(node2)

            adjacency_matrix[index1][index2] = 1
            adjacency_matrix[index2][index1] = 1
    return adjacency_matrix, nodes_list

# adj_matrix, nodes = read_interaction_file_mat(filename)
# print("Matrice d'adjacence :")
# print(adj_matrix)
# print("Liste ordonnée des sommets :")
# print(nodes)

#1.2.4 Question structure 4
def read_interaction_file(filename):
    d_int = read_interaction_file_dict(filename)
    l_int = read_interaction_file_list(filename)
    m_int, l_som = read_interaction_file_mat(filename)
    return d_int, l_int, m_int, l_som

# d_int, l_int, m_int, l_som = read_interaction_file(filename)
# print("Dictionnaire représentant le graphe :")
# print(d_int)
# print("Liste d'interactions représentant le graphe :")
# print(l_int)
# print("Matrice d'adjacence correspondant au graphe :")
# print(m_int)
# print("Liste ordonnée des sommets :")
# print(l_som)

#1.2.5 Question structure 5
#Comment optimiser read_interaction_file ? Ne charger uniquement ce qu'on a besoin pour ce qu'on veut faire. 


#2.1.1 Question exploration 1
def count_vertices(filename):
    interaction_list = read_interaction_file_list(filename)
    vertices = set()  # Utilise un ensemble pour stocker les sommets uniques

    for interaction in interaction_list:
        vertices.add(interaction[0])
        vertices.add(interaction[1])

    num_vertices = len(vertices)
    return num_vertices

#2.1.2 Question exploration 2
def count_edges(filename):
    interaction_list = read_interaction_file_list(filename)
    num_edges = len(interaction_list)
    return num_edges



#2.1.3 Question nettoyage
def clean_interactome(filein, fileout):
    interaction_list = read_interaction_file_list(filein)
    unique_interactions = set()
    
    cleaned_interactions = []
    for interaction in interaction_list:
        protein1, protein2 = interaction
        if protein1 != protein2: #Exclu les homo-dimères
            unique_interaction = tuple(sorted([protein1, protein2])) #Ensemble donc plus de doublons
            if unique_interaction not in unique_interactions:
                cleaned_interactions.append(interaction)
                unique_interactions.add(unique_interaction)

    with open(fileout, 'w') as f:
        f.write(f"{len(cleaned_interactions)}\n")
        for interaction in cleaned_interactions:
            f.write(f"{interaction[0]}\t{interaction[1]}\n")

#2.2.1 Question degré 1
def get_degree(file, prot):
    interaction_list = read_interaction_file_list(file)
    degree = 0
    for interaction in interaction_list:
        protein1, protein2 = interaction
        if protein1 == prot or protein2 == prot:
            degree += 1
    return degree

#2.2.2 Question degré 2
def get_max_degree(file):
    interaction_list = read_interaction_file_list(file)
    degree_count_dict = {}
    for interaction in interaction_list:
        protein1, protein2 = interaction
        if protein1 in degree_count_dict:
            degree_count_dict[protein1] += 1
        else:
            degree_count_dict[protein1] = 1
        if protein2 in degree_count_dict:
            degree_count_dict[protein2] += 1
        else:
            degree_count_dict[protein2] = 1
    max_protein = max(degree_count_dict, key=degree_count_dict.get)
    max_degree = degree_count_dict[max_protein]

    return max_protein, max_degree

#2.2.3 Question degré 3
def get_ave_degree(file):
    interaction_list = read_interaction_file_list(file)
    degree_count_dict = {}
    for interaction in interaction_list:
        protein1, protein2 = interaction
        if protein1 in degree_count_dict:
            degree_count_dict[protein1] += 1
        else:
            degree_count_dict[protein1] = 1
        if protein2 in degree_count_dict:
            degree_count_dict[protein2] += 1
        else:
            degree_count_dict[protein2] = 1
    
    total_degree = sum(degree_count_dict.values())
    num_proteins = len(degree_count_dict)

    if num_proteins > 0:
        average_degree = total_degree / num_proteins
    else:
        average_degree = 0

    return average_degree


#2.2.4 Question degré 4
def count_degree(file, deg):
    interaction_list = read_interaction_file_list(file)
    degree_count_dict = {}
    for interaction in interaction_list:
        protein1, protein2 = interaction
        if protein1 in degree_count_dict:
            degree_count_dict[protein1] += 1
        else:
            degree_count_dict[protein1] = 1
        if protein2 in degree_count_dict:
            degree_count_dict[protein2] += 1
        else:
            degree_count_dict[protein2] = 1
    # print(degree_count_dict)
    count=0
    for degree in degree_count_dict.values():
        if degree == deg:
            count += 1
    return count


#2.2.5 Question degré 5
def histogram_degree(file, dmin, dmax):
    interaction_list = read_interaction_file_list(file)
    degree_count_dict = {}
    for interaction in interaction_list:
        protein1, protein2 = interaction
        if protein1 in degree_count_dict:
            degree_count_dict[protein1] += 1
        else:
            degree_count_dict[protein1] = 1
        if protein2 in degree_count_dict:
            degree_count_dict[protein2] += 1
        else:
            degree_count_dict[protein2] = 1

    histogram = {}
    for degree in degree_count_dict.values():
        if dmin <= degree <= dmax:
            if degree in histogram:
                histogram[degree] += 1
            else:
                histogram[degree] = 1
    for degree in range(dmin, dmax + 1):
        count = histogram.get(degree, 0)
        stars = '*' * count
        print(f"{degree} {stars}")






