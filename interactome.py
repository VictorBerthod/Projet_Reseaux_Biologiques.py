import random
import string
import tempfile
import urllib
import numpy as np
from urllib import request
import re

from proteome import Proteome

class Interactome(Proteome):
    """
    Interactome class representing a collection of protein interactions.
    """
    def __init__(self, filename, proteome):
        """
        Initialize an Interactome object.

        Parameters:
        -----------
        filename : str or None
            Name of the file containing interaction data. If None, initializes an empty Interactome.
        proteome : Proteome
            Instance of the Proteome class.

        Attributes:
        -----------
        int_list : list
            List containing interactions between proteins.
        int_dict : dict
            Dictionary storing protein interactions as keys and their connected neighbors as values.
        proteins : set
            Set containing unique proteins present in the interactome.
        proteome : Proteome
            Instance of the Proteome class.
        dict_graph : dict
            Dictionary containing protein details like UniProt IDs, neighbors, and domains.
        """
        self.int_list = []
        self.int_dict = {}
        self.proteins = set()
        self.proteome = proteome

        self.dict_graph = {}
        if filename : 
            self.read_interaction_file(filename)
            self.xlink_Uniprot() #Initialise le dictionnaire dict_graph sans domaines
            # self.xlink_domains() #Ajoute les domaines au dictionnaire dict_graph 


    def get_int_dict(self):
        """
        Get the dictionary storing protein interactions.
        Returns:
        --------
        dict
            Dictionary with proteins as keys and their connected neighbors as values.
        """
        return self.int_dict
    def get_int_list(self):
        """
        Get the list containing interactions between proteins.
        Returns:
        --------
        list
            List of protein interactions.
        """
        return self.int_list
    def get_proteins(self):
        """
        Get the set of unique proteins present in the interactome.
        Returns:
        --------
        set
            Set containing unique proteins.
        """
        return self.proteins

    def xlink_Uniprot(self):
        """
        Create a dictionary containing protein details like UniProt IDs, neighbors, and domains.

        This method populates the 'dict_graph' attribute with protein information. Each protein in the interactome
        is represented as a key in the dictionary with the following details:
        - 'UniprotID': UniProt ID of the protein if available, otherwise 'None'.
        - 'voisins': List of neighbors connected to the protein based on stored interactions.
        - 'domains': Empty list representing the protein's domains (to be populated later).

        The 'dict_graph' attribute is used to store detailed information about each protein's interactions and UniProt IDs for further analysis or processing.
        """
        dict_graph = {}
        for protein in self.proteins:
            dict_graph[protein] = {
                'UniprotID': self.proteome.prot_to_uniprot_id[protein] if protein in self.proteome.prot_to_uniprot_id else 'None',                
                'voisins': self.int_dict[protein],
                'domains': []
                }
        self.dict_graph = dict_graph

    def get_protein_domains(self, p):
        """
        Retrieve the domains associated with a specific protein using its UniProt ID.

        Parameters:
        -----------
        p : str
            Name of the protein for which domains are requested.

        Returns:
        --------
        list
            List of domains associated with the protein. If the protein or its UniProt ID is not found,
            an empty list is returned.
        """
        # list_no_ID=[]
        if p not in proteome.proteins:
            return([])        
        if proteome.prot_to_uniprot_id[p] == None:
            return([])
        else:
            UniprotID = proteome.prot_to_uniprot_id[p]
            url = "https://rest.uniprot.org/uniprotkb/" + UniprotID + ".txt"
            # print(url)
            try:
                response = urllib.request.urlopen(url)
            except urllib.error.HTTPError as e:
                if e.code == 400:
                    print("URL Problem : ",UniprotID)
                    return([])
            with urllib.request.urlopen(url) as response:
                data = response.read()
                # print(data)
                with tempfile.NamedTemporaryFile(mode='w', delete=True) as f:
                    f.write(data.decode('utf-8'))
                    #print(f.name)
                    f.close()
            expression = r'FT\s*DOMAIN\s+[\s\S]*?/note="(.*?)"'
            res = re.findall(expression, data.decode('utf-8'))
            return(res)

    def xlink_domains(self):
        """
        Retrieve and associate protein domains from UniProt for each protein in the interactome.
        This method fetches the protein domains information from UniProt and associates them with
        the respective proteins in the interactome's dictionary representation.

        Note:
        -----
        This method requires access to the internet to fetch data from UniProt.

        Returns:
        --------
        None
        """
        total = len(self.proteins)
        i = 1
        for p in self.proteins:
            # print("Progression : " + str(i) + "/" + str(total))
            i += 1
            self.dict_graph[p]['domains'] = self.get_protein_domains(p)

    def open_dict_grah_saved(self, file):
        """
        Open and read a saved dictionary graph from a file.

        Parameters:
        -----------
        file : str
            Path to the file containing the saved dictionary graph.

        Returns:
        --------
        dict
            A dictionary with protein details such as UniProt IDs, neighbors, and domains.
        """
        with open(save_file, 'r') as file:
            dict_readed = {}
            lignes = file.readlines()
            for ligne in lignes:
                protein, donnees = ligne.split(" {'")

                debut_uniprot_id = ligne.find("UniprotID': '") + len("UniprotID': '")
                find_uniprot_id = ligne.find("', 'voisins': [")
                uniprot_id = ligne[debut_uniprot_id:find_uniprot_id]

                debut_voisins = ligne.find("'voisins': [") + len("'voisins': [")
                fin_voisins = ligne.find("]", debut_voisins)
                voisins_str = ligne[debut_voisins:fin_voisins]
                voisin = (eval("["+voisins_str+"]"))

                debut_domains = ligne.find("'domains': [") + len("'domains': [")
                fin_domains = ligne.find("]}", debut_domains)
                domains_str = ligne[debut_domains:fin_domains]
                domains = (eval("["+domains_str+"]"))

                dict_readed[protein] = {
                        'UniprotID': uniprot_id,                
                        'voisins': voisin,
                        'domains': domains
                        }
            return dict_readed

    def read_interaction_file(self, filename):
        """
        Read interaction data from a file and populate the interactome attributes.

        Parameters:
        -----------
        filename : str
            Name of the file containing interaction data.

        Notes:
        ------
        - The file should contain interactions between proteins, each line indicating an interaction between two proteins.
        - Proteins are identified by their names.
        """
        with open(filename, 'r') as file:
            lines = file.readlines()
            first_line = lines[0].strip()
            if first_line.isdigit():
                lines = lines[1:]
            for line in lines:
                protein1, protein2 = line.strip().split('\t')
                self.int_list.append((protein1, protein2))
                if protein1 in self.int_dict:
                    self.int_dict[protein1].append(protein2)
                else:
                    self.int_dict[protein1] = [protein2]
                if protein2 in self.int_dict:
                    self.int_dict[protein2].append(protein1)
                else:
                    self.int_dict[protein2] = [protein1]
                self.proteins.add(protein1)
                self.proteins.add(protein2)

    def print_interaction_list(self):
        """
        Print the list of interactions between proteins.

        This method iterates through the list of interactions and prints each interaction.
        Each interaction consists of a pair of proteins.
        """
        for interaction in self.int_list:
            print(interaction)

    def print_interaction_dict(self):
        """
        Print the dictionary representing protein interactions.

        This method iterates through the protein interaction dictionary and prints each protein along
        with its connected neighbors.
        """
        for protein, neighbors in self.int_dict.items():
            print(f"{protein}: {neighbors}")

    def get_adjacency_matrix(self):
        """
        Generate and return the adjacency matrix representing protein interactions.

        Returns:
        --------
        list
            List of lists representing the adjacency matrix.
            Each inner list corresponds to a row in the matrix, where the presence of an interaction is denoted by 1,
            and absence by 0. The order of proteins corresponds to the row and column indices.
        """
        adjacency_matrix = []
        for protein1 in self.proteins:
            row = []
            for protein2 in self.proteins:
                if protein1 in self.int_dict and protein2 in self.int_dict[protein1]:
                    row.append(1)
                else:
                    row.append(0)
            adjacency_matrix.append(row)
        return adjacency_matrix

    
    def read_interaction_file_mat(self):
        """
        Create an adjacency matrix representation of the interactome.

        Returns:
        --------
        tuple
            A tuple containing the adjacency matrix and a list of ordered proteins.
            - adjacency_matrix (numpy.ndarray): An adjacency matrix representing protein interactions.
            - ordered_proteins (list): A list of proteins ordered for the adjacency matrix.
        """
        ordered_proteins = sorted(list(self.proteins))
        adjacency_matrix = np.zeros((len(ordered_proteins), len(ordered_proteins)), dtype=int)
        for interaction in self.int_list:
            protein1, protein2 = interaction
            index1 = ordered_proteins.index(protein1)
            index2 = ordered_proteins.index(protein2)
            adjacency_matrix[index1, index2] = 1
            adjacency_matrix[index2, index1] = 1
        return adjacency_matrix, ordered_proteins

    def count_vertices(self):
        """
        Count the number of vertices (proteins) in the interactome.

        Returns:
        --------
        int
            Number of proteins in the interactome.
        """
        return len(self.proteins)

    def count_edges(self):
        """
        Count the number of edges (interactions) in the interactome.

        Returns:
        --------
        int
            Number of interactions in the interactome.
        """
        return len(self.int_list)

    def clean_interactome(self):
        """
        Clean the interactome to remove duplicate interactions.

        Returns:
        --------
        Interactome
            An Interactome object with cleaned interactions.
        """
        cleaned_interactome = Interactome(None)
        for interaction in self.int_list:
            protein1, protein2 = interaction
            if (protein1, protein2) not in cleaned_interactome.int_list and (protein2, protein1) not in cleaned_interactome.int_list:
                if protein1 != protein2:
                    cleaned_interactome.int_list.append((protein1, protein2))
                    if protein1 in cleaned_interactome.int_dict:
                        cleaned_interactome.int_dict[protein1].append(protein2)
                    else:
                        cleaned_interactome.int_dict[protein1] = [protein2]
                    if protein2 in cleaned_interactome.int_dict:
                        cleaned_interactome.int_dict[protein2].append(protein1)
                    else:
                        cleaned_interactome.int_dict[protein2] = [protein1]
                    cleaned_interactome.proteins.add(protein1)
                    cleaned_interactome.proteins.add(protein2)
        return cleaned_interactome

    def get_degree(self, prot):
        """
        Get the degree of a specified protein in the interactome.

        Parameters:
        -----------
        prot : str
            Protein name.

        Returns:
        --------
        int
            Degree of the specified protein.
        """
        if prot in self.int_dict:
            return len(self.int_dict[prot])
        else:
            return 0

    def get_max_degree(self):
        """
        Get the protein with the maximum degree in the interactome.

        Returns:
        --------
        tuple
            A tuple containing the protein with the maximum degree and its degree.
            - max_protein (str): Protein with the maximum degree.
            - max_degree (int): Degree of the protein.
        """
        max_degree = 0
        max_protein = None
        for protein in self.proteins:
            degree = self.get_degree(protein)
            if degree > max_degree:
                max_degree = degree
                max_protein = protein

        return max_protein, max_degree

    def get_ave_degree(self):
        """
        Calculate the average degree of all proteins in the interactome.

        Returns:
        --------
        float
            Average degree of proteins in the interactome.
        """
        total_degree = sum(self.get_degree(prot) for prot in self.proteins)
        return total_degree / len(self.proteins)

    def count_degree(self, deg):
        """
        Count the number of proteins with a specified degree.

        Parameters:
        -----------
        deg : int
            Degree to count.

        Returns:
        --------
        int
            Number of proteins with the specified degree.
        """
        count = 0
        for protein in self.proteins:
            if self.get_degree(protein) == deg:
                count += 1
        return count

    def histogram_degree(self, dmin, dmax):
        """
        Generate a histogram of degrees within a specified range.

        Parameters:
        -----------
        dmin : int
            Minimum degree for the histogram.
        dmax : int
            Maximum degree for the histogram.
        """
        histogram = {}
        for protein in self.proteins:
            degree = self.get_degree(protein)
            if dmin <= degree <= dmax:
                if degree in histogram:
                    histogram[degree] += 1
                else:
                    histogram[degree] = 1

        for degree in range(dmin, dmax + 1):
            count = histogram.get(degree, 0)
            stars = '*' * count
            print(f"{degree} {stars}")
    
    def density(self):
        """
        Calculate the density of the interactome.

        Returns:
        --------
        float
            Density of the interactome.
        """
        num_vertices = len(self.proteins)
        num_edges = len(self.int_list)
        if num_vertices < 2:
            return 0
        max_possible_edges = (num_vertices * (num_vertices - 1))
        density = (2 * num_edges) / max_possible_edges
        return density

    
    def clustering(self, prot):
        """
        Calculate the clustering of a specified protein.

        Parameters:
        -----------
        prot : str
            Protein name.

        Returns:
        --------
        list
            List of neighboring proteins for the specified protein.
        """
        if prot not in self.get_proteins():
            print(f"Error : no protein named {prot} doesn't exist in interactome.")
            exit()
        if self.get_degree(prot) == 0 : #sécu
            return 0
        neighbor = []
        for protein in self.int_dict[prot]:
            neighbor.append(protein)
        return(neighbor)
        
    def reccursive_search(self, prot):
        """
        Perform a recursive search for neighboring proteins.

        Parameters:
        -----------
        prot : str
            Protein name.

        Returns:
        --------
        list
            List of neighboring proteins obtained through recursive search.
        """
        list_neighbors = []
        list_neighbors.append(self.clustering(prot))
        if len(list_neighbors) == 0:
            return list_neighbors
        else:
            for neighbor in list_neighbors:
                list_neighbors.append(self.reccursive_search(neighbor))
        return list_neighbors

    def grapher(self, n, p):
        """
        Generate a random graph using the Erdos Renyi model.

        Parameters:
        -----------
        n : int
            Number of proteins in the graph.
        p : float
            Probability of edge creation between proteins.
        """
        new_int_mat = []
        new_int_dict = {}
        nbr_edges = n * (n - 1) // 2
        connection = np.random.binomial(nbr_edges, p)
        alphabet = string.ascii_lowercase
        self.set_proteins(random.sample(alphabet, n))

        relations = set()
        while len(relations) < connection:
            proteins = random.sample(self.get_proteins(), 2)
            if proteins[0] != proteins[1]:
                relation = frozenset(proteins)
                relations.add(relation)
        self.set_int_list(list(relations))

        for prot_names in self.get_proteins():
            line_mat = []
            for prot_int in self.get_proteins():
                if {prot_names, prot_int} in relations:
                    line_mat.append(1)
                else:
                    line_mat.append(0)
            new_int_mat.append(line_mat)
        self.set_int_mat(new_int_mat)

        for interaction in self.get_int_list():
            for prot in interaction:
                if prot in new_int_dict:
                    new_int_dict[prot]["voisins"] = list(set(new_int_dict[prot]["voisins"]) | {p for p in interaction if p != prot})
                else:
                    new_int_dict[prot] = {"voisins": list({p for p in interaction if p != prot})}
        self.set_int_dict(new_int_dict)    
    

    def graphba(self, n):
        """
        Generate a random graph using the Barabasi Albert model.

        Parameters:
        -----------
        n : int
            Number of proteins in the graph.
        """
        alphabet = string.ascii_lowercase
        list_prot = random.sample(alphabet, n)
        self.set_proteins(list_prot)

        tree = [[list_prot[0], list_prot[1]]]
        self.set_int_list(tree)

        int_dict = {list_prot[0]: {"voisins": [list_prot[1]]}, list_prot[1]: {"voisins": [list_prot[0]]}}

        for i in range(2, n):
            prot = list_prot[i]
            dictio = self.get_int_dict()
            for neigh in list_prot[:i]:
                tot_connections = len(tree) * 2
                proba = len(int_dict[neigh]["voisins"]) / tot_connections if neigh in int_dict else 0
                if random.random() <= proba:
                    tree.append([neigh, prot])

                    if prot not in dictio:
                        dictio[prot] = {"voisins": [neigh]}
                    else:
                        dictio[prot]["voisins"].append(neigh)

                    if neigh not in dictio:
                        dictio[neigh] = {"voisins": [prot]}
                    else:
                        dictio[neigh]["voisins"].append(prot)

            int_dict = dictio
            self.set_int_dict(int_dict)
            self.set_int_list(tree)


#--------------------------------------------------------------------------------------------------------------------------------------------------------

# filename = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\sample2.txt"
filename2 = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\human_high_quality.txt"
# file_bdd_uniprot = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\sample_uniprot_out.txt"
file_bdd_uniprot2 = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\uniprot_out.txt"

proteome = Proteome(file_bdd_uniprot2) #Création du proteome : dictionnaire protéine <-> uniprot_id
interactome = Interactome(filename2, proteome) #Création de l'interactome 

# Sauvegarde du dictionnaire dict_graph dans un fichier texte
# dict_graph_all=(interactome.dict_graph)
# for key in dict_graph_all:
#     print(key, dict_graph_all[key])
# file_to_save = "dict_graph_human_high_quality_SAVE.txt"
# with open(file_to_save, 'w') as file:
#     for key in dict_graph_all:
#         file.write('{} {}\n'.format(key, dict_graph_all[key]))

# Lecture du fichier dict_graph_human_high_quality_SAVE.txt pour récupérer le dictionnaire dict_graph
save_file = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\dict_graph_human_high_quality_SAVE.txt"
dict_readed = interactome.open_dict_grah_saved(save_file)

#--------------------------------------------------------------------------------------------------------------------------------------------------------
#Zone de test

proteins = dict_readed.keys()

#Nombre de protéines sans domaines ?
# cpt = 0
# for p in proteins:
#     if dict_readed[p]['domains'] == []:
#         cpt += 1
# print(cpt)

# Nombre de domaines par proteine ?
# nb_domaines_prot=dict()
# for p in proteins:
#     nb_domaines_prot[p] = len(dict_readed[p]['domains'])
# valeurs_triees = sorted(nb_domaines_prot, key=nb_domaines_prot.get, reverse=True)
# for v in valeurs_triees:
#     if nb_domaines_prot[v] != 0:
#         print(v, nb_domaines_prot[v])

# Pour stocker les domaines dans un fichier texte
# file_to_save = "proteine_domaine.txt"
# with open(file_to_save, 'w') as file:
#     for v in valeurs_triees:
#         text=str(v)+"    "+str(dico[v])
#         file.write(text)
#         file.write('\n')

#Nombre occurences domaines ?
# cpt=dict()
# for p in proteins:
#     for domain in dict_readed[p]['domains'] :
#         if domain in cpt:
#             cpt[domain] += 1
#         else:
#             cpt[domain] = 1
# valeurs_triees = sorted(cpt, key=cpt.get, reverse=True)[:10]
# for v in valeurs_triees:
#     print(v, cpt[v])
# print(len(cpt.keys())) #1823 domaines différents

# Voisins co-occurence pour graph domaine ? TOP10
# dict_graph_domain=dict()
# for p in proteins:
#     for domain in dict_readed[p]['domains'] :
#         if domain is not None:
#             if domain not in dict_graph_domain:
#                 dict_graph_domain[domain] = dict_readed[p]['domains']
#             else:
#                 clean_list = dict_readed[p]['domains'].remove(domain)
#                 if clean_list is not None:
#                     dict_graph_domain[domain].append(clean_list)

# print(dict_graph_domain)
# dic_nb_voisins=dict()
# for d in dict_graph_domain.keys():
#     dic_nb_voisins[d] = len(dict_graph_domain[d])
# # print(dic_nb_voisins)

# voisins_trie = sorted(dic_nb_voisins, key=dic_nb_voisins.get, reverse=True)[10:]
# ii=0
# for v in voisins_trie:
#     if dic_nb_voisins[v] == 1:
#         ii+=1
#         print(v, dic_nb_voisins[v])
# print(ii)

# valeurs_triees2 = sorted(dict_graph_domain, key=dict_graph_domain.get, reverse=True)[:100]
# for v2 in valeurs_triees2:
#     print(v2, dict_graph_domain[v2])


# print(interactome.dict_graph['INSR_HUMAN']['UniprotID'])
# print(interactome.dict_graph['Q8NEF3_HUMAN']['voisins'])

# print(interactome.get_protein_domains("INSR_HUMAN"))

#print(proteome.prot_to_uniprot_id['INSR_HUMAN'])

# interactome.xlink_Uniprot()

# list = list(interactome.proteins)
# print(interactome.int_dict)
# print(interactome.clustering('rad53'))
#print(interactome.reccursive_search('rad53'))

# interactome.print_interaction_list()
# interactome.print_interaction_dict()
# adjacency_matrix = interactome.get_adjacency_matrix()
# print(adjacency_matrix)
# interactome.read_interaction_file_mat()
# cleaned_interactome = interactome.clean_interactome()

# prot_degree = interactome.get_degree('prot1')
# max_prot, max_degree = interactome.get_max_degree()
# ave_degree = interactome.get_ave_degree()
# count_deg_3 = interactome.count_degree(3)
# interactome.histogram_degree(1, 5)
# print(interactome.density())
# interactome.clustering('rad53')

# interactome_er = Interactome(filename)
# interactome_er.graph_er(0.3)
# interactome_er.print_interaction_dict()
# interactome_er.histogram_degree(1, 10)

# interactome_ba = Interactome(None)
# interactome_ba.graph_ba(3)
# interactome_ba.print_interaction_dict()
# interactome_ba.histogram_degree(1, 10)


