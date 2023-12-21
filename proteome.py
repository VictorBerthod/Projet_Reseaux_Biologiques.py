class Proteome():
    def __init__(self, filename_bdd):
        """
        Initialize a Proteome object.

        Parameters:
        -----------
        filename_bdd : str
            Name of the file containing protein data.

        Attributes:
        -----------
        proteins : set
            Set containing unique proteins.
        prot_to_uniprot_id : dict
            Dictionary mapping proteins to their UniProt IDs.
        uniprot_id_to_prot : dict
            Dictionary mapping UniProt IDs to proteins.
        """
        self.proteins = set()
        self.prot_to_uniprot_id = {} 
        self.uniprot_id_to_prot = {}
        if filename_bdd != "":
            self.set_dict(filename_bdd)

    def get_proteins(self):
        """
        Get the set of unique proteins.

        Returns:
        --------
        set
            Set containing unique proteins.
        """
        return self.proteins  
    
    def get_prot_to_uniprot_id(self):
        """
        Get the dictionary mapping proteins to their UniProt IDs.

        Returns:
        --------
        dict
            Dictionary mapping proteins to their UniProt IDs.
        """
        return self.prot_to_uniprot_id
    
    def get_uniprot_id_to_prot(self):
        """
        Get the dictionary mapping UniProt IDs to proteins.

        Returns:
        --------
        dict
            Dictionary mapping UniProt IDs to proteins.
        """
        return self.uniprot_id_to_prot
    
    def set_dict(self, filename_bdd):
        """
        Set the dictionaries for proteins and their corresponding UniProt IDs.

        Parameters:
        -----------
        filename_bdd : str
            Name of the file containing protein data.
        """
        with open(filename_bdd, 'r') as fichier:
            for ligne in fichier:
                elements = ligne.split('|')
                uniprot_id = elements[1]
                protein_element = elements[2].split(' ')
                protein = protein_element[0]
                self.proteins.add(protein)
                self.prot_to_uniprot_id[protein] = uniprot_id
                self.uniprot_id_to_prot[uniprot_id] = protein
                # print("Protein : " + protein + " - Uniprot ID : " + uniprot_id)




# filename = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\sample.txt"
# filename_prot_id = "C:\\Users\\victo\\Documents\\Victor\\Cours\\M2 Bio-informatique - InterBio (Rennes)\\Cours\\Projet\\uniprot_out.txt" #570 155 lignes
# proteome = Proteome(filename_prot_id)
# print(proteome.proteins)
# print(proteome.prot_to_uniprot_id)
# print(proteome.uniprot_id_to_prot)


#INSR_HUMAN     P06213
#https://rest.uniprot.org/uniprotkb/P06213.txt
# Domain		624-726	Fibronectin type-III 1 	
# Domain		757-842	Fibronectin type-III 2 	
# Domain		853-947	Fibronectin type-III 3 	
# Domain		1023-1298	Protein kinase
