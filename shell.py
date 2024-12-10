#!python
import numpy as np
import pandas as pd
import re

filename="voro.txt"

def line_parser(line):
    #extract name and index in the line of the cell:
    #   "Face  1:  6 Vertices, Area= 9.6354 A^2, Neighbor is 13085 (C6H11N2[689] H9)"
    #   name = C6H11N2, idx = 689
    pattern = r'\(([\w]+)\[(\d+)\]' # finds the (name[idx]) pattern
    match = re.search(pattern,line)
    if match:
        name = match.group(1)
        idx = match.group(2)
    return name,idx 

def cleaning(cell_mols):
    #remove duplicates from cell array from voro.txt
    #use both name and index to check molecules
    #output is a [n_cells,n_mols] array with mol names of shell composition
    uniques = set(cell_mols)
    mols = []
    for mol in uniques:
        if not mol in mols: #remove douplicates, uses both idx and name
            mols.append(mol)
    return mols

def one_shot_encoding(mols):
    #take array [n_cells,n_mols]
    #output dataframe with molecules name in columns
    #number of neighbours in rows. each row is a cell
    data = []
    for cell in mols:
        cell_dict = {'Al2Cl7':0,'AlCl2':0,'AlCl3':0,'AlCl4':0,'C6F2H4':0,'CH4N2O':0,'C6H11N2':0}
        keys, counts = np.unique(cell,return_counts=True)
        for idx,key in enumerate(keys):
            cell_dict[key]=counts[idx]
        data.append(cell_dict)
    return data

class Analyzer:
    def __init__(self, filename):
        self.mols , self.ref_mols = self.load_file(filename)
        self.shells, self.ref_uniq = self.shell()
        self.df = self.dataframe()
        print(self.df)

    def __str__(self):
        return (r"Data Object: \n mols : [N_cells,N_mols]\n ref_mols : [N_cells,1]")
    def load_file(self,filename):
        # output mols:[n_cells,n_mol] each row is a shell.
        # The same molecules is composed of many cells
        # is required to reunite the cells shell compostion
        # of several cells that comes from the same molecule
        # name_idx to identify molecules
        mols = []
        ref_mols = []
        cell_mols = []
        with (open(filename) as file):
            while True:
                line = file.readline()
                line = line.strip()
                if line.startswith("Cell"):  # start of each cell
                    name, idx = line_parser(line)
                    ref_mols.append(name + "_" + idx)
                    if cell_mols:
                        tmp_mols = cleaning(cell_mols)
                        mols.append(tmp_mols)
                    cell_mols = []
                elif line == '':  # end of file
                    tmp_mols = cleaning(cell_mols)
                    mols.append(tmp_mols)
                    break
                elif line.startswith("- Face"):  # lines inside a cell
                    molname, molidx = line_parser(line)
                    cell_mols.append(molname + "_" + molidx)
        return mols, ref_mols

    def shell(self):
        # To unite all the shell compositions of cells that
        # come from the same molecules. Unique RM (reference molecules)
        # are used to find the right rows of mols array and unreavel in
        # a simple array to be further processed. RM and duplicates need to be removed
        ref_unique = list(set(self.ref_mols))
        #loop over all uniques references
        shells = []
        for ref in ref_unique:
             tmp_mols = self.mols
             ref_mols = np.array(self.ref_mols)
             #search for the same reference molecule and unite shells
             idx = ref_mols == ref
             idx=np.where(idx)[0]
             shell = [ tmp_mols[i]  for i in idx ] # a list of lists
             #if len(mols)>1:
             shell = [ mol for cell in shell for mol in cell ]  # flatten
             shell = cleaning(shell)                 # remove double
             #removing reference molecule
             #shell = np.array(shell)
             if ref in shell:
                 idx = np.array(shell) == ref
                 idx = np.invert(idx)
                 shell = [ shell[i].split("_")[0] for i in np.where(idx)[0] ] # removed RM and idx
             shells.append(shell)
        ref_unique = [ mol.split('_')[0] for mol in ref_unique] #remove idx
        return shells, ref_unique
    def dataframe(self):
        #extract occurrence of each molecular type from the shell array
        #outputs a pandas dataframe with RM in the rows and molecular types in columns

        # create dictionary to fill with shells
        dict_mol = list(set(i for j in self.shells for i in j))
        dict_mol = dict(zip(dict_mol, np.zeros(len(dict_mol), dtype='int')))
        dataframe = []
        # fill dict with counts of each molecule type
        for shell in self.shells:
            line = dict_mol.copy()
            keys, counts = np.unique(shell, return_counts=True)
            for idx, key in enumerate(keys):
                line[key] = counts[idx]
            dataframe.append(line)
        dataframe = pd.DataFrame(dataframe, index=self.ref_uniq)
        return dataframe


data_obj = Analyzer(filename)
df = data_obj.df
index = set(df.index.values.tolist())

#Select each RF type one at time
for name in index:
    tmp_df = df[df.index == name]
    counts = tmp_df.groupby(list(df.columns), as_index=False).size()
    print(counts)
    #counts = tmp_df.sort_values("size",axis=1, ascending=False)
    #counts["p"] = counts["size"] / sum(counts["size"]) * 100
    print(counts["size"].max())


