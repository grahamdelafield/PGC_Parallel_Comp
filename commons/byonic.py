import os
import re
import pandas as pd


class ByFile:
    """
    A class representing an input Byonic File.

    Functions built enable, combining multiple bf files,
    renaming columns in byonic files, cleaning the peptide sequence,
    determining glycosites, filtering results, and pulling total/unique
    glycopeptides.

    Typical usage:
        bf = ByFile(data)
        bf.fill_no_glycans()
        bf.remove_reverse(modify=?)
        bf.determine_glycosites()
        bf.categorize_glycans()
        bf.filter_hits(modify=?)
        bf.rame = bf.reduce_frame(gp_only=?)
        total = bf.total_gp()
        unique = bf.unique_gp()
    
    """

    def __init__(self, file_input):
        if isinstance(file_input, list) and len(file_input) > 1:
            self.frame = self.combine_dataframes(file_input)
        elif isinstance(file_input, list) and len(file_input) == 1:
            self.frame = pd.read_excel(file_input[0], "Spectra")
        else:
            self.frame = pd.read_excel(file_input, "Spectra")

        cols = [c.lower() for c in self.frame]
        cols = [c.replace(" ", "_") for c in cols]
        cols = ["_".join(c.split("\n")) for c in cols]
        pattern = r"[\(\)\+\:\/\|]"
        for i in range(len(cols)):
            cols[i] = re.sub(pattern, "", cols[i])
            cols[2:5] = ["peptide", "glycan", "modifications"]
        self.frame.columns = cols

        self.clean_peptides()
        self.fill_no_glycans()

    def combine_dataframes(self, file_list):
        """
        Concatenates all DataFrame objects passed in "df_list."

        attributes:
        df_list (type: list) list of pd.DataFrame objects
        """
        assert isinstance(file_list, list)

        df = pd.read_excel(file_list[0], "Spectra")
        i = 1
        while i < len(file_list):
            df2 = pd.read_excel(file_list[i], "Spectra")
            df = pd.concat([df, df2], ignore_index=True)
            i += 1

        return df

    def rename_columns(self, old_name, new_name):
        cols = [c for c in self.frame.columns]
        cols[cols.index(old_name)] = new_name
        self.frame.columns = cols

    def clean_peptides(self):
        peptides = self.frame.peptide.tolist()
        pattern = [r"^[A-Z]\.", r"\.[A-Z]$", r"\[\+\d*\.\d*\]", r"^\-\.", r"\.\-$"]
        for i in range(len(peptides)):
            for p in pattern:
                peptides[i] = re.sub(p, "", peptides[i])
        self.frame["clean_peptide"] = peptides

    def fill_no_glycans(self):
        new_list = self.frame.glycan.fillna(0)
        self.frame.loc[:, "glycan"] = new_list

    def remove_reverse(self, modify=False):
        pattern = r"Reverse"
        idxs = []
        for row in self.frame.index:
            cell = self.frame.loc[row, 'protein_name']
            if re.search(pattern, cell):
                idxs.append(row)
        if modify:
            self.frame.drop(idxs, inplace=True)
            return self.frame
        return self.frame.drop(idxs)

    def determine_glycosites(self):
        sub = self.frame[self.frame.glycan != 0]
        idxs = sub.index
        peps = sub.peptide.tolist()
        start_pos = sub.starting_position.tolist()
        pattern = [r"\[\+\d*\.\d*\]", r"\."]
        for i in range(len(peps)):
            peps[i] = peps[i].split("N[")[0]
            for p in pattern:
                peps[i] = re.sub(p, "", peps[i])
            peps[i] = len(peps[i])
            peps[i] = start_pos[i]+peps[i]
        self.frame.loc[idxs, "glycosite"] = peps
        self.frame.glycosite = self.frame.glycosite.fillna(0)
        self.frame.glycosite = self.frame.glycosite.astype(int)

        return self.frame

    def get_glycosites(self):
        site_dict = {}
        sub_frame = self.frame[self.frame.glycan != 0][['glycosite', 'protein_name']]
        for row in sub_frame.index:
            protein = sub_frame.loc[row, "protein_name"]
            site = sub_frame.loc[row, "glycosite"]
            if protein not in site_dict:
                site_dict[protein] = []
            if site not in site_dict[protein]:
                site_dict[protein].append(site)
        return site_dict

    def categorize_glycan(self):
        glycans = self.frame.glycan.tolist()
        glycan_types = []
        for s in glycans:
            s = s.replace(')', ',')
            s = s.replace('(', ' ')
            s = s.split(',')[:-1]
            d = {k:int(v) for [k, v] in [i.split(' ') for i in s]}

            if 'NeuAc' in d or 'NeuGc' in d:
                glycan_types.append('Sialylated')
            elif 'Fuc' in d:
                if d['HexNAc'] > 2:
                    glycan_types.append('Fucosylated')
                elif d['HexNAc'] == 2:
                    if 'Hex' in d:
                        if d['Hex'] > 4:
                            glycan_types.append('Complex')
                        else:
                            glycan_types.append('Paucimannose' )
            elif d['HexNAc'] > 2:
                glycan_types.append('Complex')
            if d['HexNAc'] <= 2:
                if 'Hex' in d:
                    if d['Hex'] <= 9 and d['Hex'] > 4:
                        glycan_types.append('High Mannose')
                glycan_types.append('Paucimannose')

        self.frame.loc[:, 'glycan_types'] = glycan_types

    def filter_hits(self, score=150, delta_mod=10, log_prob=1, modify=False):
        crit_1 = (self.frame.score >= score)
        crit_2 = (self.frame.delta_mod >= delta_mod)
        crit_3 = (self.frame.log_prob >= log_prob)

        if modify:
            self.frame = self.frame[crit_1 & crit_2 & crit_3]
        return self.frame[crit_1 & crit_2 & crit_3]

    def total_gp(self, cols=["peptide", "z", "glycan", "calc._mass_mh"]):
        return self.frame.drop_duplicates(subset=cols, keep='first',
                                          inplace=False)

    def unique_gp(self, cols=["peptide", "glycan", "calc._mass_mh"]):
        return self.frame.drop_duplicates(subset=cols, keep='first',
                                          inplace=False)

    def reduce_frame(self, cols=["clean_peptide", "glycan", "z", "observed_mz",
                                 "calc._mass_mh", "glycosite", "score",
                                 "delta_mod", "log_prob", "peptide",
                                 "protein_name"], gp_only=False):
        if gp_only:
            return self.frame[self.frame.glycan != 0][cols]
        return self.frame[cols]


# bf = By_File(f)
# bf.frame.columns
# bf.frame.shape
# bf.determine_glycosites().shape
# bf.remove_reverse(modify=True).shape
# bf.frame.shape
#
#
# bf.filter_hits(modify=True).shape
# print(bf.total_gp().shape)
# bf.unique_gp()
# # bf.frame[['score', 'delta_mod', 'log_prob']].describe()
# bf.frame = bf.reduce_frame(gp_only=True)
# bf.frame.shape
# bf.filter_hits(modify=True).shape
# bf.total_gp().shape
# bf.unique_gp().shape
# bf.frame
