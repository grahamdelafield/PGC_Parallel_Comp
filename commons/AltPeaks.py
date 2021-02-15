import altair as alt
alt.data_transformers.disable_max_rows()
import pandas as pd 
import numpy as np
import os
import re
from venn import venn, pseudovenn

def get_files(directory='.', exts=['-peptides.csv']):
    all_files = []
    for root, _, files in os.walk(directory, topdown=True):
        for name in files:
            file_path = os.path.join(root, name)
            for ext in exts:
                if file_path.endswith(ext):
                    all_files.append(file_path)
    return all_files

class PeaksGroup:
    '''
    Class instantiated from PEAKS DB results files.
    '''

    def __init__(self, file_list):
        assert isinstance(file_list, list), 'Class must be initiated with list of file names'
        self.files = file_list
        self._name_files()

        self.total_peptides = []
        self.unique_peptides = []
        self.peptide_dict = {}
        self._get_peptides()

        self.total_proteins = []
        self.unique_proteins = []
        self.protein_dict = {}
        self._get_proteins()


    def __repr__(self):
        files = '\n'.join(self.files)
        return f'PEAKS group comprised of the files:\n{files}'


    def _name_files(self):
        '''
        Get input from user to name each file.
        '''
        labels = []
        for file in self.files:
            name = input(f'Please name the file {file}:\t')
            labels.append(name)
        self.labels = labels
        return

    def _get_peptides(self):
        '''
        Parse each file and return total/unique proteins
        from each sample.
        '''
        for i, file in enumerate(self.files):
            df = pd.read_csv(file)
            peptides = df.Peptide
            self.total_peptides.append(len(peptides))
            peptides = df[df.Unique=='Y'].Peptide
            self.unique_peptides.append(len(set(peptides)))
            self.peptide_dict[self.labels[i]] = set(peptides)
        print('Peptides parsed successfully.\n')
        return 

    def _get_proteins(self):
        '''
        Parse each file and return total/unique proteins
        from each sample.
        '''
        for i, file in enumerate(self.files):
            df = pd.read_csv(file)
            total_prot = df['Protein Accession']
            self.total_proteins.append(len(set(total_prot)))
            df = df[df.Unique == 'Y'].drop_duplicates(subset='Protein Accession')
            proteins = df['Protein Accession']
            self.unique_proteins.append(len(set(proteins)))
            self.protein_dict[self.labels[i]] = set(proteins)
        print('Proteins parsed successfully.\n')
        return

    def _compile_dataframe(self, kind):
        if kind == 'peptides':
            data = self.total_peptides + self.unique_peptides
        
        elif kind == 'proteins':
            data = self.total_proteins + self.unique_proteins

        labels = self.labels * 2
        kind = ['Total'] * len(self.labels) + ['Unique'] * len(self.labels)

        df = pd.DataFrame({
            'Sample':labels, 'Type':kind, 'Value':data
        })
        return df


    def plot_peptides(self, save=False):
        '''
        Plot barchart of total and unique peptides.
        '''

        df = self._compile_dataframe('peptides')
        bars = alt.Chart(
            df
        ).mark_bar().encode(
            x=alt.X('Type:O', title=None), 
            y=alt.Y('Value:Q', title='Peptide Count'),
            column=alt.Column('Sample:O', title=None),
            color='Type:O'
        ).configure_title(anchor='middle')

        if save:
            bars.save('Peptide_Bar.png', scale_factor=10)
            bars.save('Peptide_Bar.svg', scale_factor=10)
            return

        return bars

  
    def plot_proteins(self, save=False):
        '''
        Plot barchart of total and unique proteins.
        '''

        df = self._compile_dataframe('proteins')
        bars = alt.Chart(
            df, title='Total and Unique Proteins'
        ).mark_bar().encode(
            x=alt.X('Type:O', title=None), 
            y=alt.Y('Value:Q', title='Protein Count'),
            column=alt.Column('Sample:O', title=None),
            color='Type:O', 
        ).configure_title(anchor='middle')

        if save:
            bars.save('Protein_Bar.png', scale_factor=10)
            bars.save('Protien_Bar.svg', scale_factor=10)
            return

        return bars

    def plot_overlay_dist(self, save=False):
        '''
        Plot overlayed lines showing binned m/z values for identified
        peptides.
        '''
        df = pd.DataFrame()
        for i, file in enumerate(self.files):
            bins = {}
            mzs = pd.read_csv(file)['m/z'].values
            binned = mzs//10 * 10
            for b in binned:
                bins[b] = bins.get(b, 0) + 1
            bins, vals = list(bins.keys()), list(bins.values())
            samples = [self.labels[i]] * len(bins)
            data = pd.DataFrame({'Sample': samples, 'Mass':bins, 'Value':vals})
            if df.empty:
                df = data
            else:
                df = pd.concat([df, data])
        
        lines = alt.Chart(df).mark_line().encode(
            x='Mass:Q', y='Value',
            color='Sample'
        ).properties(
            width=800
        )
        if save:
            lines.save('OverlayDist_Line.png', scale_factor=15)
            lines.save('OverlayDist_Line.svg', scale_factor=15)
            return
        return lines

    def plot_all(self, save=False):
        self.plot_peptides(save=save)
        self.plot_proteins(save=save)
        self.plot_overlay_dist(save=save)

if __name__=='__main__':
    files = get_files()
    p = PeaksGroup(files)
    p.plot_proteins(save=True)


