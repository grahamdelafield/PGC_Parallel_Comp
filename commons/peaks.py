'''
The following code is an OOP arrangement for working with
PEAKS DB result files. The following representations are
intended to employ the proteins-peptides result files after
database searching. 

Further functionality may be developed upon interest.
Contact: delafield@wisc.edu
'''

import os
import re
import shutil
import pandas as pd
from venn import venn, pseudovenn
import  upsetplot
import matplotlib.pyplot as plt
plt.style.use('seaborn')


def get_files(directory=".", exts=["-peptides.csv"]):
    all_files = []
    for root, dirs, files in os.walk(directory, topdown=True):
        for name in files:
            file_path = os.path.join(root, name)
            for ext in exts:
                if file_path.endswith(ext):
                    all_files.append(file_path)
    return all_files

def file_root(filename):
    filename = filename.split("\\")[:-1]
    filename = "\\".join(filename)
    return filename


def p_diff(v1, v2):
    num = abs(v1 - v2)
    denom = abs(((v1 + v2) / 2))
    return num/denom * 100

    
class peaks_group:

    def __init__(self, name, file_list):
        self.name = name
        self.files = file_list
        self.labels = []
        self.total_peptides = []
        self.unique_peptides = []
        self.peptide_dict = {}
        self.protein_dict = {}
        self.total_proteins = []
        self.unique_proteins = []
        self.prot_cov = []

        self.peptide_files = [f for f in self.files if f.endswith("-peptides.csv")]
        self.protein_files = [f for f in self.files if f.endswith("proteins.csv")]
        
        self._name_files()
        self._get_peptides()
        self._get_proteins()
        self._write_results()
        self._get_coverage()

    def __str__(self):
        files = "\n".join(self.files)
        return f'PEAKS group {self.name} comprised of the files:\n{files}'

    def _name_files(self):
        d = {}
        labels = []
        for _, file in enumerate(self.files):
            root = "\\".join(file.split("\\")[:-1])
            if root not in d:
                name = input(f"Provide label for {root}")
                d[root] = d.get(root, name)
                self.labels.append(name)
        return  

    def _get_peptides(self):
        file_list = self.peptide_files
        labels = self.labels
        for i, file in enumerate(file_list):
            label = labels[i]
            df = pd.read_csv(file)
            peptides = df.Peptide.tolist()
            self.total_peptides.append(len(peptides))
            

            peptides = df[df.Unique == "Y"].Peptide.tolist()
            self.unique_peptides.append(len(set(peptides)))
            self.peptide_dict[label] = set(peptides)
        print("Total peptides imported.")
        print("Unique peptides evaluated.")
        print("Peptide dict initialized.\n")


    def _get_proteins(self):
        file_list = self.peptide_files
        labels = self.labels
        for i, file in enumerate(file_list):
            label = labels[i]
            df = pd.read_csv(file)
            total_prot = df["Protein Accession"].tolist()
            df = df[df.Unique == "Y"].drop_duplicates(subset="Protein Accession")
            proteins = df["Protein Accession"].tolist()

            self.protein_dict[label] = set(proteins)
            self.total_proteins.append(len(set(total_prot)))
            self.unique_proteins.append(len(list(set(proteins))))
        print("Total proteins imported.")
        print("Unique peptides evaluated.")
        print("Protein dict initialized.")

    def _get_coverage(self):
        file_list = self.protein_files
        for file in file_list:
            df = pd.read_csv(file)
            cov = df["Coverage (%)"].values[0]
            self.prot_cov.append(cov)
        return


    def _write_results(self):
        with open("PEAKS_compare_output.txt", "w") as f:
            f.write("Total peptides identified:" + "\n")
            labels = self.labels
            for sample, num in zip(labels, self.total_peptides):
                f.write(str(sample)+": "+str(num)+"\n")
            f.write("\n")
            f.write("Unique peptides identified:" + "\n")
            for sample, num in zip(labels, self.unique_peptides):
                f.write(str(sample)+": "+str(num)+"\n")
            f.write("\n")
            f.write("Total proteins identified:" + "\n")
            for sample, num in zip(labels, self.total_proteins):
                f.write(str(sample)+": "+str(num)+"\n")
            f.write("\n")
            f.write("Unique proteins identified:" + "\n")
            for sample, num in zip(labels, self.unique_proteins):
                f.write(str(sample)+": "+str(num)+"\n")
        print("Results file created.")


    def plot_peptides(self, save=False):
        fig, ax = plt.subplots(figsize=(10, 5))

        plot_df = pd.DataFrame({"Total Peptides": self.total_peptides,
                                "Unique Peptides": self.unique_peptides}, index=self.labels)
        plot_df.plot.bar(rot=0, fontsize=13, ax=ax)
        ax.legend(loc='best', bbox_to_anchor=[1, 1], fontsize=13)
        ax.set_title("Total and Unique Peptides", fontsize=20)
        if save:
            plt.savefig("Total_and_unique_peptides.svg")
            plt.savefig("Total_and_unique_peptides.png", bbox_inches="tight")

    def plot_peptide_overlap(self, save=False):
        fig, ax = plt.subplots(figsize=(10, 10))
        if len(self.peptide_dict) in range(2, 6):
            venn(self.peptide_dict, cmap="viridis", ax=ax)
        elif len(self.peptide_dict) == 6:
            pseudovenn(self.peptide_dict, cmap="viridis", ax=ax)
        else:
            print("No Peptide Venn Diagram plotted due to invalid number of samples.")
            print("Venn Diagrams require between 2 and 6 samples.")

        ax.set_title("Peptide Overlap", fontsize=20)
        ax.legend(labels=self.labels, fontsize=15, loc='best', bbox_to_anchor=[1.1, 1])
        if save:
            fig.savefig("Peptide_overlap_venn.svg")
            fig.savefig("Peptide_overlap_venn.png", bbox_inches="tight")


    def plot_peptide_upset(self, save=False):
        color = '#21918cff'
        plot_df = upsetplot.from_contents(self.peptide_dict)
        upsetplot.plot(plot_df, sort_by='cardinality', subset_size='auto', facecolor=color)
        # plt.ylim(0, 400)
        plt.title("Distribution of Peptide Overlap")
        if save:
            plt.savefig("Peptide_upset.svg")
            plt.savefig("Peptide_upset.png")

    def plot_proteins(self, save=False):
        fig, ax = plt.subplots(figsize=(10, 5))
        plot_df = pd.DataFrame({"Total Proteins": self.total_proteins,
                                "Unique Proteins": self.unique_proteins}, index=self.labels)
        plot_df.plot.bar(rot=0, fontsize=13, legend=None, ax=ax)

        y_max, difference = plt.yticks()[0][-1], plt.yticks()[0][-1]-plt.yticks()[0][-2]
        ax.set_ylim(0, y_max+difference)
        ax.set_title("Total Proteins Identified", fontsize=20)
        if save:
            plt.savefig("Total_proteins.svg")
            plt.savefig("Total_proteins.png")

    def plot_protein_overlap(self, save=False):
        fig, ax = plt.subplots(figsize=(10, 10))

        if len(self.protein_dict) in range(2, 6):
            venn(self.protein_dict, cmap="viridis", ax=ax)
        elif len(self.protein_dict) == 6:
            pseudovenn(self.protein_dict, cmap="viridis", ax=ax)
        else:
            print("No Protein Venn Diagram plotted due to invalid number of samples.")
            print("Venn Diagrams require between 2 and 6 samples.")

        ax.legend(labels=self.labels, fontsize=15, loc="best", bbox_to_anchor=[1.1, 1])
        ax.set_title("Protein Overlap", fontsize=20)
        if save:
            plt.savefig("Protein_overlap_venn.svg")
            plt.savefig("Protein_overlap_venn.png", bbox_inches="tight")


    def plot_protein_upset(self, save=False):
        color = '#21918cff'
        plot_df = upsetplot.from_contents(self.protein_dict)
        upsetplot.plot(plot_df, sort_by='cardinality', subset_size='auto', facecolor=color)
        # plt.ylim(0, 60)
        plt.title("Distribution of Protein Overlap")
        if save:
            plt.savefig("Protein_upset.svg")
            plt.savefig("Protein_upset.png")

    def plot_overlay_dist(self, save=False):
        fig, ax = plt.subplots(figsize=(20, 6))
        colors = [plt.cm.viridis(i/float(len(self.labels)-1)) for i in range(len(self.labels))]


        for i in range(len(self.labels)):
            binned = {}
            observed_mz = []
            df = pd.read_csv(self.peptide_files[i])
            mz = df["m/z"].tolist()
            observed_mz.append(mz)
            observed_mz = list(set(mz))

            for mass in observed_mz:
                bin_val = (mass//10)*10
                if bin_val not in binned:
                    binned[bin_val] = 1
                else:
                    binned[bin_val] += 1

            binned = dict(sorted(binned.items(), key=lambda x: x[0]))
            xs = [x for x in binned.keys()]
            ys = [y for y in binned.values()]
            ax.scatter(xs, ys, color=colors[i])
            ax.plot(xs, ys, linestyle="--", color=colors[i], label=self.labels[i])
            ax.legend(fontsize=15)
            ax.xaxis.set_tick_params(labelsize=13)
            ax.yaxis.set_tick_params(labelsize=13)
            ax.set_title("Overlaid Distribution of detected m/z", fontsize=20)
        if save:
            plt.savefig("Line_overlay_dist.svg")
            plt.savefig("Line_overlay_dist.png")

    def plot_mult_dist(self, save=False):
        fig, axs = plt.subplots(len(self.labels), 1, sharex=True, figsize=(20, 20))
        colors = [plt.cm.viridis(i/float(len(self.labels)-1)) for i in range(len(self.labels))]

        y_max = 0

        for i in range(len(self.labels)):
            binned = {}
            observed_mz = []
            df = pd.read_csv(self.peptide_files[i])
            mz = df["m/z"].tolist()
            observed_mz.append(mz)

            observed_mz = list(set(mz))

            for mass in observed_mz:
                bin_val = (mass//10)*10
                if bin_val not in binned:
                    binned[bin_val] = 1
                else:
                    binned[bin_val] += 1

            binned = dict(sorted(binned.items(), key=lambda x: x[0]))
            if max(binned.values()) > y_max:
                y_max = max(binned.values())//10*10+10

            axs[i].bar(binned.keys(), binned.values(), width=8, color=colors[i])
            axs[i].set_title("Distribution of m/z for sample " + self.labels[i], fontsize=20)
            axs[i].set_ylim(0, y_max)
            axs[i].set_ylabel("Peptide Count", fontsize=12)
            axs[i].yaxis.set_major_locator(plt.MaxNLocator(5))
            axs[i].yaxis.set_tick_params(labelsize=14)

        ax = plt.gca()
        ax.xaxis.set_tick_params(labelsize=14)
        fig.text(0.5, 0.08, 'm/z', ha='center', size=20)
        if save:
            plt.savefig("Multi_dist.svg")
            plt.savefig("Multi_dist.png")

    def plot_retention(self, save=False):
        hists = []
        for i in range(len(self.peptide_files)):
            file = self.peptide_files[i]
            df = pd.read_csv(file)
            rt = df.RT.tolist()

            buckets = {}
            if max(rt) < 20:
                window = 2  # minutes
            window = 5  # minutes
            for item in rt:
                time = (item//window) * window
                if time not in buckets:
                    buckets[time] = 1
                else:
                    buckets[time] += 1

            buckets = dict(sorted(buckets.items(), key=lambda x: x[0], reverse=False))
            hists.append(buckets)

        grad_info = str(input("Would you like to enter gradient information? [y/n]\t")).lower()
        if grad_info == "y":
            grad_x = []
            grad_y = [None]
            while len(grad_x) != len(grad_y):
                grad_x = str(input("Enter time points for gradient (ex. 0, 5, 10...):\t")).split(",")
                grad_x = [float(x) for x in grad_x]
                grad_y = str(input("Enter %B values for gradient (ex. 0, 5, 10...):\t")).split(",")
                grad_y = [float(y) for y in grad_y]
                if len(grad_x) != len(grad_y):
                    print("Error: please make the number of points for time and %B are equal.")
        if len(hists) <= 2:
            n_cols = 2
            n_rows = 1
        else:
            n_cols = 3
            n_rows = len(hists) // n_cols
            if len(hists) > n_cols * n_rows:
                n_rows += 1
        fig, ax = plt.subplots(n_rows, n_cols, figsize=(15, 10), sharey=True, sharex=True)
        colors = [plt.cm.viridis(i/float(25-1), alpha=0.75) for i in range(25)]

        if n_rows == 1:
            i = 0
            for c in range(n_cols):
                x = hists[i].keys()
                y = hists[i].values()
                ax1 = ax[c]
                ax1.bar(x, y, width=window, facecolor=colors[8], edgecolor=colors[8])
                ax1.set_title("Dist of RTs, "+self.labels[i])
                if grad_info == 'y':
                    ax2 = ax1.twinx()
                    ax2.plot(grad_x, grad_y, color=colors[6])
                    ax2.grid(visible=False)
                i += 1

        else:
            i = 0
            for r in range(n_rows):
                for c in range(n_cols):
                    if i > len(hists)-1:
                        break
                    x = hists[i].keys()
                    y = hists[i].values()
                    ax1 = ax[r][c]
                    ax1.bar(x, y, width=5,  facecolor=colors[8], edgecolor=colors[8])
                    ax1.set_title("Dist of RTs, "+self.labels[i])
                    if grad_info == 'y':
                        ax2 = ax1.twinx()
                        ax2.plot(grad_x, grad_y, color=colors[6])
                        ax2.grid(visible=False)
                    i += 1
        if save:
            plt.savefig("RetentionOverlay.svg")
            plt.savefig("RetentionOverlay.png")

    def plot_coverage(self, save=False):
        fig, ax = plt.subplots(figsize=(10, 5))
        plot_df = pd.DataFrame({"Unique Peptides": self.unique_peptides,
                                "Coverage":self.prot_cov}, 
                                index=self.labels)            
        plot_df["Coverage"].plot(secondary_y=True, kind="line", color="k", ax=ax, rot=0)
        # plot_df["Coverage"].plot(secondary_y=True, kind="scatter", color="k", ax=ax, rot=0)
        plt.grid(visible=False)
        plot_df[["Unique Peptides"]].plot.bar(ax=ax, rot=0)
        ax.legend(bbox_to_anchor=[1.2, 1])
        if save:
            plt.savefig("ProteinCoverage.png")
            plt.savefig("ProtienCoverage.svg")


    def plot_all(self, flag=False):
        self.plot_peptides(save=flag)
        self.plot_peptide_overlap(save=flag)
        self.plot_peptide_upset(save=flag)
        self.plot_proteins(save=flag)
        self.plot_protein_overlap(save=flag)
        self.plot_protein_upset(save=flag)
        self.plot_mult_dist(save=flag)
        self.plot_overlay_dist(save=flag)
        self.plot_retention(save=flag)


