import matplotlib.pyplot as plt
import altair as alt
alt.data_transformers.disable_max_rows()
import pandas as pd 
import numpy as np
import os
from scipy.ndimage import gaussian_filter


def smooth_chrom(xs=[], ys=[], smooth_factor=1, source=None, filename=None, save_as=None):
    '''
    Function to smooth chromatogram from MS data.
    Values can either be read from 
    
    :param xs: (array-like) x values from chromatogram
    :param ys: (array-like) y values from chromatogram
    :param source: (string) dictates from where the function should read data.
        None --> implies data is passed as an argument
        'clip' --> implies pandas should read data from clipboard
        'excel' --> implies pandas should read data from excel file 
    '''
    if source is None:
        asrt_text = 'If no source is provided, data arrays must be passed as arguments'
        assert xs != [] and ys != [], asrt_text
    elif source=='clip':
        df = pd.read_clipboard()
    elif source=='excel':
        df = pd.read_excel(filename)
    else:
        assert xs is not None, 'Found no valid X values'
        assert ys is not None, 'Found no valid Y values'
    xs = df.iloc[3:, 0].astype(float)
    ys = df.iloc[3:, 1].astype(float)
    ys = gaussian_filter(ys, smooth_factor)
    plt.plot(xs, ys)
    plt.fill_between(xs, ys, alpha=0.3)
    if save_as:
        plt.savefig(save_as)

def get_files(directory='.', exts=['.']):
    '''
    Function that searches the defined directory and reutrns list
    of all files with the specified extension or ending.

    :param directory: (str) raw string of directory to be searched
    :param exts: (list) list of extensions or endings to be returned
    '''
    all_files = []
    for root, _, files in os.walk(directory, topdown=True):
        if exts == ['.']:
            all_files.extend([os.path.join(root, name) for name in files])
        for name in files:
            file_path = os.path.join(root, name)
            for ext in exts:
                if file_path.endswith(ext):
                    all_files.append(file_path)
    return all_files

def find_nearest(array, value):
    '''
    Function that searches an array and returns the value nearest to the
    one passed.

    :param array: (array-like) array to be searched
    :param value: (int or float) experimental value
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def plot_ms2_data(xs, ys, peptide, frag_dict, mods=None, show_error=False, tolerance=25):
    '''
    Function to return altair plot of identified fragments for a theoretical 
    peptide.

    :param xs: (array) x/time data
    :param ys: (array) intensity data
    :param peptide: (string) peptide sequence
    :param frag_dict: (dict) output returned from data_processing.fragments func
    '''
    df = pd.DataFrame({
        'x':xs, 'y':ys,
        'fragment':['None']*len(xs),
        'label':['']*len(xs)
    })
    df.loc[:, 'y'] = df.y / np.max(df.y) * 100
    df['label position'] = df.y + 5

    err_mass, err_dist, err_kind= [], [], []

    if mods is not None:
        assert isinstance(mods, dict), 'modifications must enter as dictionary'
        for k in mods:
            frag_dict[k] = mods[k]

    for k, v in frag_dict.items():
        for frag in v:
            nearest = find_nearest(df.x, frag)
            error = mass_error(frag, nearest)
            if abs(error) <= tolerance:
                err_mass.append(nearest)
                err_dist.append(error)
                err_kind.append(k)
                df.loc[(df.x==nearest), 'fragment'] = k
                df.loc[(df.x==nearest), 'label'] = k+f'{v.index(frag)+1}'
    
    dom = {
        'b': '#3b5bad',
        'y': '#c42e23'
    }
    
    bars = alt.Chart(df).mark_bar(size=2).encode(
        x=alt.X('x', title='m/z', axis=alt.Axis(grid=False)),
        y=alt.Y('y', title='Relative Abundance',
                axis=alt.Axis(grid=False, tickCount=1),
                scale=alt.Scale(domain=(0, 100))),
        color=alt.Color('fragment', scale=alt.Scale(domain=list(dom.keys()),
                        range=list(dom.values())), legend=None)
    ).properties(
        title=peptide
    )

    text = alt.Chart(df).mark_text().encode(
        y=alt.Y('label position'),
        x=alt.X('x'),
        text='label'
    )
    chart = alt.vconcat()
    chart &= alt.layer(bars, text)
    if show_error:
        err_df = pd.DataFrame({
            'mass':err_mass,
            'error':err_dist,
            'kind':err_kind
        })
        dots = alt.Chart(err_df).mark_circle().encode(
            x=alt.X('mass:Q', title='m/z', axis=alt.Axis(grid=False)),
            y=alt.Y('error:Q', title='error (ppm)', axis=alt.Axis(grid=True, tickCount=3),
                    scale=alt.Scale(domain=(-tolerance, tolerance))),
            color=alt.Color('kind:O', scale=alt.Scale(domain=list(dom.keys()), range=list(dom.values()))
            )).properties(height=100)

        line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule(
            strokeDash=[10, 10]).encode(y='y')


        chart &= (dots + line)
    return chart.configure_view(
        strokeWidth=0
    ).resolve_scale(x='shared')

def mass_error(measured, exact):
    '''
    Returns mass error between measured and theoretical values.
    '''
    dif = measured - exact
    quo = dif / exact
    return quo * 10**6