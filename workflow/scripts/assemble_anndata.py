import argparse
import pandas as pd
import anndata as an
import os

def assemble_anndata(counts_dir, output_path, output_txt=None):
    print(f'Assembling AnnData from counts in {counts_dir} to {output_path}')

    # list all files in the directory
    files = os.listdir(counts_dir)
    files = [f for f in files if f.endswith('.txt')]
    print(files)

    # read all files into a dictionary
    data = {}
    var = []
    for f in files:
        df = pd.read_csv(os.path.join(counts_dir, f), sep='\t', index_col=0, low_memory=False)
        counts_cols = [col for col in df.columns if (col.split('.')[-1] in ['count', 'counts'])]
        data[f] = df[counts_cols]
        var.append(df[['Gene_name', 'Gene_Length']])

    # get union of all indices
    all_indices = set()
    for k, v in data.items():
        all_indices = all_indices.union(set(v.index))

    # fill missing values with 0
    for k, v in data.items():
        data[k] = v.reindex(all_indices).fillna(0)

    # concatenate all dataframes
    df = pd.concat(data.values(), axis=1)
    
    def obs_from_col_name(col_name):
        parts = col_name.split('_')
        # indexing from the back, because some samples have additional prefixes
        cell_type = parts[-3]
        treatment = parts[-2]
        replicate = parts[-1].split('.')[0]

        return {
            'cell_type': cell_type,
            'treatment': treatment,
            'replicate': replicate
        }

    obs = pd.DataFrame([obs_from_col_name(col) for col in df.columns])
    obs['label'] = obs.apply(lambda x: 0 if x['treatment'] in ['P', 'EV'] else 1, axis=1)
    obs.index = obs.apply(lambda x: f"{x['cell_type']}_{x['treatment']}_{x['replicate']}", axis=1)

    var = pd.concat(var)
    var = var[~var.index.duplicated(keep='first')]

    adata = an.AnnData(X=df.values.T, obs=obs, var=var)
    adata.write(output_path)

    if output_txt is not None:
        # reindex to gene names using mapping from var
        df.index = var.loc[df.index, 'Gene_name']
        # drop nan indices
        df = df[~df.index.isna()]
        print(f"Keeping {df.shape[0]} genes for txt output")
        df.to_csv(output_txt, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('counts_dir', type=str)
    parser.add_argument('output_h5ad', type=str)
    parser.add_argument('--output_txt', type=str, default=None)
    args = parser.parse_args()

    assemble_anndata(args.counts_dir, args.output_h5ad, args.output_txt)

if __name__ == '__main__':
    main()
