import pandas as pd

FILES = {
    "expression": "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "crispr":     "CRISPRGeneEffect.csv",
    "model":      "Model.csv",
}

for label, path in FILES.items():
    df = pd.read_csv(path, index_col=0, nrows=3)
    print("=" * 70)
    print(f"{label.upper()}  ({path})")
    print("index name  :", repr(df.index.name))
    print("index values:", list(df.index))
    print("n columns   :", df.shape[1])
    print("first 5 cols:", list(df.columns[:5]))
    print("last 3 cols :", list(df.columns[-3:]))
    print("dtype counts:", df.dtypes.value_counts().to_dict())
    print("sample row  :", df.iloc[0, :3].tolist())