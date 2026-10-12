"""Save a research output (figure, table or text) to a file. One small function, used by every research script.

    from research.common.save import save
    save(fig, P.FIG_MEDIAN_CURVES, dpi=110)          # matplotlib figure -> savefig
    save(df, P.TAILS_INDEX, index=False)             # DataFrame -> .csv / .csv.gz, or .parquet by the file name
    save(report_text, P.REPORT)                      # str -> UTF-8 text file

The folder is created if needed. Extra keyword arguments go to the writer (dpi=..., index=False, compression=...).
"""
from pathlib import Path


def save(obj, path, **kwargs):
    """Write `obj` to `path` and return the path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(obj, "savefig"):                                          # matplotlib figure
        obj.savefig(path, **kwargs)
    elif hasattr(obj, "to_csv"):                                         # pandas DataFrame / Series
        if path.suffix == ".parquet":
            obj.to_parquet(path, **kwargs)
        else:
            obj.to_csv(path, **kwargs)
    elif type(obj).__module__.startswith("pyarrow"):                     # pyarrow Table
        import pyarrow.parquet as pq
        pq.write_table(obj, path, **kwargs)
    elif isinstance(obj, str):
        path.write_text(obj, encoding=kwargs.get("encoding", "utf-8"))
    else:
        raise TypeError(f"save: cannot write a {type(obj).__name__} to {path}")
    return path
