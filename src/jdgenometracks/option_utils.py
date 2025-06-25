# Option mapping utility for jdgenometracks
from .option_mapping import option_mapping


def set_in_hierarchical_dict(d, subdict, key, value):
    """
    Set value in a nested dict at d[subdict][key] (create subdict as needed).
    If subdict is None, set at top level.
    """
    if subdict is None:
        d[key] = value
    else:
        if subdict not in d:
            d[subdict] = {}
        d[subdict][key] = value


def map_options_for_backend(config: dict, backend: str) -> dict:
    """
    Map unified config keys to hierarchical backend-specific keys using option_mapping.
    Args:
        config (dict): The unified config dict for a track.
        backend (str): 'matplotlib' or 'plotly'.
    Returns:
        dict: A hierarchical dict with backend-specific keys and values.
    """
    mapped = {}
    for key, value in config.items():
        if key in option_mapping and backend in option_mapping[key]:
            subdict, mapped_key = option_mapping[key][backend]
            set_in_hierarchical_dict(mapped, subdict, mapped_key, value)
        else:
            mapped[key] = value
    return mapped
