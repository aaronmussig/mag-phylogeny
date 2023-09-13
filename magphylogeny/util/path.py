import os


def get_module_dir():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    return base_dir


def get_snakemake_file_path():
    return os.path.join(get_module_dir(), 'workflow', "Snakefile")


def get_snakemake_config_path():
    return os.path.join(get_module_dir(), 'workflow', 'config', 'config.yaml')


def get_snakemake_rsrq_path():
    return os.path.join(get_module_dir(), 'workflow', 'rsrq')
