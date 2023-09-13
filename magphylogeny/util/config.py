import os

import yaml

from magphylogeny.util.genome import get_number_of_input_genomes
from magphylogeny.util.path import get_module_dir


def read_config_file():
    config_path = os.path.join(get_module_dir(), 'workflow', 'config', 'config.yaml')
    with open(config_path) as f:
        return yaml.safe_load(f)


def create_config_file(genome_dir: str, output_dir: str, genome_ext: str):
    # Determine the number of genomes in the genome directory
    n_gids = get_number_of_input_genomes(genome_dir, genome_ext)

    # Read the config file yaml
    config = read_config_file()

    # Append the config file with the information read
    config['genome']['dir'] = genome_dir
    config['genome']['ext'] = genome_ext
    config['genome']['qty'] = n_gids

    # Update specific threads with the maximum number of CPUs they should be able to use
    config['gtdbtk']['identify']['cpus'] = min(config['gtdbtk']['identify']['cpus'], n_gids + 1)
    config['gtdbtk']['align']['cpus'] = min(config['gtdbtk']['align']['cpus'], n_gids + 1)

    # Write the changes
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'config.yaml')
    with open(output_path, 'w') as f:
        yaml.dump(config, f)

    return output_path
