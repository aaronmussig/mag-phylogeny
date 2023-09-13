import os


def get_number_of_input_genomes(path: str, extension: str) -> int:
    return len([name for name in os.listdir(path) if name.endswith(extension)])
