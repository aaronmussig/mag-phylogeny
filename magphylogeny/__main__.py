import subprocess

import typer

from magphylogeny.util.cmd import get_snakemake_command, get_snakemake_command_rsrq
from magphylogeny.util.config import create_config_file
from magphylogeny.util.path import get_module_dir, get_snakemake_file_path, get_snakemake_config_path

app = typer.Typer()


# @app.command()
# def bootstrap_msa():
#     cmd = get_snakemake_command()
#     cmd.extend(('-r', 'all_bootstrap_msa_ar53'))
#     try:
#         subprocess.run(cmd, check=True)
#     except Exception as e:
#         print(f'snakemake failed: {e}')

@app.command()
def local(genome_dir: str, output_dir: str, genome_ext: str):
    """Run all steps on the local machine."""
    config_path = create_config_file(genome_dir, output_dir, genome_ext)
    cmd = get_snakemake_command(output_dir, config_path)
    try:
        subprocess.run(cmd, check=True)
        print("Done")
    except Exception as e:
        print(f'snakemake failed: {e}')


@app.command()
def rsrq(genome_dir: str, output_dir: str, genome_ext: str):
    """Run all steps using RSRQ as the cluster scheduler."""
    config_path = create_config_file(genome_dir, output_dir, genome_ext)
    cmd = get_snakemake_command_rsrq(output_dir, config_path)
    try:
        subprocess.run(cmd, check=True)
        print("Done")
    except Exception as e:
        print(f'snakemake failed: {e}')


# @app.command()
# def plot():
#     cmd = [
#         'snakemake',
#         '--cores', '1',
#         '--snakefile', get_snakemake_file_path(),
#         '--configfile', get_snakemake_config_path(),
#         '--directory', '/tmp/aaron/pipeline',
#         '--dag'
#     ]
#     try:
#         subprocess.run(cmd, check=True)
#     except Exception as e:
#         print(f'snakemake failed: {e}')

# https://github.com/akcorut/kGWASflow/blob/bd28c56ad92a711e51a527501e28356ed291c803/workflow/rules/common.smk#L317
def main():
    app()


if __name__ == '__main__':
    main()
