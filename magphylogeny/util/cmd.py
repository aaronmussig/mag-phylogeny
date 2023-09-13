from magphylogeny.util.path import get_snakemake_file_path, get_snakemake_rsrq_path


def get_snakemake_command(output_dir: str, config_path: str):
    snake_path = get_snakemake_file_path()
    cmd = [
        'snakemake',
        '--cores', '90',
        '--snakefile', snake_path,
        '--configfile', config_path,
        '--directory', output_dir,
        '--use-conda',
        '--keep-incomplete'
    ]
    return cmd


def get_snakemake_command_rsrq(output_dir: str, config_path: str):
    snake_path = get_snakemake_file_path()
    rsrq_profile = get_snakemake_rsrq_path()
    cmd = [
        'snakemake',
        '--jobs', '200',
        '--snakefile', f'{snake_path}',
        '--configfile', config_path,
        '--directory', output_dir,
        '--use-conda',
        '--profile', rsrq_profile,
        '--latency-wait', '60',
        '--rerun-incomplete',
        '--max-status-checks-per-second', '5',

        # '--cluster-status', 'rsrq snakemake status',
        # '--cluster', 'rsrq snakemake submit',
        # '--cluster-cancel', 'rsrq snakemake cancel',
    ]
    print(' '.join(cmd))
    return cmd
