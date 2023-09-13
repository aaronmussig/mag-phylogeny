import os

from snakemake.io import expand
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from magphylogeny.model.base import BaseModel

HTTP = HTTPRemoteProvider()


class MarkerGenesRepDownload(BaseModel):

    def __init__(self, config: dict):
        super().__init__(config)

    def url(self, domain):
        return f'data.gtdb.ecogenomic.org/releases/release214/214.1/genomic_files_reps/{domain}_msa_marker_genes_reps_r214.tar.gz'

    def dir(self, domain):
        return os.path.join(
            self.config["marker_genes"]["dir"],
            f'reps_{domain}'
        )

    def archive(self, domain):
        return HTTP.remote(self.url(domain), keep_local=False)

    def markers(self, domain):
        path = os.path.join(self.dir(domain), 'individual', f'{domain}_r214_reps_{{marker}}.faa')
        return expand(path, marker=self.config["gtdbtk"]["markers"][domain])

    def marker(self, domain, marker):
        return os.path.join(self.dir(domain), 'individual', f'{domain}_r214_reps_{marker}.faa')
