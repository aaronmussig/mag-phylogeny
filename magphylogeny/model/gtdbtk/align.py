import os

from snakemake.io import expand

from magphylogeny.model.base import BaseModel


class GTDBTkAlign(BaseModel):

    def __init__(self, config):
        super().__init__(config)
        self.dir = os.path.join(self.config['gtdbtk']['root_dir'], config["gtdbtk"]["align"]["dir"])

    def msa(self, domain):
        return os.path.join(self.dir, 'align', f'gtdbtk.{domain}.msa.fasta.gz')

    def markers(self, domain):
        root = os.path.join(self.dir, 'align', 'intermediate_results', 'markers')
        path = os.path.join(root, f'gtdbtk.{domain}.{{marker}}.faa')
        return expand(path, marker=self.config['gtdbtk']['markers'][domain])

    def marker(self, domain, marker):
        return os.path.join(
            self.dir,
            'align',
            'intermediate_results',
            'markers',
            f'gtdbtk.{domain}.{marker}.faa'
        )
