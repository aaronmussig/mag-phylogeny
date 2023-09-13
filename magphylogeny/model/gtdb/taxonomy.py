import os

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from magphylogeny.model.base import BaseModel

HTTP = HTTPRemoteProvider()


class GTDBTaxonomyFile(BaseModel):

    def __init__(self, config: dict):
        super().__init__(config)

    def url(self, domain):
        return f'data.gtdb.ecogenomic.org/releases/release214/214.1/{domain}_taxonomy_r214.tsv'

    def get_input_tsv(self, domain):
        return HTTP.remote(self.url(domain), keep_local=False)

    def get_output_tsv(self, domain):
        gtdb_dir = self.config['gtdb']['root_dir']
        return os.path.join(gtdb_dir, f'{domain}_taxonomy_r214.tsv')
