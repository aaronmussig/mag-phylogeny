import os

from snakemake.io import expand

from magphylogeny.model.base import BaseModel


class GTDBTkBootstrapMsaReplicate(BaseModel):

    def __init__(self, config):
        super().__init__(config)

    def dir(self, domain):
        return os.path.join(
            self.config['gtdbtk']['root_dir'],
            self.config["gtdbtk"]["bootstrap"]["dir"],
            'msa',
            domain
        )

    def faa(self, domain, replicate):
        return os.path.join(self.dir(domain), f'{domain}_{replicate}.faa')


class GTDBTkBootstrapTreeReplicate(BaseModel):

    def __init__(self, config):
        super().__init__(config)

    def dir(self, domain):
        return os.path.join(
            self.config['gtdbtk']['root_dir'],
            self.config["gtdbtk"]["bootstrap"]["dir"],
            'tree',
            domain
        )

    def tree(self, domain, replicate):
        return os.path.join(self.dir(domain), f'{domain}_{replicate}.tree')

    def trees(self, domain):
        return expand(
            self.tree(domain, '{replicate}'),
            replicate=range(self.config["gtdbtk"]["bootstrap"]["n_boot"])
        )

    def log(self, domain, replicate):
        return os.path.join(self.dir(domain), f'{domain}_{replicate}.log')


class GTDBTkBootstrapMergeTrees(BaseModel):

    def __init__(self, config):
        super().__init__(config)

    def dir(self, domain):
        return os.path.join(
            self.config['gtdbtk']['root_dir'],
            self.config["gtdbtk"]["bootstrap"]["dir"],
            'tree_merged',
            domain
        )

    def tree(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.bootstrap.tree')
