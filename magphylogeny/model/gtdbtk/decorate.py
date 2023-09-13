import os

from magphylogeny.model.base import BaseModel


class GTDBTkDecorate(BaseModel):

    def __init__(self, config):
        super().__init__(config)

    def dir(self, domain):
        return os.path.join(self.config['gtdbtk']['root_dir'], self.config["gtdbtk"]["decorate"]["dir"], domain)

    def tree(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.tree')

    def tree_table(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.tree-table')

    def taxonomy(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.tree-taxonomy')
