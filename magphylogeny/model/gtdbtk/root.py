import os

from magphylogeny.model.base import BaseModel


class GTDBTkRoot(BaseModel):

    def __init__(self, config):
        super().__init__(config)

    def dir(self, domain):
        return os.path.join(self.config['gtdbtk']['root_dir'], self.config["gtdbtk"]["root"]["dir"], domain)

    def tree(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.tree')

    def outgroup(self, domain):
        return self.config['gtdbtk']['root']['outgroup'][domain]
