import os

from magphylogeny.model.base import BaseModel


class GTDBTkInfer(BaseModel):

    def __init__(self, config):
        super().__init__(config)

    def dir(self, domain):
        return os.path.join(self.config['gtdbtk']['root_dir'], self.config["gtdbtk"]["infer"]["dir"], domain)

    def tree(self, domain):
        return os.path.join(self.dir(domain), 'infer', 'intermediate_results', 'gtdbtk.unrooted.tree')
