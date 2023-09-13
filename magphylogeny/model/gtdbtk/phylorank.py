import os

from magphylogeny.model.base import BaseModel


class GTDBTkPhylorankTrustedTaxa(BaseModel):

    def __init__(self, config):
        super().__init__(config)

        self.dir = os.path.join(self.config['gtdbtk']['root_dir'], self.config["gtdbtk"]["trusted_taxa"]["dir"])

    def trusted_taxa(self, domain):
        return os.path.join(self.dir, f'{domain}_trusted_taxa.txt')


class GTDBTkPhylorankOutliers(BaseModel):

    def __init__(self, config):
        super().__init__(config)

    def dir(self, domain):
        return os.path.join(
            self.config['gtdbtk']['root_dir'],
            self.config["gtdbtk"]["phylorank"]["dir"],
            domain
        )

    def node_rd(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.bootstrap.node_rd.tsv')

    def red_scaled_tree(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.bootstrap.scaled.tree')

    def red_scaled_decorated_tree(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.bootstrap.red_decorated.tree')

    def red_html_plot(self, domain):
        return os.path.join(self.dir(domain), f'gtdbtk.{domain}.rooted.decorated.bootstrap.html')
