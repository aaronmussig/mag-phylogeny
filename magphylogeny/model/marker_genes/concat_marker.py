import os

from magphylogeny.model.base import BaseModel


class MarkerGenesConcatMarker(BaseModel):

    def __init__(self, config: dict):
        super().__init__(config)

        self.dir = os.path.join(
            self.config["marker_genes"]["dir"],
            f'concat_markers'
        )

    def marker(self, domain, marker):
        return os.path.join(self.dir, f'{domain}_{marker}.faa')


class MarkerGenesConcatMarkerTrimmed(BaseModel):

    def __init__(self, config: dict):
        super().__init__(config)

        self.dir = os.path.join(
            self.config["marker_genes"]["dir"],
            f'concat_markers_trimmed'
        )

    def marker(self, domain, marker):
        return os.path.join(self.dir, f'{domain}_{marker}.faa')


class MarkerGenesConcatGeneTree(BaseModel):

    def __init__(self, config: dict):
        super().__init__(config)

        self.dir = os.path.join(
            self.config["marker_genes"]["dir"],
            f'gene_trees'
        )

    def tree(self, domain, marker):
        return os.path.join(self.dir, f'{domain}_{marker}.tree')

    def log(self, domain, marker):
        return os.path.join(self.dir, f'{domain}_{marker}.log')

    def mem_mb(self, domain):
        return int(self.config["marker_genes"]["fasttree"]["mem_mb"][domain])

    def runtime(self, domain):
        return int(self.config["marker_genes"]["fasttree"]["runtime_min"][domain])
