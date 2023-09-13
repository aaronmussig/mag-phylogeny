import os

from magphylogeny.model.base import BaseModel


class SummaryNoveltyOfGenomes(BaseModel):

    def __init__(self, config: dict):
        super().__init__(config)
        self.dir = config["summary"]["dir"]

    def path(self, domain):
        return os.path.join(self.dir, f'{domain}_novelty_of_genomes.tsv')
