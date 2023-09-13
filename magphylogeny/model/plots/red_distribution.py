import os

from magphylogeny.model.base import BaseModel


class PlotRedDistribution(BaseModel):

    def __init__(self, config: dict):
        super().__init__(config)
        self.dir = config["plots"]["dir"]

    def path(self, domain, ext):
        return os.path.join(self.dir, f'{domain}_red_distribution.{ext}')
