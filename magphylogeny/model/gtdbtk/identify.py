import os

from magphylogeny.model.base import BaseModel


class GTDBTkIdentify(BaseModel):

    def __init__(self, config):
        super().__init__(config)
        self.dir = os.path.join(self.config['gtdbtk']['root_dir'], config["gtdbtk"]["identify"]["dir"])

    def markers_summary(self, domain):
        return os.path.join(self.dir, 'identify', f'gtdbtk.{domain}.markers_summary.tsv')
