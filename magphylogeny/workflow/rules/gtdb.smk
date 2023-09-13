from pathlib import Path
import sys

sys.path.insert(0,Path(workflow.basedir).parent.parent.as_posix())
from magphylogeny.model.gtdb.taxonomy import GTDBTaxonomyFile


tax_file = GTDBTaxonomyFile(config)

rule download_gtdb_taxfile:
    resources:
        queue='quick'
    # benchmark:
    #     "benchmarks/gtdb/download_gtdb_taxfile_{domain}.tsv"
    input:
        tsv=tax_file.get_input_tsv('{domain}')
    output:
        tsv=tax_file.get_output_tsv('{domain}')
    shell:
        """
        mv {input.tsv} {output.tsv}
        """
