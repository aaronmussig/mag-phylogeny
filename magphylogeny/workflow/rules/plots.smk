import sys
from pathlib import Path


sys.path.insert(0,Path(workflow.basedir).parent.parent.as_posix())

from magphylogeny.model.gtdbtk.identify import GTDBTkIdentify
from magphylogeny.model.gtdbtk.decorate import GTDBTkDecorate
from magphylogeny.model.plots.marker_summary import PlotMarkerSummary
from magphylogeny.model.plots.red_distribution import PlotRedDistribution
from magphylogeny.model.gtdbtk.phylorank import GTDBTkPhylorankOutliers


GTDBTK_IDENTIFY = GTDBTkIdentify(config)
GTDBTK_DECORATE = GTDBTkDecorate(config)
PLOT_MARKER_SUMMARY = PlotMarkerSummary(config)
PLOT_RED_DIST = PlotRedDistribution(config)
GTDBTK_PHYLORANK_OUTLIERS = GTDBTkPhylorankOutliers(config)

rule plot_marker_summary:
    input:
        markers=GTDBTK_IDENTIFY.markers_summary('{domain}'),
        taxonomy=GTDBTK_DECORATE.taxonomy('{domain}'),
    resources:
        queue='quick'
    output:
        plot=PLOT_MARKER_SUMMARY.path('{domain}','{ext}'),
    shell:
        """
        python "{workflow.basedir}/scripts/plot_marker_summary.py" "{input.markers}" "{input.taxonomy}" "{output.plot}"
        """

rule plot_red_distribution:
    input:
        tree=GTDBTK_PHYLORANK_OUTLIERS.red_scaled_tree('{domain}'),
        red=GTDBTK_PHYLORANK_OUTLIERS.node_rd('{domain}'),
        red_html=GTDBTK_PHYLORANK_OUTLIERS.red_html_plot('{domain}'),
    resources:
        queue='quick'
    output:
        plot=PLOT_RED_DIST.path('{domain}','{ext}'),
    params:
        steps=1000
    shell:
        """
        python "{workflow.basedir}/scripts/plot_novelty_at_red_ranks.py" "{input.tree}" "{input.red}" "{input.red_html}" "{output.plot}" {params.steps}
        """
