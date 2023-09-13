from phylodm import PhyloDM


def get_pdm(tree):
    pdm = PhyloDM.load_from_dendropy(tree)
    dm = pdm.dm(norm=True)
    labels = pdm.taxa()
    return dm, labels
