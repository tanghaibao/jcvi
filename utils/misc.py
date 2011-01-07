"""
Random ad-hoc functions
"""

# helper functions in the BLAST filtering to get rid alternative splicings
def gene_name(st):
    # this is ugly, but different annotation groups are inconsistent
    # with how the alternative splicings are named;
    # mostly it can be done by removing the suffix
    # except for papaya (evm...) and maize (somewhat complicated)
    if st.startswith("ev"):
        return st
    if st.startswith("Os"):
        return st.rsplit("-",1)[0]
    return st.rsplit(".", 1)[0]
