"""
Taxonomy diversity visualization

Authors:
  ???? (load_taxonomy_lineage)
  Mrunali Manjrekar
"""

import pandas as pd
import plotly.express as px

from ete3 import NCBITaxa # another import? can I do this with biopython, which might be lighter
# makes more sense to have the database since you're doing thousands of lookups, though
# TODO: see if BioPython has similar capability/ is just as reliable
# TODO: see if timing is similar for code. How does it scale?

COLOR_DISCRETE_MAP =    {'Bacteria':'#56B4E9',  # blue
                        'Eukaryota':'#D53500',  # red
                        'Archaea':'#E69F00',    # orange
                        'Viruses': '#AB63FA',   # purple
                        'Other': '#189e3c'}     # green

SUNBURST_HIERARCHY = ['superkingdom', "phylum", "order"]

PATH_TO_NCBI_TAXA_DATABASE = "/n/groups/marks/databases/etetoolkit/taxa.sqlite" # TODO: add to O2


def load_taxonomy_lineage(tax_ids, ncbi):
    """
    Using NCBITaxa, querying all the taxonomic information on the 
    species' proteins included in the alignment.


    Parameters
    ----------
    tax_ids : Python list
        1D list of NCBI Taxonomy IDs.

    ncbi : NCBITaxa() instance
        An instance that only gets created if get_taxa gets called; that is, the
        the user wants to query the taxonomic ranks for a set of sequences 
        to visualize species diversity or for other purposes. 


    Returns
    -------
    rank_sequencevalue_hm : pd.DataFrame
        dataframe with columns making up all taxonomic ranks
        covered by NCBI. These ranks are: 
        'superkingdom': 
        'phylum': 
        'genus': 
        'class':
        'subphylum': 
        'family': 
        'order': 
        'species': 

    """
    # TODO: update the docstring

    rank_sequencevalue_hm = { # columns of a DataFrame that will get added to the annotation.
        'superkingdom': [],
        'phylum': [],
        'genus': [],
        'class': [],
        'subphylum': [],
        'family': [],
        'order': [],
        'species': [],
    }
    for tax_id in tax_ids:
        try:
            lineage = ncbi.get_lineage(int(tax_id))
            lineageid_name_dict = ncbi.get_taxid_translator(lineage) 
            # dict: key=lineageid, value=sequence value
            
            lineageid_rank_dict = ncbi.get_rank(lineage)
            # flipping the keys and entries of the dictionary.
            rank_lineageid_dict = dict((v,k) for k,v in lineageid_rank_dict.items())

            for rank in rank_sequencevalue_hm:
                sequence_value_for_rank = None
                if rank in rank_lineageid_dict:
                    lineageid = rank_lineageid_dict[rank]
                    sequence_value_for_rank = lineageid_name_dict[lineageid]
                rank_sequencevalue_hm[rank].append(
                    sequence_value_for_rank
                )
        except ValueError as e:
            print('Warning: {0}'.format(str(e)))
            # TODO: consider whether you should adjust this depending on database type.
            # TODO: create test cases? hm
    return pd.DataFrame.from_dict(rank_sequencevalue_hm)


def get_taxa(annotation, format, database_file=PATH_TO_NCBI_TAXA_DATABASE):
    """
    Helper function for loading taxa from an DataFrame of annotations.


    Parameters
    ----------
    tax_ids : Python list
        1D list of NCBI Taxonomy IDs.
    format : {"fasta", "stockholm", None}
        Format of alignment, None if not detectable
        if stockholm, annotation file should include taxids
        if fasta: #TODO: what happens when it's fasta? is
            file difference enough to figure out whether or not taxids will be included
    database_file : string
        TODO: # dbfile="/path/to/taxa.sqlite" 

    Returns
    -------
    annotation: pd.DataFrame
        Original annotations alignment but with taxanomic 
        information for each entry of the alignment. 
    """
    if format != "stockholm":
        pass
        # TODO: Fix the formatting

    annotation['tax_ID'] = annotation['Tax'].str.split('=').str[-1]
    # TODO: check whether these columns are the same for 
    # uniprot vs uniref vs metagenomics: 
    # do the Tax IDs always get represented in this format?

    # "Extract Uniprot/Uniref sequence annotation from Stockholm file
    #(as output by jackhmmer). This function may not work for other
    #formats."

    # TODO: try/except for error checks - are these all numbers, basically? etc 
    # TODO: have column name be pulled properly based on the type of pipeline being used
    # pull in that name properly so as to be consistent with database used and 
    # other config settings, based on `extract_header_annotation` function.

    ncbi = NCBITaxa(taxdump=database_file)  
    # TODO: figure out where this ends up getting downloaded, and
    # make sure it downloads once! can make it similar to SIFTS.py
    # TODO: set path with global variable.

    # TODO: should this be integrated w/ update_database.py?
    # doesn't have to be called in there but can be stored there

    tax_ids = annotation['tax_ID'].unique().tolist()

    taxs = load_taxonomy_lineage(tax_ids, ncbi)
    
    taxs = pd.concat([pd.Series(tax_ids, name='tax_ID'), taxs], axis=1)
    annotation = annotation.merge(taxs, on='tax_ID', how='left')
    
    return annotation


def sunburst(annotation, title, hier=SUNBURST_HIERARCHY, color_map=COLOR_DISCRETE_MAP):

    # keyword argument for hier if confident?
    # other category: very specific colors, you should allow them to pass in colors
    # colormap=...
    # e.g. for contact maps, at the top - global variables w/ color defaults, but can be modified
    # set defaults as global variables/dictionaries

    """
    Generates sunburst plot from annotation dataframe.

    Parameters
    ----------
    annotation : pd.DataFrame
        Original annotations alignment but with taxanomic 
        information for each entry of the alignment. 
    title : Python string 
        Name for plot, to be passed to plotly. 

    hier : Python list
        1D list of desired ranks to include in the plot. 
        These should be ordered from highest to lowest rank desired, 
        and must begin with `superkingdom`, so as to color the top rank
        entries systematically for straightforward comparison between alignments. 
        The ordering of ranks provided by NCBI:
        ['superkingdom', "phylum", "genus", "class", "subphylum", "family", "order"]
    color_map : dictionary
        A mapping to colors for the top rank. 

    Returns
    -------
    fig : Plotly figure
        Sunburst plot instance. Can visualize figure with fig.show() in a Jupyter 
        notebook, or saving to HTML using fig.write_html(filepath) and then 
        opening HTML output in a browser of your choice. 
    """

    # plotly will throw an error if any intermediate rank entries are empty, so 
    # we must fill in the empty intermediate ranks so as not to lose any hits. 
    annotation.fillna("Other", inplace = True) # for purposes of filling in empty intermediate ranks.
    
    annotation["count"] = 1 # a helper column for providing "counts" to plotting function.
    # my alignments were comparatively small. what happens when it gets big?
    # TODO: try running ur code on bigger alignments (think 1-2 GB).
    
    fig = px.sunburst(annotation, path=hier, values='count', 
                      title = title, color=hier[0],
                      color_discrete_map=color_map) 
                                        
    return fig
       


