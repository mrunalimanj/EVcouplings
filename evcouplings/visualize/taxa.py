"""
Taxonomy diversity visualization

Authors:
"""

import pandas as pd
import plotly.express as px

from ete3 import NCBITaxa # another import? can I do this with plotly


def load_taxonomy_lineage(tax_ids):
	"""
	Using NCBITaxa, querying all the taxonomic information on the 
	species' proteins included in the alignment.


    Parameters
    ----------
    tax_ids : Python list
        1D list of NCBI Taxonomy IDs.


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

    rank_sequencevalue_hm = { #will become part of the dataframe
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
    return pd.DataFrame.from_dict(rank_sequencevalue_hm)


def get_taxa(annotation):
	"""
	Helper function for loading taxa from an DataFrame of annotations.


    Parameters
    ----------
    tax_ids : Python list
        1D list of NCBI Taxonomy IDs.


    Returns
    -------
    annotation: pd.DataFrame
    	Original annotations alignment but with taxanomic 
    	information for each entry of the alignment. 
    """
    
    annotation['tax_ID'] = annotation['Tax'].str.split('=').str[-1]
    ncbi = NCBITaxa()

    taxs = loadTaxonomyLineage(annotation['tax_ID'].unique().tolist())
    
    taxs = pd.concat([pd.Series(annotation['tax_ID'].unique(), name='tax_ID'),
                  taxs], axis=1)
    annotation = annotation.merge(taxs, on='tax_ID',how='left')
    
    return annotation



def sunburst(annotation, title, hier=['superkingdom', "phylum", "order"]):
	"""
	Generates sunburst plot from annotation dataframe.

    Parameters
    ----------
    annotation: pd.DataFrame
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
        ['superkingdom', "phylum","genus","class", "subphylum", "family", "order"]

    Returns
    -------
    annotation: pd.DataFrame
    	Original annotations alignment DataFrame but with taxanomic 
    	information for each entry of the alignment. 
    """
    
    # plotly will throw an error if any intermediate rank entries are empty, so 
    # we must filter the annotation DataFrame to include only rows which 
    # have entries for all the intermediate ranks. 
    for rank in hier[:-1]:
        annotation = annotation[(annotation[rank].notnull())]
    
    fig = px.sunburst(annotation, path=hier, values='n', 
    				  title = title, color="superkingdom",
                      color_discrete_map={'Bacteria':'#3366CC', 
                      					  'Eukaryota':'#EF553B', 
                      					  'Archaea':'#FF6692', 
                      					  'Viruses': '#AB63FA'})
    
    return fig
       


