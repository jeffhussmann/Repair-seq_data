import logging
import warnings
from pathlib import Path

import matplotlib.colors
import numpy as np
import pandas as pd
from jinja2 import Template

import repair_seq as rs
from repair_seq.visualize.interactive import alphabetical, UMAP
from repair_seq.visualize import gene_significance

import repair_seq.base_editor as be

import hits.visualize.interactive

import templates

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                   )


def make_index():
    logging.info('Starting make_index()')
    
    template = Template(Path('base.jinja').read_text())
    
    data = {
        'navbar': Path('navbar.html').read_text(),
        'contents': Path('about.html').read_text(),
    }
    
    Path('../index.html').write_text(template.render(data))
    
def make_BE_scatters():
    logging.info('Starting make_BE_scatters()')
    
    base_dir = '/lab/solexa_weissman/jah/projects/base_editing_screens'

    pools = rs.pooled_screen.get_all_pools(base_dir)
    
    gl = rs.guide_library.GuideLibrary(base_dir, 'DDR_library')
    
    pns = [
        '2018_12_24_base_editor_BE1_1',
        '2018_12_24_base_editor_BE1_2',
        '2018_12_24_base_editor_BE4_1',
        '2018_12_24_base_editor_BE4_2',
    ]

    pool_and_outcomes = []

    for pn in pns:
        pool = pools[pn]
        pool_and_outcomes.append((pool, ('Any C→T edit', be.at_least_n_Bs(pool, 1, 'T'))))
        pool_and_outcomes.append((pool, ('Any C→G edit', be.at_least_n_Bs(pool, 1, 'G'))))

        pool_and_outcomes.append((pool, ('Exactly 1 C→T edit', be.exactly_n_Bs(pool, 1, 'T'))))
        pool_and_outcomes.append((pool, ('Exactly 2 C→T edits', be.exactly_n_Bs(pool, 2, 'T'))))
        pool_and_outcomes.append((pool, ('Exactly 3 C→T edits', be.exactly_n_Bs(pool, 3, 'T'))))

        pool_and_outcomes.append((pool, ('Deletion', ['deletion'])))
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
    
        data = gene_significance.compute_table(pool_and_outcomes,
                                               initial_dataset='HeLa BE4B rep1',
                                               initial_outcome='Any C→G edit',
                                               pool_name_aliases={
                                                   '2018_12_24_base_editor_BE1_1': 'HeLa BE1 rep1',
                                                   '2018_12_24_base_editor_BE1_2': 'HeLa BE1 rep2',
                                                   '2018_12_24_base_editor_BE4_1': 'HeLa BE4B rep1',
                                                   '2018_12_24_base_editor_BE4_2': 'HeLa BE4B rep2',
                                               },
                                              )
    
    df = data['guides_df'].xs('log2_fold_change', axis='columns', level=1).copy()

    df.columns = pd.MultiIndex.from_tuples(df.columns.map(lambda s: tuple(s.rsplit('_', 1))))

    df.columns.names = ['Screen condition', 'Outcome category']

    df['color'] = 'grey'

    df.loc[gl.gene_guides('RFWD3'), ('color', '')] = 'tab:purple'
    df.loc[gl.gene_guides('UNG'), ('color', '')] = 'tab:orange'
    df.loc[gl.gene_guides('ASCC3'), ('color', '')] = 'tab:green'

    df['color'] = df['color'].map(matplotlib.colors.to_hex)

    df.loc[gl.non_targeting_guides, ('color', '')] = 'black'

    df['gene'] = df.index.map(gl.guide_to_gene_with_non_targeting_guide_sets)

    df.index.name = 'CRISPRi sgRNA'
    
    description_data = {
        'title': 'Scatter plots of log₂ fold changes in base editing outcome frequencies',
        'details': 'Comparisons of the effects of CRISPRi sgRNAs on log₂ fold changes in the frequencies of different base editing outcomes in different screen conditions.',
    }
    
    scatter_layout = hits.visualize.interactive.scatter(df,
                                                        initial_data_lims=(-1.5, 1.5),
                                                        two_level_index=True,
                                                        return_layout=True,
                                                        identical_bins=True,
                                                        size=600,
                                                        menu_width=(200, 300),
                                                        color_by='color ',
                                                        grid='diagonal+axes',
                                                        initial_xy_names=(('HeLa BE4B rep1', 'Any C→G edit'), ('HeLa BE4B rep1', 'Any C→T edit')),
                                                        initial_selection=gl.gene_guides(['UNG', 'RFWD3', 'ASCC3']),
                                                        initial_alpha=0.7,
                                                        hide_widgets=['annotation', 'color_by', 'grid_radio_buttons'],
                                                        table_keys=['gene '],
                                                        level_order={
                                                           'Screen condition': [
                                                               'HeLa BE1 rep1',
                                                               'HeLa BE1 rep2',
                                                               'HeLa BE4B rep1',
                                                               'HeLa BE4B rep2',
                                                           ],
                                                            'Outcome category': [
                                                                'Any C→G edit',
                                                                'Any C→T edit',
                                                                'Exactly 1 C→T edit',
                                                                'Exactly 2 C→T edits',
                                                                'Exactly 3 C→T edits',
                                                                'Deletion',
                                                           ],
                                                       },
                                                      )

    templates.save_bokeh_html(scatter_layout,
                              '../BE_scatter_l2fcs.html',
                              description_data=description_data,
                              modal_data='scatter_tutorial.yaml',
                             )
    
    description_data = {
        'title': 'Base editing screens with 1,573 sgRNA CRISPRi library',
        'details': '''\
Scatter plots of the effect of each CRISPRi sgRNA on select phenotypes in BE1 and BE4B screens.
sgRNAs are ordered alphabetically along the x-axis to highlight when multiple sgRNAs targetting the same gene produce consistent phenotypes.''',
    }
    
    alphabetical_layout = alphabetical.scatter(data, plot_width=1000, plot_height=500, initial_genes=['RFWD3', 'UNG'], save_as='layout')

    templates.save_bokeh_html(alphabetical_layout,
                              '../BE_1573.html',
                              description_data=description_data,
                              modal_data='alphabetical_tutorial.yaml',
                             )    

    return scatter_layout, alphabetical_layout

def make_HDR():
    logging.info('Starting make_HDR()')
    
    base_dir = '/lab/solexa_weissman/jah/projects/ddr/dist'

    pools = rs.pooled_screen.get_all_pools(base_dir)

    outcomes_list_forward = [
        ('capture of donor fragment', [('donor misintegration', 'left unintended, right unintended')]),
        ('scarless HDR', [('donor', 'collapsed')]),
        ('half-HDR with 5\' homology unpaired', [('donor misintegration', 'left unintended, right intended')]),
        ('capture of human genomic sequence ≤75 nts', [('genomic insertion', 'hg19', '<=75 nts')]),
        ('capture of human genomic sequence >75 nts', [('genomic insertion', 'hg19', '>75 nts')]),
    ]

    outcomes_list_reverse = [
        ('capture of donor fragment', [('donor misintegration', 'left unintended, right unintended')]),
        ('scarless HDR', [('donor', 'collapsed')]),
        ('half-HDR with 5\' homology unpaired', [('donor misintegration', 'left intended, right unintended')]),
        ('capture of human genomic sequence ≤75 nts', [('genomic insertion', 'hg19', '<=75 nts')]),
        ('capture of human genomic sequence >75 nts', [('genomic insertion', 'hg19', '>75 nts')]),
    ]

    outcomes_list_NH = [
        ('capture of donor fragment', [('donor misintegration', 'left unintended, right unintended (no SNVs)')]),
        #('intended HDR', [('donor', 'collapsed')]),
        #('half-HDR with 5\' HA unpaired', [('donor misintegration', 'left unintended, right unintended (no SNVs)')]),
    ]
    
    logging.info('Making HDR_1573')

    pool_and_outcomes = [
        (pools['K562_SpCas9_target-1_oBA701_AX227_1'], outcomes_list_forward),
        (pools['K562_SpCas9_target-1_oBA701_AX227_2'], outcomes_list_forward),
        (pools['K562_SpCas9_target-1_oBA701-PCR_AX227_1'], outcomes_list_forward),
        (pools['K562_SpCas9_target-2_oBA701_AX227_1'], outcomes_list_forward),
    ]

    pool_and_outcomes = [(pool, (outcome_name, outcomes)) for pool, outcomes_list in pool_and_outcomes for outcome_name, outcomes in outcomes_list]
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        data = gene_significance.compute_table(pool_and_outcomes, 
                                               use_high_frequency_counts=True,
                                               pool_name_aliases={
                                                   'K562_SpCas9_target-1_oBA701_AX227_1': 'Cas9 target site 1, ssODN with uneven SNVs, + strand, rep 1',
                                                   'K562_SpCas9_target-1_oBA701_AX227_2': 'Cas9 target site 1, ssODN with uneven SNVs, + strand, rep 2',
                                                   'K562_SpCas9_target-1_oBA701-PCR_AX227_1': 'Cas9 target site 1, dsODN with uneven SNVs, + strand',
                                                   'K562_SpCas9_target-2_oBA701_AX227_1': 'Cas9 target site 2, ssODN with uneven SNVs, + strand, rep 1',
                                               },
                                               initial_dataset='Cas9 target site 1, ssODN with uneven SNVs, + strand, rep 1',
                                               initial_outcome='scarless HDR',
                                              )
    
    layout_1573 = alphabetical.scatter(data,
                                       outcome_names=[
                                           'scarless HDR',
                                           'half-HDR with 5\' homology unpaired',
                                           'capture of donor fragment',
                                           'capture of human genomic sequence ≤75 nts',
                                           'capture of human genomic sequence >75 nts',
                                       ],
                                       plot_width=1000,
                                       plot_height=500,
                                       initial_genes=['DNA2', 'BRCA2', 'PALB2', 'RAD51'],
                                       save_as='layout',
                                   )
    
    description_data = {
        'title': 'HDR screens with 1,573 sgRNA CRISPRi library',
        'details': '''\
Scatter plots of the effect of each CRISPRi sgRNA on select phenotypes in DSB screens performed with oligonucleotide HDR donors.
sgRNAs are ordered alphabetically along the x-axis to highlight when multiple sgRNAs targetting the same gene produce consistent phenotypes.''',
    }
    
    templates.save_bokeh_html(layout_1573,
                              '../HDR_1573.html',
                              description_data=description_data,
                              modal_data='alphabetical_tutorial.yaml',
                             )
    
    logging.info('Making HDR_366')
    
    pool_and_outcomes = [
        (pools['K562_SpCas9_target-1_oBA701_AC001_1'], outcomes_list_forward),
        (pools['K562_SpCas9_target-1_oBA701_AC001_2'], outcomes_list_forward),
        (pools['K562_SpCas9_target-1_oBA701-PAGE_AC001_1'], outcomes_list_forward),
        (pools['K562_SpCas9_target-1_oJAH158_AC001_1'], outcomes_list_reverse),
        (pools['K562_SpCas9_target-1_oJAH160_AC001_1'], outcomes_list_forward),
        (pools['K562_SpCas9_target-1_oJAH159_AC001_1'], outcomes_list_reverse),
        #(pools['K562_SpCas9_target-1_oJAH165_AC001_1'], outcomes_list_NH),
    ]

    pool_and_outcomes = [(pool, (outcome_name, outcomes)) for pool, outcomes_list in pool_and_outcomes for outcome_name, outcomes in outcomes_list]
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
    
        data = gene_significance.compute_table(pool_and_outcomes,
                                               use_high_frequency_counts=True,
                                               pool_name_aliases={
                                                   'K562_SpCas9_target-1_oBA701-PAGE_AC001_1': 'Cas9 target site 1, ssODN with uneven SNVs, + strand, PAGE-purified',
                                                   'K562_SpCas9_target-1_oBA701_AC001_1': 'Cas9 target site 1, ssODN with uneven SNVs, + strand, rep 1',
                                                   'K562_SpCas9_target-1_oBA701_AC001_2': 'Cas9 target site 1, ssODN with uneven SNVs, + strand, rep 2',
                                                   'K562_SpCas9_target-1_oJAH158_AC001_1': 'Cas9 target site 1, ssODN with uneven SNVs, - strand',
                                                   'K562_SpCas9_target-1_oJAH160_AC001_1': 'Cas9 target site 1, ssODN with even SNVs, + strand',
                                                   'K562_SpCas9_target-1_oJAH159_AC001_1': 'Cas9 target site 1, ssODN with even SNVs, - strand',
                                               },
                                              )
    
    description_data = {
        'title': 'HDR screens with 366 sgRNA CRISPRi library',
        'details': '''\
Scatter plots of the effect of each CRISPRi sgRNA on select phenotypes in DSB screens performed with oligonucleotide HDR donors.
sgRNAs are ordered alphabetically along the x-axis to highlight when multiple sgRNAs targetting the same gene produce consistent phenotypes.''',
    }
    
    layout_366 = alphabetical.scatter(data,
                                      outcome_names=[
                                          'scarless HDR',
                                          'half-HDR with 5\' homology unpaired',
                                          'capture of donor fragment',
                                          'capture of human genomic sequence ≤75 nts',
                                          'capture of human genomic sequence >75 nts',
                                      ],
                                      plot_width=900,
                                      plot_height=500,
                                      initial_genes=['DNA2', 'LIG4'],
                                      save_as='layout',
                                  )
    
    templates.save_bokeh_html(layout_366,
                              '../HDR_366.html',
                              description_data=description_data,
                              modal_data='alphabetical_tutorial.yaml',
                             )
    
    return layout_1573, layout_366
    
def make_UMAPs():
    logging.info('Starting make_UMAPs()')
    
    base_dir = '/lab/solexa_weissman/jah/projects/ddr/dist'

    pools = rs.pooled_screen.get_all_pools(base_dir)
    
    pns = [
        'K562_SpCas9_target-1_none_AC001_1',
        'K562_SpCas9_target-1_none_AC001_2',
        'K562_SpCas9_target-4_none_AC001_1',
        'K562_SpCas9_target-4_none_AC001_2',
        'K562_SpCas9_target-3_none_AC001_1',
        'K562_SpCas9_target-3_none_AC001_2',
        'K562_SpCas9_target-2_none_AC001_1',
    ]

    clusterer = rs.cluster.MultiplePoolClusterer([pools[pn] for pn in pns],
                                                 outcomes_method='HDBSCAN',
                                                 outcomes_selection_method='above_frequency_threshold',
                                                 outcomes_selection_kwargs=dict(threshold=2e-3),
                                                 outcomes_min_dist=0.5,
                                                 use_high_frequency_counts=True,
                                                )

    with warnings.catch_warnings():
        # Suppress a warning about TBB_INTERFACE_VERSION
        warnings.simplefilter('ignore')

        # Compute memoized outcome_embedding in warning context
        clusterer.outcome_embedding
        
    layout = UMAP.plot(clusterer)
    
    description_data = {
        'title': 'UMAP embedding of genetic dependencies of repair outcomes at Cas9-induced DSBs',
        'details': '''\
The upper left panel summarizes the genetic dependencies of individual repair outcomes from screens at four different Cas9 target sites, as in <a href="https://www.sciencedirect.com/science/article/pii/S0092867421011764#fig3">Figure 3G here</a>.
Outcomes are initially colored by the effect of a POLQ-targeting CRISPRi sgRNA on their frequencies.
Use the menu to the right to select different CRISRPi sgRNAs, or the tabs above to color outcomes by other properties.
The bottom panel shows the correlation in outcome frequency redistribution between the currently selected CRISPRi sgRNA and all other sgRNAs.
Mouse over points in the bottom panel to compare the effects of other sgRNAs.
''',
    }

    templates.save_bokeh_html(layout,
                              '../Cas9_UMAP.html',
                              description_data=description_data,
                              modal_data='UMAP_tutorial.yaml',
                             )
    
    return layout
    
def make_PE_scatters():
    logging.info('Starting make_PE_scatters()')
    
    base_dir = '/lab/solexa_weissman/jah/projects/prime_editing_screens/dist'

    pools = rs.pooled_screen.get_all_pools(base_dir)
    
    gl = rs.guide_library.GuideLibrary(base_dir, 'AX227')
    
    pool_names = [
        'K562_PE2_rep1',
        'K562_PE2_rep2',
        'K562_PE3+50_rep1',
        'K562_PE3+50_rep2',
        'HeLa_PE2_rep1',
        'HeLa_PE2_rep2',
        'HeLa_PE3+50_rep1',
        'HeLa_PE3+50_rep2',
        'HeLa_PE3-50_rep1',
        'HeLa_PE3-50_rep2',
    ]

    categories = [
        'intended edit',
        'unintended annealing of RT\'ed sequence',
        'extension from intended annealing',
        'deletion',
        'duplication',
    ]

    cat_aliases = {
        'intended edit': 'Intended edit',
        'deletion': 'Deletion',
        'duplication': 'Tandem duplication',
        'unintended annealing of RT\'ed sequence': 'Joining of RT\'ed sequence at unintended location',
        'extension from intended annealing': 'Installation of additional edits from nearly-matched scaffold sequence',
    }

    edit_percentages = {}

    for pn in pool_names:
        for cat in categories:
            pool = pools[pn]
            edit_percentages[pn, cat] = pool.category_fractions.loc[cat] * 100

    edit_percentages = pd.DataFrame(edit_percentages)

    edit_percentages.index.name = 'CRISPRi sgRNA'
    edit_percentages.columns.names = ['Screen condition', 'Outcome category']
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        
        fcs = edit_percentages / edit_percentages.loc['all_non_targeting']
        l2fcs = np.log2(fcs)
        
    colors = pd.Series('grey', index=edit_percentages.index)

    gene_to_color = {
        'MSH2': 'tab:orange',
        'MSH6': 'tab:orange',
        'PMS2': 'tab:green',
        'MLH1': 'tab:green',
        'FEN1': 'tab:purple',
        'LIG1': 'tab:purple',
        'HLTF': 'tab:red',
        'EXO1': 'tab:blue',
    }

    for nt_guide_set in gl.non_targeting_guide_sets:
        gene_to_color[nt_guide_set] = 'black'

    for gene, color in gene_to_color.items():
        colors[gl.gene_guides(gene)] = color

    colors[rs.pooled_screen.ALL_NON_TARGETING] = 'black'

    colors = colors.map(matplotlib.colors.to_hex)
    
    for df in [edit_percentages, l2fcs]:
        df['color'] = colors

        df.drop(rs.pooled_screen.ALL_NON_TARGETING, inplace=True)

        df['gene'] = df.index.map(gl.guide_to_gene_with_non_targeting_guide_sets)

        df.rename(columns=lambda s: s.replace('_', ' '), level=0, inplace=True)
        df.rename(columns=cat_aliases, level=1, inplace=True)
        
    condition_order = [
        'K562 PE2 rep1',
        'K562 PE2 rep2',
        'K562 PE3+50 rep1',
        'K562 PE3+50 rep2',
        'HeLa PE2 rep1',
        'HeLa PE2 rep2',
        'HeLa PE3+50 rep1',
        'HeLa PE3+50 rep2',
        'HeLa PE3-50 rep1',
        'HeLa PE3-50 rep2',
    ]

    category_order = [cat_aliases.get(cat, cat) for cat in categories]

    common_kwargs = dict(
        table_keys=['gene '],
        size=600,
        identical_bins=True,
        two_level_index=True,
        color_by='color ',
        hide_widgets=['annotation', 'color_by', 'grid_radio_buttons'],
        initial_alpha=0.7,
        menu_width=(200, 300),
        level_order={
            'Outcome category': category_order,
            'Screen condition': condition_order,
        },
        return_layout=True,
    )

    layout_percentages = hits.visualize.interactive.scatter(edit_percentages,
                                                initial_xy_names=(('K562 PE2 rep1', 'Intended edit'), ('K562 PE2 rep2', 'Intended edit')),
                                                data_bounds=(0, 100),
                                                grid='grid',
                                                **common_kwargs,
                                               )

    description_data = {
        'title': 'Scatter plots of prime editing outcome frequencies',
        'details': 'Comparisons of the effects of CRISPRi sgRNAs on the frequencies of different prime editing outcomes in different screen conditions.',
    }

    templates.save_bokeh_html(layout_percentages,
                              '../PE_scatter_frequencies.html',
                              description_data=description_data,
                              modal_data='scatter_tutorial.yaml',
                             )
    
    layout_l2fcs = hits.visualize.interactive.scatter(l2fcs.replace([np.inf, -np.inf], np.nan),
                                                      initial_xy_names=(('K562 PE2 rep1', 'Intended edit'), ('K562 PE3+50 rep1', 'Intended edit')),
                                                      initial_data_lims=(-2.25, 3),
                                                      data_bounds=(-5, 5),
                                                      grid='diagonal+axes',
                                                      **common_kwargs,
                                                     )

    description_data = {
        'title': 'Scatter plots of log₂ fold changes in prime editing outcome frequencies',
        'details': 'Comparisons of the effects of CRISPRi sgRNAs on log₂ fold changes in the frequencies of different prime editing outcomes in different screen conditions.',
    }

    templates.save_bokeh_html(layout_l2fcs,
                              '../PE_scatter_l2fcs.html',
                              description_data=description_data,
                              modal_data='scatter_tutorial.yaml',
                             )
    
    return layout_percentages, layout_l2fcs

def make_PE_PCs():
    logging.info('Starting make_PE_PCs()')
    
    base_dir = '/lab/solexa_weissman/jah/projects/prime_editing_screens/dist'

    pools = rs.pooled_screen.get_all_pools(base_dir)
    
    gl = rs.guide_library.GuideLibrary(base_dir, 'AX227')
    
    pool_names = [
        'K562_PE2_rep1',
        'K562_PE2_rep2',
        'K562_PE3+50_rep1',
        'K562_PE3+50_rep2',
        'HeLa_PE2_rep1',
        'HeLa_PE2_rep2',
        'HeLa_PE3+50_rep1',
        'HeLa_PE3+50_rep2',
        'HeLa_PE3-50_rep1',
        'HeLa_PE3-50_rep2',
    ]

    edit_percentages = {}

    for pn in pool_names:
        pool = pools[pn]
        edit_percentages[pn] = pool.category_fractions.loc['intended edit'] * 100

    edit_percentages = pd.DataFrame(edit_percentages)
    edit_percentages.rename(columns=lambda s: s.replace('_', ' '), inplace=True)
    
    fcs = edit_percentages / edit_percentages.loc['all_non_targeting']

    l2fcs = np.log2(fcs)
    
    gene_to_color = {
        'MSH2': 'tab:orange',
        'MSH6': 'tab:orange',
        'PMS2': 'tab:green',
        'MLH1': 'tab:green',
        'FEN1': 'tab:purple',
        'LIG1': 'tab:purple',
        'HLTF': 'tab:red',
        'EXO1': 'tab:blue',
    }
    
    for nt_guide_set in gl.non_targeting_guide_sets:
            gene_to_color[nt_guide_set] = 'black'

    for df in [l2fcs, edit_percentages]:
        df['color'] = 'grey'
        for gene, color in gene_to_color.items():
            df.loc[gl.gene_guides(gene), 'color'] = color
            
        df.loc[rs.pooled_screen.ALL_NON_TARGETING, 'color'] = 'black'
        
        df.index.name = 'CRISPRi sgRNA'
        
    description_data = {
        'title': 'Parallel coordinates plot of intended prime edit frequencies',
        'details': 'Lines represent individual CRISPRi sgRNAs, and columns represent distinct screens. Lines plots the frequency of the intended substitution prime edit in the presence of each sgRNA across 10 different screen conditions.',
    }

    templates.save_PC_html(edit_percentages,
                           '../PE_PC_frequencies.html',
                           description_data=description_data,
                           modal_data='PC_tutorial.yaml',
                           initial_limits=(0, 50),
                           y_axis_label='Percentage of sequencing reads with intended edit',
                           min_start=0, max_end=50,
                           label_by='CRISPRi sgRNA',
                          )
    
    description_data = {
        'title': 'Parallel coordinates plot of log₂ fold changes in intended prime edit frequencies',
        'details': 'Lines represent individual CRISPRi sgRNAs, and columns represent distinct screens. Lines plots the log₂ fold change in frequency of the intended substitution prime edit in the presence of each sgRNA across 10 different screen conditions.',
    }

    templates.save_PC_html(l2fcs,
                           '../PE_PC_l2fcs.html',
                           description_data=description_data,
                           modal_data='PC_tutorial.yaml',
                           initial_limits=(-2, 3),
                           y_axis_label='Log₂ fold change from non-targeting in intended edit',
                           min_start=-2, max_end=3,
                           label_by='CRISPRi sgRNA',
                          )
    

if __name__ == '__main__':
    make_index()
    #make_UMAPs()
    make_HDR()
    make_BE_scatters()
    make_PE_scatters()
    make_PE_PCs()
