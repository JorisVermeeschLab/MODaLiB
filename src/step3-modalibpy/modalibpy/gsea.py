from typing import List, Optional

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import gseapy
import kaleido


class GSEA(object):

    def __init__(self, gene_names: List[str], background: Optional[List[str]] = None):
        self.gene_names: List[str] = [gene_name for gene_name in gene_names if isinstance(gene_name, str)]
        self.background: Optional[List[str]] = [gene_name for gene_name in background if isinstance(gene_name, str)]
        self.enr = gseapy.enrichr(
            gene_list=self.gene_names,
            gene_sets=['KEGG_2021_Human'],
            #gene_sets=['Cancer_Cell_Line_Encyclopedia'],
            background=self.background,
            organism='human',
            outdir=None
        )

    def plot(self) -> go.Figure:
        df = self.enr.results[['Term', 'P-value', 'Adjusted P-value', 'Genes']][:10]
        df['P-value'] = df['P-value'].apply(lambda x: f'{x:.3f}')
        df['Adjusted P-value'] = df['Adjusted P-value'].apply(lambda x: f'{x:.3f}')

        fig = go.Figure(data=[go.Table(
          columnorder = [1, 2, 3, 4],
          columnwidth = [150, 50, 50, 200],
          header = dict(
            values = [
                ['<b>Term</b>'],
                ['<b>p-value</b>'],
                ['<b>Adjusted p-value</b>'],
                ['<b>Genes</b>'],
            ],
            line_color='darkslategray',
            fill_color='royalblue',
            align=['left','center'],
            font=dict(color='white', size=12),
            height=40
          ),
          cells=dict(
            values=df.to_numpy().T,
            line_color='darkslategray',
            fill=dict(color=['paleturquoise', 'white']),
            align=['left', 'center'],
            font_size=12,
            height=30)
            )
        ])
        fig.update_layout(autosize=False, width=1200, height=600)
        plt.tight_layout()
        return fig
