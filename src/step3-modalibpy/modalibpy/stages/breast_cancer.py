import numpy as np

from modalibpy.stages.base import Annotation, Grade


class BreastCancerStage(Annotation):

    def __init__(self, *args, **kwargs) -> None:
        Annotation.__init__(self, *args, cancer_type='breast', **kwargs)

        sf = self.extra_info.strip().lower()

        if sf == 'her2':
            er, pr, her2 = 0., 0., 1.
        elif sf == 'triple negative':
            er, pr, her2 = 0., 0., 0.
        elif sf in {'luminal a like', 'luminal a'}:
            er, pr, her2 = 1., 1., 0.
        elif sf in {'luminal b like', 'luminal b'}:
            er, pr, her2 = 1., 1., 1.
        elif sf == '':
            er, pr, her2 = np.nan, np.nan, np.nan
        else:

            er = np.nan
            if 'er+' in sf:
                er = 1.
            if 'er8/8' in sf:
                er = 1.
            if 'er-' in sf:
                er = 0.
            
            pr = np.nan
            if 'pr+' in sf:
                pr = 1.
            if 'pr8/8' in sf:
                pr = 1.
            if 'pr-' in sf:
                pr = 0.

            her2 = np.nan
            if 'her2+' in sf:
                her2 = 1.
            elif 'HER21+' in sf:
                her2 = 0.
            elif 'her2-' in sf:
                her2 = 0.
            elif 'her2' in sf:
                her2 = 1.
            elif 'her2negatief' in sf:
                her2 = 0.

            if 'luminalb' in sf:
                er, pr, her2 = 1., 1., 1.
            elif 'luminala' in sf:
                er, pr, her2 = 1., 1., 0.
            elif 'triplepositief' in sf:
                er, pr, her2 = 1., 1., 1.
            elif 'triplenegative' in sf:
                er, pr, her2 = 0., 0., 0.

        self.er = er
        self.pr = pr
        self.her2 = her2

        # Some more hardcoding
        if 'grII' in self.extra_info:
            self.grade = Grade('II')

    def triple_negative(self) -> float:
        if (self.er == 0) and (self.pr == 0) and (self.her2 == 0):
            return 1.
        elif np.isnan(self.er + self.pr + self.her2):
            return np.nan
        else:
            return 0.

    def luminal_a(self) -> float:
        return (self.er or self.pr) and (not self.her2)

    def luminal_b(self) -> float:
        return (self.er or self.pr) and (self.her2)

    def er_positive(self) -> float:
        return self.er

    def pr_positive(self) -> float:
        return self.pr

    def her2_positive(self) -> float:
        return self.her2

    def subtype(self) -> str:
        if self.triple_negative() == 1.:
            return 'Triple negative'
        elif self.luminal_a() == 1.:
            return 'Luminal A'
        elif self.luminal_b() == 1.:
            return 'Luminal B'
        elif self.her2_positive() == 1.:
            return 'HER2-positive'
        elif float(self) == 0.:
            return 'DCIS'
        else:
            return 'Unknown'

    def __str__(self) -> str:
        if self.is_cancer:
            annotations = []
            annotations.append(str(self.grade))
            annotations.append(f'T{self.tumor_size}N{self.lymph_node_involvement}M{self.metastasis}')
            for attr, name in zip([self.er, self.pr, self.her2], ['ER', 'PR', 'HER2']):
                if not np.isnan(attr):
                    if attr:
                        annotations.append(f'{name}+')
                    else:
                        annotations.append(f'{name}-')
            return f'BreastCancerStage({", ".join(annotations)})'
        else:
            return f'NotACancer()'
