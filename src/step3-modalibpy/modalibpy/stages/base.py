import re
import math
from typing import Union, Optional

import numpy as np


class Grade:

    def __init__(self, s: Union[str, float, int]) -> None:

        self.s: str = str(s)

        if (not isinstance(s, str)) and np.isnan(s):
            s = 'X'
        if s == '?':
            s = 'X'

        if isinstance(s, str):
            s = s.replace(' ', '').upper()
            self.g1: float = np.nan
            self.g2: float = np.nan
            res = re.findall(re.compile('(X|0|1|2|3|4|IV|III|II|I)([ABC]?)'), s)
            if len(res) == 0:
                raise ValueError(f'Unknown cancer grade "{s}"')
            if res[0][0] == 'X':
                self.g1 = np.nan
            elif res[0][0] == '0':
                self.g1 = 0.
            elif res[0][0] in {'I', '1'}:
                self.g1 = 1.
            elif res[0][0] in {'II', '2'}:
                self.g1 = 2.
            elif res[0][0] in {'III', '3'}:
                self.g1 = 3.
            elif res[0][0] in {'IV', '4'}:
                self.g1 = 4.
            else:
                raise ValueError(f'Unknown cancer grade "{s}"')
            if res[0][1] == 'a':
                self.g2 = 0.1
            elif res[0][1] == 'b':
                self.g2 = 0.2
            elif res[0][1] == 'c':
                self.g2 = 0.3
            else:
                self.g2 = 0.
        else:
            g1 = int(round(s - s % 1.0))
            if g1 == 0:
                self.g1 = 0.
            elif g1 == 1:
                self.g1 = 1.
            elif g1 == 2:
                self.g1 = 2.
            elif g1 == 3:
                self.g1 = 3.
            elif g1 == 4:
                self.g1 = 4.
            else:
                raise ValueError(f'Unknown cancer grade "{s}"')
            self.g2 = 0.
        self.g: float = self.g1 + self.g2

    def unknown(self) -> bool:
        return bool(np.isnan(self.g1))

    @staticmethod
    def ensure_grade(o: object) -> 'Grade':
        if isinstance(o, Grade):
            return o
        else:
            return Grade(o)

    def __float__(self) -> float:
        return self.g

    def __eq__(self, other: 'Grade') -> bool:
        return self.g == Grade.ensure_grade(other).g

    def __lt__(self, other: 'Grade') -> bool:
        return self.g < Grade.ensure_grade(other).g

    def __ge__(self, other: 'Grade') -> bool:
        return self.g >= Grade.ensure_grade(other).g

    def __str__(self) -> str:
        return self.s

    def __repr__(self) -> str:
        return self.__str__()


class Annotation:

    TN_STAGING_SYSTEM_REGEX = r'[cp]T([0-9X][abc]?)N([0-3X][ab]?)'
    TNM_STAGING_SYSTEM_REGEX = r'[cp]T([0-9X][abc]?)N([0-3X][ab]?)M([0-1X])'

    def __init__(self, is_cancer: bool=False, grade: str='', extra_info: str='', cancer_type: Optional[str] = None) -> None:

        self.is_cancer: bool = is_cancer
        self.cancer_type: Optional[str] = cancer_type
        self.extra_info: str = extra_info
        self.msi: float = np.nan
        self.subtype_: str = ""

        # TNM system
        self.tumor_size: Grade = Grade('X')
        self.lymph_node_involvement: Grade = Grade('X')
        self.metastasis: Grade = Grade('X')

        if self.is_cancer:

            self.grade: Grade = Grade(grade)

            # If multiple annotations, consider only the most advanced cancer stage
            for res in re.findall(Annotation.TNM_STAGING_SYSTEM_REGEX, extra_info):
                if self.tumor_size.unknown() or (Grade(res[0]) > self.tumor_size):
                    self.tumor_size = Grade(res[0])
                if self.lymph_node_involvement.unknown() or (Grade(res[1]) > self.lymph_node_involvement):
                    self.lymph_node_involvement = Grade(res[1])
                if self.metastasis.unknown() or (Grade(res[2]) > self.metastasis):
                    self.metastasis = Grade(res[2])
            for res in re.findall(Annotation.TN_STAGING_SYSTEM_REGEX, extra_info):
                if self.tumor_size.unknown() or (Grade(res[0]) > self.tumor_size):
                    self.tumor_size = Grade(res[0])
                if self.lymph_node_involvement.unknown() or (Grade(res[1]) > self.lymph_node_involvement):
                    self.lymph_node_involvement = Grade(res[1])

            # Check metastases
            if (self.grade >= 4) and ('metastasis' in extra_info):
                self.metastasis = Grade('1')

        else:
            self.grade: Grade = Grade('X')

    def metastasized(self) -> float:
        if self.metastasis.unknown():
            return np.nan
        else:
            return float(self.metastasis)

    def __str__(self) -> str:
        if self.is_cancer:
            annotations = [
                str(self.grade),
                f'T{self.tumor_size}N{self.lymph_node_involvement}M{self.metastasis}'
            ]
            return f'CancerStage({", ".join(annotations)})'
        else:
            return f'NotACancer()'

    def __repr__(self) -> str:
        return self.__str__()

    @staticmethod
    def ensure_stage(o: object):
        if isinstance(o, Annotation):
            return o
        elif isinstance(o, int):
            return Annotation(True, str(o), '')
        elif isinstance(o, float):
            return Annotation(True, str(o), '')
        elif isinstance(o, str):
            return Annotation(True, o, '')
        else:
            raise ValueError(f'Invalid cancer stage"{o}"')

    def __float__(self) -> float:
        return float(self.grade)

    def __eq__(self, other: 'Annotation') -> bool:
        return float(self.grade) == float(Annotation.ensure_stage(other).grade)

    def __le__(self, other: 'Annotation') -> bool:
        return float(self.grade) <= float(Annotation.ensure_stage(other).grade)

    def __lt__(self, other: 'Annotation') -> bool:
        return float(self.grade) < float(Annotation.ensure_stage(other).grade)

    def __ge__(self, other: 'Annotation') -> bool:
        return float(self.grade) >= float(Annotation.ensure_stage(other).grade)

    def __gt__(self, other: 'Annotation') -> bool:
        return float(self.grade) > float(Annotation.ensure_stage(other).grade)

    def subtype(self) -> str:
        return self.subtype_
    
