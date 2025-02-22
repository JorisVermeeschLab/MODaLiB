import enum

from .cell_type_deconvolution import FragCellTypeDeconvolutionModality, MethCellTypeDeconvolutionModality
from .coverage import CoverageModality, CAFFModality
from .cpg_methylation import CpGMethylationModality, FullyMethylatedFragModality
from .end_motif_frequency import EndMotifFrequencyModality
from .frag_length_counts import FragLengthCountModality
from .jaggedness import JaggednessModality
from .long_frag_ratio import LongFragmentRatioModality
from .mutational import MutationalSignatureModality
from .nucleosome_positioning import NucleosomePositioningModality
from .context import TrinucleotideContext


class Modalities(enum.Enum):

    NUCLEOSOME_POSITIONING = enum.auto()
    LONG_FRAG_RATIO = enum.auto()
    FRAG_CELL_TYPE_DECONVOLUTION = enum.auto()
    METH_CELL_TYPE_DECONVOLUTION = enum.auto()
    FRAG_LENGTH_COUNTS = enum.auto()
    CPG_METHYLATION = enum.auto()
    FULLY_METHYLATED_FRAG = enum.auto()
    FULLY_METH_FRAG = enum.auto()
    COVERAGE = enum.auto()
    CAFF = enum.auto()
    END_MOTIFS = enum.auto()
    JAGGEDNESS = enum.auto()
    MUTATIONS = enum.auto()
    TRINUCLEOTIDE_CONTEXT = enum.auto()

    @classmethod
    def list(cls):
        return list(map(lambda c: c.value, cls))


MODALITY_CLASSES = {
    Modalities.NUCLEOSOME_POSITIONING: NucleosomePositioningModality,
    Modalities.LONG_FRAG_RATIO: LongFragmentRatioModality,
    Modalities.FRAG_CELL_TYPE_DECONVOLUTION: FragCellTypeDeconvolutionModality,
    Modalities.METH_CELL_TYPE_DECONVOLUTION: MethCellTypeDeconvolutionModality,
    Modalities.FRAG_LENGTH_COUNTS: FragLengthCountModality,
    Modalities.CPG_METHYLATION: CpGMethylationModality,
    Modalities.FULLY_METH_FRAG: FullyMethylatedFragModality,
    Modalities.COVERAGE: CoverageModality,
    Modalities.CAFF: CAFFModality,
    Modalities.END_MOTIFS: EndMotifFrequencyModality,
    Modalities.JAGGEDNESS: JaggednessModality,
    Modalities.MUTATIONS: MutationalSignatureModality,
    Modalities.TRINUCLEOTIDE_CONTEXT: TrinucleotideContext
}
