# Chromosome sizes from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
CHROM_LENGTHS = {
        'chr1': 248_956_422,
        'chr2': 242_193_529,
        'chr3': 198_295_559,
        'chr4': 190_214_555,
        'chr5': 181_538_259,
        'chr6': 170_805_979,
        'chr7': 159_345_973,
        'chrX': 156_040_895,
        'chr8': 145_138_636,
        'chr9': 138_394_717,
        'chr11': 135_086_622,
        'chr10': 133_797_422,
        'chr12': 133_275_309,
        'chr13': 114_364_328,
        'chr14': 107_043_718,
        'chr15': 101_991_189,
        'chr16': 90_338_345,
        'chr17': 83_257_441,
        'chr18': 80_373_285,
        'chr20': 64_444_167,
        'chr19': 58_617_616,
        'chrY': 57_227_415,
        'chr22': 50_818_468,
        'chr21': 46_709_983
    }


# PhyloP scores from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP30way/hg38.30way.phyloP/
PHYLOP_BASE_URL = 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP30way/hg38.30way.phyloP'
PHYLOP_WIGS = [
    'chr1.phyloP30way.wigFix',
    'chr2.phyloP30way.wigFix',
    'chr3.phyloP30way.wigFix',
    'chr4.phyloP30way.wigFix',
    'chr5.phyloP30way.wigFix',
    'chr6.phyloP30way.wigFix',
    'chr7.phyloP30way.wigFix',
    'chr8.phyloP30way.wigFix',
    'chr9.phyloP30way.wigFix',
    'chr10.phyloP30way.wigFix',
    'chr11.phyloP30way.wigFix',
    'chr12.phyloP30way.wigFix',
    'chr13.phyloP30way.wigFix',
    'chr14.phyloP30way.wigFix',
    'chr15.phyloP30way.wigFix',
    'chr16.phyloP30way.wigFix',
    'chr17.phyloP30way.wigFix',
    'chr18.phyloP30way.wigFix',
    'chr19.phyloP30way.wigFix',
    'chr20.phyloP30way.wigFix',
    'chr21.phyloP30way.wigFix',
    'chr22.phyloP30way.wigFix',
    'chrX.phyloP30way.wigFix',
    'chrY.phyloP30way.wigFix'
]

SNP_REGION_SIDE_LENGTHS = [10_000, 5_000, 2_500, 1_000, 100, 10]


LABEL_COLOR_MAP = {
    'Promoter': (1.0, 0.0, 0.0),
    'PromoterFlanking': (1.0, 0.26666666666666666, 0.0),
    'Enhancer': (1.0, 0.7647058823529411, 0.30196078431372547),
    'EnhancerLow': (1.0, 1.0, 0.0),
    'Bivalent': (0.7411764705882353, 0.7176470588235294, 0.4196078431372549),
    'CTCF': (0.7686274509803922, 0.8823529411764706, 0.0196078431372549),
    'Transcribed': (0.0, 0.5019607843137255, 0.0),
    'K9K36': (0.4, 0.803921568627451, 0.6666666666666666),
    'FacultativeHet': (0.5019607843137255, 0.0, 0.5019607843137255),
    'ConstitutiveHet': (0.5411764705882353, 0.5686274509803921, 0.8156862745098039),
    'Quiescent': (1.0, 1.0, 1.0)
}