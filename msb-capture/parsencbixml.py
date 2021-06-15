#!/usr/bin/python
"""
Parses an NCBI INSDSeq XML (-f) file and creates a fasta file (-o) with
either the full source sequences (-t s); gene sequences (-t g); coding sequnces
(CDS) (-t c); mature peptide sequences (-t m); 5'-UTR sequences (-t 5) or
3'-UTR sequences (-t 3). For genes, CDS and mature peptides a name (-n) is
required that also can be a list of names separated by blank spaces. If the
names contain blank spaces enclose in quotation marks. Besides the name (s),
which are required an alias file (-a) can be prepared with ALIASES for each
feature on single line separated by commas. The name given by (-n) must be
included among the ALIASES. The alias file can be combined with with a list of
names given with the name option and only one of the names need to to be in the
alias file
"""


class SeqFeatures:
    """General information about the serquence entries"""

    def __init__(self, seq):
        self.seq_locus = seq.find('INSDSeq_locus').text
        self.seq_length = seq.find('INSDSeq_length').text
        self.seq_strandedness = seq.find('INSDSeq_strandedness').text
        self.seq_moltype = seq.find('INSDSeq_moltype').text
        self.seq_topology = seq.find('INSDSeq_topology').text
        self.seq_division = seq.find('INSDSeq_division').text
        self.seq_update = seq.find('INSDSeq_update-date').text
        self.seq_create_date = seq.find('INSDSeq_create-date').text
        self.seq_definition = seq.find('INSDSeq_definition').text
        self.seq_primary_accession = seq.find('INSDSeq_primary-accession').text
        self.seq_accession_ver = seq.find('INSDSeq_accession-version').text
        self.seq_other_ids = list()
        for item in seq.findall('INSDSeq_other-seqids/INSDSeqid'):
            self.seq_other_ids.append(item.text)
        self.seq_keywords = list()
        if not seq.find('INSDSeq_project') is None:
            self.seq_project = seq.find('INSDSeq_project').text
        for item in seq.findall('INSDSeq_keywords/INSDKeyword'):
            self.seq_keywords.append(item.text)
        self.seq_source = seq.find('INSDSeq_source').text
        self.seq_organism = seq.find('INSDSeq_organism').text
        self.seq_taxonomy = seq.find('INSDSeq_taxonomy').text
        self.seq_references = list()
        reference = dict()
        for item in seq.findall('INSDSeq_references/INSDReference'):
            reference['reference'] = item.find('INSDReference_reference').text
            reference['reference_position'] = \
                item.find('INSDReference_position').text
            if not item.find('INSDReference_consortium') is None:
                reference['reference_consortium'] = \
                    item.find('INSDReference_consortium').text
            reference['reference_title'] = item.find('INSDReference_title').text
            reference['reference_journal'] = \
                item.find('INSDReference_journal').text
            reference['reference_title'] = item.find('INSDReference_title').text
            if not item.find('INSDReference_authors') is None:
                authors = list()
                for author in item.findall('INSDReference_authors/INSDAuthor'):
                    authors.append(author.text)
                reference['reference_authors'] = authors
            self.seq_references.append(reference)
            reference = dict()
        if not seq.find('INSDSeq_comment') is None:
            self.seq_comment = seq.find('INSDSeq_comment').text
        if not seq.find('INSDSeq_sequence') is None:
            self.seq_sequence = seq.find('INSDSeq_sequence').text


class FivepUTR:
    """Class defining the genetic region and sequence of an 5'UTR"""
    def __init__(self, ft):
        self.feature_key = ft.find('INSDFeature_key').text
        self.feature_location = ft.find('INSDFeature_location').text
        interval = dict()
        intervallst = list()
        for item in ft.findall('INSDFeature_intervals/INSDInterval'):
            if not item.find('INSDInterval_from') is None:
                interval['interval_from'] = item.find('INSDInterval_from').text
                interval['interval_to'] = item.find('INSDInterval_to').text
                interval['interval_acc'] = item.find(
                    'INSDInterval_accession').text
            elif not item.find('INSDInterval_point') is None:
                interval['interval_from'] = item.find('INSDInterval_point').text
                interval['interval_to'] = item.find('INSDInterval_point').text
                interval['interval_acc'] = item.find(
                    'INSDInterval_accession').text
            else:
                interval['interval_from'] = ""
                interval['interval_to'] = ""
                interval['interval_acc'] = ""
            intervallst.append(interval)
            interval = dict()
        self.feature_intervals = intervallst


class Source(FivepUTR):
    """Class defining the source sequence"""
    def __init__(self, ft):
        super().__init__(ft)
        qualifiers = dict()
        for qualifier in ft.findall('INSDFeature_quals/INSDQualifier'):
            if not qualifier.find('INSDQualifier_name') is None:
                name = qualifier.find('INSDQualifier_name').text
            else:
                continue
            if not qualifier.find('INSDQualifier_value') is None:
                value = qualifier.find('INSDQualifier_value').text
            else:
                value = ""
            qualifiers[name] = value
        self.qualifiers = qualifiers


class CDS(Source):
    """Class defining the coding region"""
    def __init__(self, ft):
        super().__init__(ft)
        if not ft.find('INSDFeature_operator') is None:
            self.feature_operator = ft.find('INSDFeature_operator').text

#
# class MiscFeature(Source):
#     """TBD"""
#     def __init__(self, ft):
#         super().__init__(ft)
#
#
# class MiscRNA(Source):
#     """TBD"""
#     def __init__(self, ft):
#         super().__init__(ft)
#
#
# class StemLoop(FivepUTR):
#     """TBD"""
#     def __init__(self, ft):
#         super().__init__(ft)
#
#
# class Regulatory(Source):
#     """TBD"""
#     def __init__(self, ft):
#         super().__init__(ft)


def extract_seq():
    """Extract the sequences and create sequence fasta-type id string and write
    to a fasta file (-o)"""
    seq_count = 0
    seq_found = 0
    outfile = open(ARGS.o, 'w')
    for sequence in ROOT.iter('INSDSeq'):
        seq_count += 1
        seq_name = ""
        seq = ""
        int_from = set()
        int_to = set()
        seqfeatures = SeqFeatures(sequence)
        for feature in sequence.findall(
                'INSDSeq_feature-table/INSDFeature'):
            acc = seqfeatures.seq_primary_accession
            organism = seqfeatures.seq_organism.replace(' ', '_')
            if ARGS.t == 's' and feature.find('INSDFeature_key').text == \
                    'source':
                source = Source(feature)
                seq_found += 1
                for interval in source.feature_intervals:
                    int_from.add(int(interval['interval_from']))
                    int_to.add(int(interval['interval_to']))
                min_int = min(int_from)
                max_int = max(int_to)
                seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + '\n'
            if ARGS.t == 'g' and feature.find('INSDFeature_key').text == 'gene':
                gene = Source(feature)
                if 'gene' not in gene.qualifiers.keys() or \
                        gene.qualifiers['gene'].casefold() not in FEATURE_NAMES:
                    continue
                gene_name = gene.qualifiers['gene'].replace(' ', '_')
                for interval in gene.feature_intervals:
                    int_from.add(int(interval['interval_from']))
                    int_to.add(int(interval['interval_to']))
                min_int = min(int_from)
                max_int = max(int_to)
                seq_found += 1
                seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';gene=' + gene_name +\
                           '\n'
                break
            if ARGS.t == 'c' and feature.find('INSDFeature_key').text == 'CDS':
                cds = CDS(feature)
                if 'product' not in cds.qualifiers.keys() or\
                        cds.qualifiers['product'].casefold()\
                        not in FEATURE_NAMES:
                    continue
                cds_name = cds.qualifiers['product'].replace(' ', '_')
                if ARGS.s == 'a':
                    seq = cds.qualifiers['translation']
                else:
                    for interval in cds.feature_intervals:
                        int_from.add(int(interval['interval_from']))
                        int_to.add(int(interval['interval_to']))
                    min_int = min(int_from)
                    max_int = max(int_to)
                    seq_found += 1
                    seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                    seq_name = '>' + acc + ';' + organism + ';CDS=' +\
                               cds_name + '\n'
                    break
            if ARGS.t == 'm' and feature.find('INSDFeature_key').text ==\
                    'mat_peptide':
                mat_peptide = Source(feature)
                if 'product' not in mat_peptide.qualifiers.keys() or\
                        mat_peptide.qualifiers['product'].casefold()\
                        not in FEATURE_NAMES:
                    continue
                mp_name = mat_peptide.qualifiers['product']
                if ARGS.s == 'a':
                    seq = mat_peptide.qualifiers['peptide']
                else:
                    for interval in mat_peptide.feature_intervals:
                        int_from.add(int(interval['interval_from']))
                        int_to.add(int(interval['interval_to']))
                    min_int = min(int_from)
                    max_int = max(int_to)
                    seq_found += 1
                    seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                    seq_name = '>' + acc + ';' + organism + ';peptide=' +\
                               mp_name + '\n'
            if ARGS.t == '5' and feature.find('INSDFeature_key').text ==\
                    '5\'UTR':
                fivep_utr = FivepUTR(feature)
                for interval in fivep_utr.feature_intervals:
                    int_from.add(int(interval['interval_from']))
                    int_to.add(int(interval['interval_to']))
                min_int = min(int_from)
                max_int = max(int_to)
                seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';5\'UTR\n'
            if ARGS.t == '3' and feature.find('INSDFeature_key').text ==\
                    '3\'UTR':
                threep_utr = FivepUTR(feature)
                print(threep_utr.feature_intervals)
                for interval in threep_utr.feature_intervals:
                    int_from.add(int(interval['interval_from']))
                    int_to.add(int(interval['interval_to']))
                min_int = min(int_from)
                max_int = max(int_to)
                seq = seqfeatures.seq_sequence[min_int - 1: max_int]
                seq_name = '>' + acc + ';' + organism + ';3\'UTR\n'
        if seq != "":
            outfile.write(seq_name + seq + '\n')
    outfile.close()
    return seq_count, seq_found


def extract_info(xmlroot, inff):
    """Extract information about all feature names"""
    genes = dict()
    cdss = dict()
    mps = dict()
    for sequence in xmlroot.iter('INSDSeq'):
        for feature in sequence.findall(
                'INSDSeq_feature-table/INSDFeature'):
            if feature.find('INSDFeature_key').text == 'gene':
                gene = Source(feature)
                if 'gene' not in gene.qualifiers.keys():
                    continue
                if gene.qualifiers['gene'].casefold() in genes.keys():
                    genes[gene.qualifiers['gene'].casefold()] += 1
                else:
                    genes[gene.qualifiers['gene'].casefold()] = 1
            elif feature.find('INSDFeature_key').text == 'CDS':
                cds = CDS(feature)
                if cds.qualifiers['product'].casefold() in cdss.keys():
                    cdss[cds.qualifiers['product'].casefold()] += 1
                else:
                    cdss[cds.qualifiers['product'].casefold()] = 1
            elif feature.find('INSDFeature_key').text == 'mat_peptide':
                mat_pep = Source(feature)
                if mat_pep.qualifiers['product'].casefold() in mps.keys():
                    mps[mat_pep.qualifiers['product'].casefold()] += 1
                else:
                    mps[mat_pep.qualifiers['product'].casefold()] = 1
    inff.write('GENES:\n{}\n\nCDS:\n{}\n\nMATURE PEPTIDES:\n{}'.format(
        {k: v for k, v in sorted(genes.items(), key=lambda item: item[1],
                                 reverse=True)},
        {k: v for k, v in sorted(cdss.items(), key=lambda item: item[1],
                                 reverse=True)},
        {k: v for k, v in sorted(mps.items(), key=lambda item: item[1],
                                 reverse=True)}))


if __name__ == "__main__":
    import sys
    import argparse
    import xml.etree.ElementTree as Et

    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input xml-file', required=True)
    PARSER.add_argument('-o', type=str, help='output file', required=True)
    PARSER.add_argument('-t', choices=['s', 'g', 'c', 'm', '5', '3'], type=str,
                        help='type of feature (source[s] = full nt-sequence;'
                             'gene[g]; CDS[c]; mature peptide[m]; 5\'UTR[5];'
                             '3\'UTR[3])', default='s')
    PARSER.add_argument('-n', type=str, help='Feature name', nargs='+')
    PARSER.add_argument('-s', choices=['a', 'n'], type=str,
                        help='Sequence type', default='n')
    PARSER.add_argument('-l', action="store_true", default=False,
                        help='Switch for output of info-file')
    PARSER.add_argument('-a', type=str, help='input feature alias file',
                        required=False)
    ARGS = PARSER.parse_args()
    if ARGS.t in ['g', 'c', 'm'] and not ARGS.n:
        sys.exit('No feature name. Use the -n option to give the name of the '
                 'selected feature. Note: this can be a list of names separated'
                 'byspaces')
    if ARGS.t in ['g', '5', '3'] and ARGS.s == 'a':
        sys.exit('Genes and UTRs can only have nt-sequence type. Leave out the'
                 '-s option or set -s n')
    if ARGS.a and not ARGS.n:
        sys.exit('A feature name (-n) is required to use an feature alias-file'
                 ' (-a)')
    if ARGS.n:
        FEATURE_NAMES = [item.strip().casefold() for item in ARGS.n]
        if ARGS.a:
            # creates a union set of name(s) given at -n and the alias-file row
            # containing either of theses name(s) and store in 'FEATURE_NAMES'
            ALIASES = list()
            with open(ARGS.a) as f:
                for line in f.readlines():
                    alias_lst = line.split(',')
                    alias_lst = [x.strip().casefold() for x in alias_lst]
                    ALIASES.append(alias_lst)
            for alias in ALIASES:
                if set(FEATURE_NAMES).isdisjoint(set(alias)):
                    continue
                FEATURE_NAMES = list(set(alias).union(set(FEATURE_NAMES)))
                break
    TREE = Et.parse(ARGS.f)
    ROOT = TREE.getroot()
    if ARGS.l:
        INFO_FILE = open(ARGS.f.split('.')[0] + '.info', 'w')
        extract_info(ROOT, INFO_FILE)
    SEQ_COUNT, SEQ_FOUND = extract_seq()
    if ARGS.t == 's':
        LBL = 'source sequences'
    elif ARGS.t == 'g':
        LBL = 'genes'
    elif ARGS.t == 'c':
        LBL = 'coding regions (CDS)'
    elif ARGS.t == '5':
        LBL = '5\'UTR regions'
    elif ARGS.t == 'm':
        LBL = 'mature peptides'
    else:
        LBL = '3\'UTR regions'
    print('{} {} found from {} sequences ({} %)'.
          format(SEQ_FOUND, LBL, SEQ_COUNT, round(100 * SEQ_FOUND/SEQ_COUNT)))
