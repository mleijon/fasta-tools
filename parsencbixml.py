#!/usr/bin/python


class SeqFeatures:
    def __init__(self, seq):
        self.seq_locus = seq.find('INSDSeq_locus').text
        self.seq_length = seq.find('INSDSeq_length').text
        self.seq_strandedness = seq.find('INSDSeq_strandedness').text
        self.seq_moltype = seq.find('INSDSeq_moltype').text
        self.seq_topology = seq.find('INSDSeq_topology').text
        self.seq_division = seq.find('INSDSeq_division').text
        self.seq_updateDate = seq.find('INSDSeq_update-date').text
        self.seq_createDate = seq.find('INSDSeq_create-date').text
        self.seq_definition = seq.find('INSDSeq_definition').text
        self.seq_primaryAccession = seq.find('INSDSeq_primary-accession').text
        self.seq_accessionVersion = seq.find('INSDSeq_accession-version').text
        self.seq_OtherSeqIds = list()
        for item in seq.findall('INSDSeq_other-seqids/INSDSeqid'):
            self.seq_OtherSeqIds.append(item.text)
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
        self.seq_comment = seq.find('INSDSeq_comment').text
        self.seq_sequence = seq.find('INSDSeq_sequence').text


class FivepUTR:
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


class Source (FivepUTR):
    def __init__(self, ft):
        super().__init__(ft)
        qualifiers = dict()
        qualifierlst = list()
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
        qualifierlst.append(qualifiers)
        self.qualifiers = qualifierlst


class Gene (Source):
    def __init__(self, ft):
        super().__init__(ft)


class CDS(Source):
    def __init__(self, ft):
        super().__init__(ft)
        if not ft.find('INSDFeature_operator') is None:
            self.feature_operator = ft.find('INSDFeature_operator').text


class MatPeptide(CDS):
    def __init__(self, ft):
        super().__init__(ft)


class ThreepUTR(FivepUTR):
    def __init__(self, ft):
        super().__init__(ft)


class MiscFeature(Source):
    def __init__(self, ft):
        super().__init__(ft)


class MiscRNA(Source):
    def __init__(self, ft):
        super().__init__(ft)


class StemLoop(FivepUTR):
    def __init__(self, ft):
        super().__init__(ft)


class Regulatory(Source):
    def __init__(self, ft):
        super().__init__(ft)


def extract_seq(xmlroot):
    outfile = open(ARGS.o, 'w')
    for sequence in xmlroot.iter('INSDSeq'):
        seqfeatures = SeqFeatures(sequence)
        if ARGS.t == 's':
            seq = seqfeatures.seq_sequence
            acc = seqfeatures.seq_primaryAccession
            organism = seqfeatures.seq_organism
            seq_name = '>' + acc + ';' + organism + '\n'
            outfile.write(seq_name + seq + '\n')
        if ARGS.t == 'g':
            pass
        if ARGS.t == 'c':
            pass
        if ARGS.t == 'm':
            pass
        if ARGS.t == '5':
            pass
        if ARGS.t == '3':
            pass
    outfile.close()


if __name__ == "__main__":
    import argparse
    import xml.etree.ElementTree as ET

    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='input xml-file', required=True)
    PARSER.add_argument('-o', type=str, help='output file', required=True)
    PARSER.add_argument('-t', choices=['s', 'g', 'c', 'm', '5', '3'], type=str,
                        help='type of feature (source[s] = full nt-sequence;'
                             'gene[g]; CDS[c]; mature peptide[m]; 5\'UTR[5];'
                             '3\'UTR[3])', default='s')
    PARSER.add_argument('-n', type=str, help='Feature name', default=None)
    PARSER.add_argument('-s', choices=['s', 'n'], type=str,
                        help='Sequence type', default='n')
    ARGS = PARSER.parse_args()
    tree = ET.parse(ARGS.f)
    root = tree.getroot()
    extract_seq(root)
    # for sequence in root.iter('INSDSeq'):
    #
    #
    #   test = SeqFeatures(sequence)
    #     for i, feature in enumerate(
    #             sequence.findall('INSDSeq_feature-table/INSDFeature')):
    #         #print('Feature: {}'.format(i + 1))
    #         if feature.find('INSDFeature_key').text == 'source':
    #             test2 = Source(feature)
    #             attrs = vars(test2)
    #             #print(attrs)
    #         elif feature.find('INSDFeature_key').text == '5\'UTR':
    #             test2 = FivepUTR(feature)
    #             attrs = vars(test2)
    #             #print(attrs)
    #         elif feature.find('INSDFeature_key').text == 'gene':
    #             test2 = Gene(feature)
    #             attrs = vars(test2)
    #             #print(attrs)
    #         elif feature.find('INSDFeature_key').text == 'CDS':
    #             test2 = CDS(feature)
    #             attrs = vars(test2)
    #             #print(attrs)
    #         elif feature.find('INSDFeature_key').text == 'mat_peptide':
    #             test2 = MatPeptide(feature)
    #             attrs = vars(test2)
    #             #print(attrs)
    #         elif feature.find('INSDFeature_key').text == '3\'UTR':
    #             test2 = ThreepUTR(feature)
    #             attrs = vars(test2)
    #             #print(attrs)
    #         elif feature.find('INSDFeature_key').text == 'misc_feature':
    #             test2 = MiscFeature(feature)
    #             attrs = vars(test2)
    #             #print(attrs)
    #         elif feature.find('INSDFeature_key').text == 'misc_RNA':
    #             test2 = MiscRNA(feature)
    #             attrs = vars(test2)
    #             print(attrs)
    #         elif feature.find('INSDFeature_key').text == 'stem_loop':
    #             test2 = StemLoop(feature)
    #             attrs = vars(test2)
    #             print(attrs)
    #         elif feature.find('INSDFeature_key').text == 'regulatory':
    #             test2 = Regulatory(feature)
    #             attrs = vars(test2)
    #             print(attrs)