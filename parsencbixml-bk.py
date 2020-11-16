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
        self.seq_sequence = seq.find('INSDSeq_sequence').text


class FivepUTR:
    def __init__(self, ft):
        self.feature_key = ft.find('INSDFeature_key').text
        self.feature_location = ft.find('INSDFeature_location').text
        interval = dict()
        intervallst = list()
        for item in ft.findall('INSDFeature_intervals/INSDInterval'):
            interval['interval_from'] = item.find('INSDInterval_from').text
            interval['interval_to'] = item.find('INSDInterval_to').text
            interval['interval_acc'] = item.find('INSDInterval_accession').text
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


class Gene:
    pass


class CDs:
    pass


class MatPeptide:
    pass


class ThreepUTR:
    pass


if __name__ == "__main__":
    import argparse
    import xml.etree.ElementTree as ET

    PARSER = argparse.ArgumentParser(description='TBD')
    PARSER.add_argument('-f', type=str, help='xml-filename', required=True)
    ARGS = PARSER.parse_args()
    tree = ET.parse(ARGS.f)
    root = tree.getroot()
    for sequence in root.iter('INSDSeq'):
        test = SeqFeatures(sequence)
        for feature in sequence.findall('INSDSeq_feature-table/INSDFeature'):
            if feature.find('INSDFeature_key').text == 'source':
                test2 = Source(feature)
                attrs = vars(test2)
                print(attrs)
    # path = 'INSDFeature_quals/INSDQualifier/'
    # for features in root.iter('INSDSeq_feature-table'):
    #     featurenr = 0
    #     for feature in features.iter('INSDFeature'):
    #         featurenr += 1
    #         print('feature nr: {}'.format(featurenr))
    #         for qualifier in feature.iter('INSDQualifier'):
    #             if not qualifier.find('INSDQualifier_name') is None:
    #                 name = qualifier.find('INSDQualifier_name').text
    #             else:
    #                 name = ""
    #             if not qualifier.find('INSDQualifier_value') is None:
    #                 value = qualifier.find('INSDQualifier_value').text
    #             else:
    #                 value = ""
    #             print('{}: {}'.format(name, value))
