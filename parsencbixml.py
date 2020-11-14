#!/usr/bin/python

import xml.etree.ElementTree as ET

tree = ET.parse('sequence.gbc_1.xml')
root = tree.getroot()


class SeqFeature:
    def __init__(self):
        self.seq_locus = root.find('INSDSeq/INSDSeq_locus').text
        self.seq_length = root.find('INSDSeq/INSDSeq_length').text
        self.seq_strandedness = root.find('INSDSeq/INSDSeq_strandedness').text
        self.seq_moltype = root.find('INSDSeq/INSDSeq_moltype').text
        self.seq_topology = root.find('INSDSeq/INSDSeq_topology').text
        self.seq_division = root.find('INSDSeq/INSDSeq_division').text
        self.seq_updateDate = root.find('INSDSeq/INSDSeq_update-date').text
        self.seq_createDate = root.find('INSDSeq/INSDSeq_create-date').text
        self.seq_definition = root.find('INSDSeq/INSDSeq_definition').text
        self.seq_primaryAccession = root.find('INSDSeq/'
                                              'INSDSeq_primary-accession').text
        self.seq_accessionVersion = root.find('INSDSeq/'
                                              'INSDSeq_accession-version').text
        self.seq_OtherSeqIds = list()
        for item in root.findall('INSDSeq/INSDSeq_other-seqids/INSDSeqid'):
            self.seq_OtherSeqIds.append(item.text)


if __name__ == "__main__":
    test = SeqFeature()
    attr = vars(test)
    print(attr)
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
