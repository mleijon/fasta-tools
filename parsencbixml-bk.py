#!/usr/bin/python

import xml.etree.ElementTree as ET
all_features = [[()]]
tree = ET.parse('sequence.gbc_1.xml')
root = tree.getroot()
path = 'INSDFeature_quals/INSDQualifier/'
for features in root.iter('INSDSeq_feature-table'):
    featurenr = 0
    for feature in features.iter('INSDFeature'):
        featurenr += 1
        print('feature nr: {}'.format(featurenr))
        for qualifier in feature.iter('INSDQualifier'):
            if not qualifier.find('INSDQualifier_name') is None:
                name = qualifier.find('INSDQualifier_name').text
            else:
                name = ""
            if not qualifier.find('INSDQualifier_value') is None:
                value = qualifier.find('INSDQualifier_value').text
            else:
                value = ""
            print('{}: {}'.format(name, value))



        # name = feature.find(path + 'INSDQualifier_name')
        #
        # for name, value in zip(feature.find(path + 'INSDQualifier_name'),
        #                        feature.find(path + 'INSDQualifier_value')):
        #     print('{}: {}'.format(name.text, value.text))
