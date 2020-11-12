#!/usr/bin/python

import xml.etree.ElementTree as ET

tree = ET.parse('sequence.gbc_ratCoV.xml')
root = tree.getroot()
featurenr = 0
for child in root.iter('INSDFeature'):
    featurenr += 1
    print('feature nr: {}'.format(featurenr))
    for item in child[1].iter():
        print('{}: {}'.format(item.tag, item.text))
