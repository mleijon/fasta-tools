#!/usr/bin/python

import xml.etree.ElementTree as ET

tree = ET.parse('sequence.gbc_ratCoV.xml')
root = tree.getroot()
print(root.iterfind('INSDQualifier_value'))