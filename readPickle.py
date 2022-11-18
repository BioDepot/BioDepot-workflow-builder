import base64
import sys
import warnings
import ast
from ast import literal_eval

from xml.etree.ElementTree import TreeBuilder, Element, ElementTree, parse

from collections import defaultdict, namedtuple
from itertools import chain, count

import pickle as pickle
import json
import pprint

def loads(string, format):
    if format == "literal":
        return literal_eval(string)
    elif format == "json":
        return json.loads(string)
    elif format == "pickle":
        return pickle.loads(base64.decodebytes(string.encode("ascii")))
    else:
        raise ValueError("Unknown format")

xmlFile='/home/lhhung/ssh_repos/workflows/DNA/GATK_germline_variant/GATK_germline_variant.ows'
tree = parse(xmlFile)
doc = parse(xmlFile)
scheme_el = doc.getroot()
for property_el in scheme_el.findall("node_properties/properties"):
        node_id = property_el.attrib.get("node_id")
        print(node_id)

        format = property_el.attrib.get("format", "pickle")
        print(format)

        if "data" in property_el.attrib:
            # data string is 'encoded' with 'repr' i.e. unicode and
            # nonprintable characters are \u or \x escaped.
            # Could use 'codecs' module?
            data = string_eval(property_el.attrib.get("data"))
        else:
            data = property_el.text
        print(data)
        properties = loads(data, format)
        print(properties)

